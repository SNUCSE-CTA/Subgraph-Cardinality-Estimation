#include <vector>
#include <map>
#include <random>
#include "include/daf_candidate_space.h"
#include "global/timer.h"

namespace daf {

    BipartiteMaximumMatching BPSolver, BPTSolver;
    int num_edges[MAX_QUERY_VERTEX][MAX_QUERY_VERTEX], num_cur_edges[MAX_QUERY_VERTEX][MAX_QUERY_VERTEX];;

    CandidateSpace::CandidateSpace(DataGraph *data, Option filter_option) {
        data_ = data;
        BitsetCS = new bool*[MAX_QUERY_VERTEX];
        for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
            BitsetCS[i] = new bool[data->GetNumVertices()];
        }
        BitsetEdgeCS = new bool*[MAX_QUERY_EDGE];
        for (int i = 0; i < MAX_QUERY_EDGE; i++) {
            BitsetEdgeCS[i] = new bool[data->edge_info_.size()];
            memset(BitsetEdgeCS[i], false, data->edge_info_.size());
        }
        in_neighbor_cs = new bool[data->GetNumVertices()];
        neighbor_label_frequency.resize(data->GetNumVertices());
        num_visit_cs_ = new QueryDegree[data_->GetNumVertices()];
        candidate_set_.resize(MAX_QUERY_VERTEX);
        opt = filter_option;
    }

    CandidateSpace::~CandidateSpace() {
        for (int i = 0; i < MAX_QUERY_EDGE; i++) {
            delete[] BitsetEdgeCS[i];
        }
        delete[] BitsetEdgeCS;
        delete[] num_visit_cs_;
        delete[] in_neighbor_cs;
    }

    bool CandidateSpace::BuildCS(QueryGraph *query, DAG *dag) {
        if (data_->num_four_cycles_indexed == 0) {
            four_cycle_memo_old.clear();
        }
        query_ = query;
        dag_ = dag;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            memset(BitsetCS[i], false, data_->GetNumVertices());
        }
        for (int i = 0; i < query_->edge_info_.size(); i++) {
            memset(BitsetEdgeCS[i], false, data_->edge_info_.size());
        }
        memset(num_edges, 0, sizeof(num_edges));
        memset(num_cur_edges, 0, sizeof(num_cur_edges));
        memset(num_visit_cs_, 0, data_->GetNumVertices());
        if (opt.egonet_filter >= NEIGHBOR_BIPARTITE_SAFETY) {
            BPSolver.global_initialize(query_->GetMaxDegree(), data_->GetMaxDegree());
        }
        if (opt.structure_filter >= TRIANGLE_BIPARTITE_SAFETY) {
            BPTSolver.global_initialize(query_->GetMaxDegree(), data_->max_num_trigs);
        }
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            candidate_set_[i].clear();
        }
        dag_->BuildDAG(-1);
        if (!FilterByTopDownWithInit()) {
            fprintf(stderr, "INITIALIZATION FAILED\n");
            exit(4);
        }
        Timer csfilter_timer; csfilter_timer.Start();
        if (opt.refinement_order == PRIORITY_FIRST) {
            Filter(true);
        }
        else {
            Filter(false);
        }
        csfilter_timer.Stop();
        fprintf(stdout, "CS prune time : %.02lf ms\n",csfilter_timer.GetTime());
        ConstructCS();
        return true;
    }


    bool CandidateSpace::FilterByTopDownWithInit() {
        bool result = true;
        result = InitRootCandidates();
        if (!result) return false;

        for (Size i = 1; i < query_->GetNumVertices(); ++i) {
            Vertex cur = dag_->GetVertexOrderedByBFS(i);
            Label cur_label = query_->GetLabel(cur);
            QueryDegree num_parent = 0, needed = dag_->GetNumParents(cur);

            for (Size p = 0; p < dag_->GetNumParents(cur); p++) {
                Vertex parent = dag_->GetParent(cur, p);
                int query_edge_idx = query_->GetEdgeIndex(parent, cur);
                for (Vertex parent_cand : candidate_set_[parent]) {
                    for (int data_edge_idx : data_->GetIncidentEdges(parent_cand, cur_label)) {
                        Vertex cand = data_->opposite(data_edge_idx, parent_cand);
                        if (num_visit_cs_[cand] < num_parent) continue;
                        if (data_->GetDegree(cand) < query_->GetDegree(cur)) break;
                        if (data_->GetCoreNum(cand) < query_->GetCoreNum(cur)) continue;
                        if (data_->GetELabel(data_edge_idx) != query_->GetELabel(query_edge_idx)) continue;
                        for (int l = 0; l < data_->GetNumLabels(); l++) {
                            if (data_->incident_edges_[cand][l].size() < query_->incident_edges_[cur][l].size()) {
                                goto nxt_candidate;
                            }
                        }
                        if (opt.structure_filter >= TRIANGLE_SAFETY and
                            data_->GetNumLocalTriangles(data_edge_idx) < query_->GetNumLocalTriangles(query_edge_idx)) continue;
                        if (opt.structure_filter >= FOURCYCLE_SAFETY and data_->num_four_cycles_indexed > 0 and
                            data_->GetNumLocalFourCycles(data_edge_idx) < query_->GetNumLocalFourCycles(query_edge_idx)) continue;
                        if (num_visit_cs_[cand] == num_parent) {
                            num_visit_cs_[cand] += 1;
                            if (num_visit_cs_[cand] == 1) {
                                candidate_set_[cur].emplace_back(cand);
                                BitsetCS[cur][cand] = true;
                            }
                        }
                        nxt_candidate:
                        continue;
                    }
                }
                num_parent += 1;
            }
//            fprintf(stdout, "Candidate for %d : %d elements\n", cur, candidate_set_[cur].size());
            for (Size j = 0; j < candidate_set_[cur].size(); ++j) {
                Vertex cand = candidate_set_[cur][j];
//                printf("  %d...",cand);
                BitsetCS[cur][cand] = false;
                if (num_visit_cs_[cand] == num_parent) {
//                    printf("ok\n");
                    BitsetCS[cur][cand] = true;
                }
                else {
//                    printf("has %d parents only (%d needed)\n", num_visit_cs_[cand], num_parent);
                    candidate_set_[cur][j] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    j--;
                }
                num_visit_cs_[cand] = 0;
            }
//            fprintf(stdout, "====> %d elements left\n",candidate_set_[cur].size());
            if (candidate_set_[cur].empty()) {
                result = false;
                exit(2);
            }
        }

        int cs_edge = 0, cs_vertex = 0;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            for (int q_edge_idx : query_->all_incident_edges_[i]) {
                int q_nxt = query_->to_[q_edge_idx];
                for (int j : candidate_set_[i]) {
                    for (int d_edge_idx : data_->GetIncidentEdges(j, query_->GetLabel(q_nxt))) {
                        if (opt.structure_filter >= TRIANGLE_SAFETY and
                            data_->GetNumLocalTriangles(d_edge_idx) < query_->GetNumLocalTriangles(q_edge_idx)) continue;
                        if (opt.structure_filter >= FOURCYCLE_SAFETY and data_->num_four_cycles_indexed > 0 and
                            data_->GetNumLocalFourCycles(d_edge_idx) < query_->GetNumLocalFourCycles(q_edge_idx)) continue;
                        int d_nxt = data_->to_[d_edge_idx];
                        if (BitsetCS[q_nxt][d_nxt]) {
                            BitsetEdgeCS[q_edge_idx][d_edge_idx] = true;
                            num_edges[i][q_nxt]++;
                            cs_edge++;
                        }
                    }
                }
            }
            cs_vertex += candidate_set_[i].size();
        }
//        fprintf(stderr, "Initial CS : %d vertex, %d edges\n",cs_vertex, cs_edge);
//        printCS();
        return result;
    }

    bool CandidateSpace::BipartiteSafety(Vertex cur, Vertex cand) {
        if (query_->adj_list[cur].size() == 1) {
            int uc = query_->adj_list[cur][0];
            int query_edge_index = query_->GetEdgeIndex(cur, uc);
            int label = query_->GetLabel(uc);
            return std::any_of(data_->GetIncidentEdges(cand, label).begin(),
                               data_->GetIncidentEdges(cand, label).end(),
                               [&](int data_edge_index) {
                                   return BitsetEdgeCS[query_edge_index][data_edge_index];
                               });
        }
        BPSolver.reset();
        int i = 0, j = 0;
        for (int query_edge_index : query_->all_incident_edges_[cur]) {
            j = 0;
            for (int edge_id : data_->all_incident_edges_[cand]) {
                if (BitsetEdgeCS[query_edge_index][edge_id]) {
                    BPSolver.add_edge(i, j);
                }
                j++;
            }
            i++;
        }
        bool ok = (BPSolver.solve() == (query_->GetDegree(cur)));
        return ok;
    }
    /* Legacy code ---- */
    Size CandidateSpace::GetDAGNextCount(Vertex cur, bool topdown) {
        return topdown? dag_->GetNumParents(cur) : dag_->GetNumChildren(cur);
    }
    Vertex CandidateSpace::GetDAGNextVertex(Vertex cur, Size idx, bool topdown) {
        return topdown? dag_->GetParent(cur, idx) : dag_->GetChild(cur, idx);
    }
    /* ---- Legacy code */

    inline bool CandidateSpace::EdgeCandidacy(int query_edge_id, int data_edge_id) {
        if (query_edge_id == -1 || data_edge_id == -1) {
            return false;
        }
        return BitsetEdgeCS[query_edge_id][data_edge_id];
    }
    bool CandidateSpace::TriangleSafety(int query_edge_id, int data_edge_id) {
        if (query_->local_triangles[query_edge_id].empty()) return true;
        auto &candidate_triangles = data_->local_triangles[data_edge_id];
        if (query_->local_triangles[query_edge_id].size() > candidate_triangles.size()) return false;
        for (auto qtv: query_->local_triangles[query_edge_id]) {
            bool found = std::any_of(candidate_triangles.begin(),
                        candidate_triangles.end(),
                        [&](auto tv) {
                            return BitsetEdgeCS[qtv.fst_edge][tv.fst_edge] and  BitsetEdgeCS[qtv.snd_edge][tv.snd_edge];
                        });
            if (!found) return false;
        }
        return true;
    }
    bool CandidateSpace::TriangleBipartiteSafety(int query_edge_id, int data_edge_id) {
        if (query_->local_triangles[query_edge_id].empty()) return true;
        auto &candidate_triangles = data_->local_triangles[data_edge_id];
        if (query_->local_triangles[query_edge_id].size() > candidate_triangles.size()) return false;
        if (query_->local_triangles[query_edge_id].size() == 1) {
            auto qtv = query_->local_triangles[query_edge_id].back();
            return std::any_of(candidate_triangles.begin(),
                               candidate_triangles.end(),
                               [&](auto tv) {
                                   return BitsetEdgeCS[qtv.fst_edge][tv.fst_edge] and  BitsetEdgeCS[qtv.snd_edge][tv.snd_edge];
                               });
        }
        BPTSolver.reset();
        int triangle_index = -1;
        for (auto qtv: query_->local_triangles[query_edge_id]) {
            triangle_index++;
            int tv_idx = -1;
            bool found = false;
            for (auto tv: candidate_triangles) {
                tv_idx++;
                if (!BitsetEdgeCS[qtv.fst_edge][tv.fst_edge]) continue;
                if (!BitsetEdgeCS[qtv.snd_edge][tv.snd_edge]) continue;
                found = true;
                BPTSolver.add_edge(triangle_index, tv_idx);
            }
            if (!found) {
                return false;
            }
        }
        bool ok = (BPTSolver.solve() == query_->local_triangles[query_edge_id].size());
        return ok;
    }
    bool CandidateSpace::FourCycleSafetyOnline(int query_edge_id, int data_edge_id) {
        int cand, nxt_cand; std::tie(cand, nxt_cand) = data_->edge_info_[data_edge_id].vp;
        int cur, nxt; std::tie(cur, nxt) = query_->edge_info_[query_edge_id].vp;
        std::pair<int, int> edgePair = {query_edge_id, data_edge_id};
        if (four_cycle_memo_old.find(edgePair) == four_cycle_memo_old.end()) {
            four_cycle_memo_old[edgePair] = std::vector<online_cycle_information>(query_->four_cycles[query_edge_id].size());
        }
        std::vector<online_cycle_information> &four_cycle_answers = four_cycle_memo_old[edgePair];
        for (int i = 0; i < query_->four_cycles[query_edge_id].size(); i++) {
            int opp_edge_idx = query_->four_cycles[query_edge_id][i].opp_edge_idx;
            auto &info = four_cycle_answers[i];
            bool validity = true;
            if (info.d_opp_edge_idx == -1) {
                validity = false;
                VertexPair opp_edge = query_->edge_info_[opp_edge_idx].vp;
                std::tie(info.third, info.fourth) = opp_edge;
                info.q_opp_edge_idx = opp_edge_idx;
                info.q_fourth_edge_idx = query_->GetEdgeIndex(cur, info.fourth);
                info.q_third_edge_idx = query_->GetEdgeIndex(nxt, info.third);
                info.third_inc_idx = 0;
                info.fourth_inc_idx = 0;
                if (query_->GetEdgeIndex(cur, info.third) != -1) counters[17]++;
                if (query_->GetEdgeIndex(nxt, info.fourth)!= -1) counters[17]++;
            }
            else {
                if (validity) validity &= EdgeCandidacy(info.q_opp_edge_idx, info.d_opp_edge_idx);
                if (validity) validity &= EdgeCandidacy(info.q_third_edge_idx, info.d_third_edge_idx);
                if (validity) validity &= EdgeCandidacy(info.q_fourth_edge_idx, info.d_fourth_edge_idx);
            }
            if (validity) {
                continue;
            }

//            steptimer_3.Start();
            Label third_label = query_->GetLabel(info.third);
            Label fourth_label = query_->GetLabel(info.fourth);
            std::vector <int> &third_cands  = data_->GetIncidentEdges(nxt_cand, third_label);
            std::vector <int> &fourth_cands = data_->GetIncidentEdges(cand, fourth_label);
            if (info.third_inc_idx >= third_cands.size()) return false;

            while (info.third_inc_idx < third_cands.size()) {
                bool cand_validity = true;
                int third_cand = data_->opposite(third_cands[info.third_inc_idx], nxt_cand);
                int fourth_cand = data_->opposite(fourth_cands[info.fourth_inc_idx], cand);
                if (cand_validity) cand_validity &= ((fourth_cand != nxt_cand) && (third_cand != cand));
                if (cand_validity) cand_validity &= EdgeCandidacy(info.q_third_edge_idx, third_cands[info.third_inc_idx]);
                if (cand_validity) cand_validity &= EdgeCandidacy(info.q_fourth_edge_idx, fourth_cands[info.fourth_inc_idx]);
                if (cand_validity) {
                    int d_opp_idx = data_->GetEdgeIndex(third_cand, fourth_cand);
                    cand_validity &= EdgeCandidacy(info.q_opp_edge_idx, d_opp_idx);
                    info.d_opp_edge_idx = d_opp_idx;
                    info.d_third_edge_idx = third_cands[info.third_inc_idx];
                    info.d_fourth_edge_idx = fourth_cands[info.fourth_inc_idx];
                }

                info.fourth_inc_idx++;
                if (info.fourth_inc_idx >= fourth_cands.size()) {
                    info.fourth_inc_idx = 0;
                    info.third_inc_idx++;
                }
                if (cand_validity) goto ok;
            }
            if (info.third_inc_idx >= third_cands.size()) {
                return false;
            }
            ok:
            continue;
        }
        return true;
    }
    bool CandidateSpace::FourCycleSafety(int query_edge_id, int data_edge_id) {
        if (data_->num_four_cycles_indexed == 0) return FourCycleSafetyOnline(query_edge_id, data_edge_id);
        if (query_->four_cycles[query_edge_id].size() > data_->four_cycles[data_edge_id].size()) return false;
        for (int i = 0; i < query_->four_cycles[query_edge_id].size(); i++) {
            auto &q_info = query_->four_cycles[query_edge_id][i];
            bool found = false;
            for (int j = 0; j < data_->four_cycles[data_edge_id].size(); j++) {
                auto &d_info =  data_->four_cycles[data_edge_id][j];
                bool validity = true;
                validity &= BitsetEdgeCS[q_info.third_edge_idx][d_info.third_edge_idx];
                if (!validity) continue;
                validity &= BitsetEdgeCS[q_info.fourth_edge_idx][d_info.fourth_edge_idx];
                if (!validity) continue;
                validity &= BitsetEdgeCS[q_info.opp_edge_idx][d_info.opp_edge_idx];
                if (!validity) continue;
                if (validity and q_info.one_three_idx != -1) validity &= EdgeCandidacy(q_info.one_three_idx, d_info.one_three_idx);
                if (validity and q_info.two_four_idx != -1) validity &= EdgeCandidacy(q_info.two_four_idx, d_info.two_four_idx);
                if (validity) {
                    found = true;
                    goto nxt_cycle;
                }
            }
            if (!found) return false;
            nxt_cycle:
            continue;
        }
//        bool ok = (BPQsolver.solve() == (query_->four_cycles[query_edge_id].size()));
        return true;
    }
    bool CandidateSpace::StructureSafety(int query_edge_id, int data_edge_id) {
        switch (opt.structure_filter) {
            case NO_STRUCTURE_FILTER:
                return true;
            case TRIANGLE_SAFETY:
                return TriangleSafety(query_edge_id, data_edge_id);
            case TRIANGLE_BIPARTITE_SAFETY:
                return TriangleBipartiteSafety(query_edge_id, data_edge_id);
            case FOURCYCLE_SAFETY:
                return TriangleSafety(query_edge_id, data_edge_id) and FourCycleSafety(query_edge_id, data_edge_id);
        }
        return true;
    }
    bool CandidateSpace::StructureFilter(int cur, int cand, int direction) {
        std::vector<int>* query_neighbors = &(query_->all_incident_edges_[cur]);
        if (direction == -1) {
            // for bottom-up
            query_neighbors = &(dag_->dag_child_edges_[cur]);
        }
        else if (direction == 1) {
            // for top-down
            query_neighbors = &(dag_->dag_parent_edges_[cur]);
        }
        for (int query_edge_idx : *query_neighbors) {
            int nxt = query_->to_[query_edge_idx];
            int nxt_label = query_->GetLabel(nxt);
            bool found = false;
            for (int data_edge_idx : data_->GetIncidentEdges(cand, nxt_label)) {
                Vertex nxt_cand = data_->to_[data_edge_idx];
                if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) continue;
                if (!StructureSafety(query_edge_idx, data_edge_idx)) {
                    num_cur_edges[cur][nxt]--;
                    num_cur_edges[nxt][cur]--;
                    BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                    BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]] = false;
                    continue;
                }
                found = true;
            }
            if (!found) {
                return false;
            }
        }
        return true;
    }


    void CandidateSpace::PrepareNeighborSafety(Vertex cur) {
        for (Vertex q_neighbor : query_->adj_list[cur]) {
            neighbor_label_frequency[query_->GetLabel(q_neighbor)]++;
            for (Vertex d_neighbor : candidate_set_[q_neighbor]) {
                in_neighbor_cs[d_neighbor] = true;
            }
        }
    }

    bool CandidateSpace::CheckNeighborSafety(Vertex cur, Vertex cand) {
        for (Vertex d_neighbor : data_->adj_list[cand]) {
            if (in_neighbor_cs[d_neighbor]) {
                neighbor_label_frequency[data_->GetLabel(d_neighbor)]--;
            }
        }
        bool valid = true;
        for (int l = 0; l < data_->GetNumLabels(); ++l) {
            if (neighbor_label_frequency[l] > 0) {
                valid = false;
                break;
            }
        }
        for (Vertex d_neighbor : data_->adj_list[cand]) {
            if (in_neighbor_cs[d_neighbor]) {
                neighbor_label_frequency[data_->GetLabel(d_neighbor)]++;
            }
        }
        return valid;
    }

    bool CandidateSpace::EdgeBipartiteSafety(Vertex cur, Vertex cand) {
        if (query_->all_incident_edges_[cur].size() == 1) {
            int q_edge_id = query_->all_incident_edges_[cur][0];
            for (int edge_id : data_->all_incident_edges_[cand]) {
                if (BitsetEdgeCS[q_edge_id][edge_id])
                    return true;
            }
            return false;
        }
        std::vector<std::pair<int, int>> edge_pairs;
        int ii = 0, jj = 0;
        BPSolver.reset();
        for (int query_edge_index : query_->all_incident_edges_[cur]) {
            Vertex uc = query_->opposite(query_edge_index, cur);
            jj = 0;
            for (int edge_id : data_->all_incident_edges_[cand]) {
                Vertex vc = data_->opposite(edge_id, cand);
                if (data_->GetDegree(vc) < query_->GetDegree(uc)) break;
                if (BitsetEdgeCS[query_edge_index][edge_id]) {
                    BPSolver.add_edge(ii, jj);
                    edge_pairs.emplace_back(ii, jj);
                }
                jj++;
            }
            ii++;
        }
        bool b = BPSolver.FindUnmatchableEdges(query_->all_incident_edges_[cur].size());
        if (!b) return false;
        for (auto &[i, j] : edge_pairs) {
            if (!BPSolver.matchable[i][j]) {
//                printf("Edge (%d, %d) is unmatchable\n", i, j);
                int left_unmatch = query_->all_incident_edges_[cur][i];
                int right_unmatch = data_->all_incident_edges_[cand][j];
                BitsetEdgeCS[left_unmatch][right_unmatch] = false;
                BitsetEdgeCS[query_->opposite_edge[left_unmatch]][data_->opposite_edge[right_unmatch]] = false;
            }
        }
//        printf("\n");
        return true;
    }


//    bool CandidateSpace::EdgeBipartiteSafety(Vertex cur, Vertex cand) {
//        int ii = 0, jj = 0;
//        BPSolver.reset();
//        for (int query_edge_index : query_->all_incident_edges_[cur]) {
//            Vertex uc = query_->opposite(query_edge_index, cur);
//            jj = 0;
//            for (int edge_id : data_->all_incident_edges_[cand]) {
//                Vertex vc = data_->opposite(edge_id, cand);
//                if (data_->GetDegree(vc) < query_->GetDegree(uc)) break;
//                if (BitsetEdgeCS[query_edge_index][edge_id])
//                    BPSolver.add_edge(ii, jj);
//                jj++;
//            }
//            ii++;
//        }
////        if (BPSolver.solve(ii) < query_->adj_list[cur].size()) return false;
////        BPSolver.print();
//        ii = 0;
//        for (int query_edge_idx : query_->all_incident_edges_[cur]) {
//            int nxt = query_->to_[query_edge_idx];
//            bool found = false;
//            BPSolver.reset(false);
//            if (BPSolver.solve(ii) < query_->adj_list[cur].size() - 1) {
//                return false;
//            }
//            jj = -1;
//            for (int data_edge_idx : data_->all_incident_edges_[cand]) {
//                Vertex nxt_cand = data_->to_[data_edge_idx];
//                jj++;
//                if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
//                if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
//                    continue;
//                }
//                if (BPSolver.right[jj] == -1) {
//                    found = true;
//                    continue;
//                }
//                std::memset(BPSolver.used, false, sizeof(bool) * query_->adj_list[cur].size());
//                BPSolver.used[ii] = true;
//                if (BPSolver.single_dfs(BPSolver.right[jj])) {
//                    found = true;
//                }
//                else {
////                    printf("Edge (%d, %d) is unmatchable\n", ii, jj);
//                    num_cur_edges[cur][nxt]--;
//                    num_cur_edges[nxt][cur]--;
//                    BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
//                    BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]] = false;
//                    if(!BPSolver.remove_edge(ii, jj)){
//                        return false;
//                    }
//                }
//            }
//            if (!found) {
//                return false;
//            }
//            ii++;
//        }
//        return true;
//    }
    bool CandidateSpace::EgonetFilter(int cur, int cand) {
        switch (opt.egonet_filter) {
            case NEIGHBOR_SAFETY:
                return CheckNeighborSafety(cur, cand);
            case NEIGHBOR_BIPARTITE_SAFETY:
                return BipartiteSafety(cur, cand);
            case EDGE_BIPARTITE_SAFETY:
                return EdgeBipartiteSafety(cur, cand);
        }
    }
    bool CandidateSpace::Filter(bool priority_first) {
        std::vector<int> local_stage(query_->GetNumVertices(), 0);
        std::vector<double> priority(query_->GetNumVertices(), 0);
        std::priority_queue<variable> candidate_queue;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            priority[i] = 0.5 + 0.05 * (i == dag_->GetRoot()) - 0.05 * (query_->adj_list[i].size() == 1);
//            priority[i] = 0.5;
            candidate_queue.push({priority[i], 0, i});
        }
        int queue_pop_count = 0;
        int maximum_queue_cnt = 5 * query_->GetNumVertices();
        int current_stage = 0;
        if (!priority_first) goto dag_dp;
        fprintf(stderr, "Priority-First Cand Filtering\n");
        while (!candidate_queue.empty()) {
            memcpy(num_cur_edges, num_edges, sizeof(num_edges));
            auto [pri, stage, cur] = candidate_queue.top();
            candidate_queue.pop();
            if (stage < local_stage[cur]) continue;
            current_stage++;
            queue_pop_count++;
            int bef_cand_size = candidate_set_[cur].size();

            if (opt.egonet_filter == NEIGHBOR_SAFETY) {
                std::fill(neighbor_label_frequency.begin(), neighbor_label_frequency.end(), 0);
                memset(in_neighbor_cs, false, data_->GetNumVertices());
                PrepareNeighborSafety(cur);
            }

            for (int i = 0; i < candidate_set_[cur].size(); i++) {
                Vertex cand = candidate_set_[cur][i];
                bool valid = true;
                if (opt.structure_filter > NO_STRUCTURE_FILTER) {
                    valid = StructureFilter(cur, cand, 0);
                }
                if (valid) valid = EgonetFilter(cur, cand);
                if (!valid) {
                    int removed = candidate_set_[cur][i];
                    for (int query_edge_idx : query_->all_incident_edges_[cur]) {
                        int nxt = query_->to_[query_edge_idx];
                        for (int data_edge_idx : data_->GetIncidentEdges(removed, query_->GetLabel(nxt))) {
                            int nxt_cand = data_->to_[data_edge_idx];
                            if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                            if (BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
                                num_cur_edges[cur][nxt]--;
                                num_cur_edges[nxt][cur]--;
                                BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                                BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]] = false;
                            }
                        }
                    }
                    candidate_set_[cur][i] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    --i;
                    BitsetCS[cur][cand] = false;
                }
            }

            if (candidate_set_[cur].empty()) {
                exit(2);
            }

            int aft_cand_size = candidate_set_[cur].size();
            if (aft_cand_size == bef_cand_size) {
                priority[cur] = 0;
                continue;
            }
            double out_prob = 1 - aft_cand_size * 1.0 / bef_cand_size;
            priority[cur] = 0;
//            fprintf(stderr, "Reduced CS[%d] from %d to %d\n",cur,bef_cand_size,aft_cand_size);
            for (Vertex nxt : query_->adj_list[cur]) {
                double removed_edge_ratio = 1 - (num_cur_edges[cur][nxt] * 1.0 / num_edges[cur][nxt]);
//                priority[nxt] = 1 - (1 - removed_edge_ratio) * (1 - priority[nxt]);
                priority[nxt] = 1 - (1 - out_prob) * (1 - priority[nxt]);
                if (priority[nxt] < opt.cutoff) continue;
                local_stage[nxt] = current_stage;
//                fprintf(stderr, "    Cur=[%d] pushes Nxt=[%d] with priority %lf\n",cur,nxt,priority[nxt]);
                candidate_queue.push({priority[nxt], current_stage, nxt});
            }
        }
        return true;
        dag_dp:
        fprintf(stderr, "DAG-DP Cand Filtering\n");
        for (int iteration = 0; iteration < 3; iteration++) {
            bool topdown = iteration % 2 == 0;
            for (Size id = 0; id < query_->GetNumVertices(); id++) {
                Size query_idx = topdown ? id : query_->GetNumVertices() - id - 1;
                Vertex cur = dag_->GetVertexOrderedByBFS(query_idx);
                if (GetDAGNextCount(cur, topdown) == 0) continue;
                for (int i = 0; i < candidate_set_[cur].size(); i++) {
                    Vertex cand = candidate_set_[cur][i];
                    bool valid = true;
                    if (opt.structure_filter > NO_STRUCTURE_FILTER) {
                        valid = StructureFilter(cur, cand, topdown ? 1 : -1);
                    }
                    if (valid) valid = EgonetFilter(cur, cand);
                    if (!valid) {
                        int removed = candidate_set_[cur][i];
                        for (int query_edge_idx : query_->all_incident_edges_[cur]) {
                            int nxt = query_->to_[query_edge_idx];
                            for (int data_edge_idx : data_->GetIncidentEdges(removed, query_->GetLabel(nxt))) {
                                int nxt_cand = data_->to_[data_edge_idx];
                                if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                                if (BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
                                    num_cur_edges[cur][nxt]--;
                                    num_cur_edges[nxt][cur]--;
                                    BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                                    BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]] = false;
                                }
                            }
                        }
                        candidate_set_[cur][i] = candidate_set_[cur].back();
                        candidate_set_[cur].pop_back();
                        --i;
                        BitsetCS[cur][cand] = false;
                    }
                }

                if (candidate_set_[cur].empty()) {
                    exit(2);
                }
            }
        }
        return true;
    }

    void CandidateSpace::ConstructCS() {
        Timer csbuild_timer; csbuild_timer.Start();
        cs_edge_.clear();
        cs_edge_.resize(query_->GetNumVertices());
        std::vector <int> CandidateIndex(data_->GetNumVertices());
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            cs_edge_[i].resize(GetCandidateSetSize(i));
        }
        Size CandidateSpaceSize = 0, CandidateEdges = 0;
        for (Vertex u = 0; u < query_->GetNumVertices(); u++) {
            Label u_label = query_->GetLabel(u);
            Size u_degree = query_->GetDegree(u);

            CandidateSpaceSize += GetCandidateSetSize(u);
            for (Size idx = 0; idx < GetCandidateSetSize(u); idx++) {
                CandidateIndex[candidate_set_[u][idx]] = idx;
                cs_edge_[u][idx].resize(query_->GetNumVertices());
            }

            for (Vertex u_adj : query_->adj_list[u]) {
                int query_edge_idx = query_->GetEdgeIndex(u_adj, u);

                for (Size v_adj_idx = 0; v_adj_idx < candidate_set_[u_adj].size(); ++v_adj_idx) {
                    Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

                    for (int data_edge_idx : data_->GetIncidentEdges(v_adj, u_label)) {
                        Vertex v = data_->opposite(data_edge_idx, v_adj);
                        if (data_->GetDegree(v) < u_degree) break;
                        if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) continue;
                        if (!BitsetCS[u][v]) continue;
                        CandidateEdges++;
//                        printf("CS_EDGE_[%d][%d][%d].push_back(%d)\n",u,CandidateIndex[v],u_adj,v_adj_idx);
                        cs_edge_[u][CandidateIndex[v]][u_adj].emplace_back(v_adj_idx);
                    }
                }
            }
        }

        CandidateEdges /= 2;
        csbuild_timer.Stop();
//        printCS();

        fprintf(stderr, "FunctionCallCounter = %d\n",functionCallCounter);
        fprintf(stdout, "CS build time : %.02lf ms\n",csbuild_timer.GetTime());
        fprintf(stderr, "CS-SZ : %u %u\n", CandidateSpaceSize, CandidateEdges);
        fprintf(stdout, "#CandidateSetSize : %u\n", CandidateSpaceSize);
        fprintf(stdout, "#CandidateSetEdges : %u\n", CandidateEdges);
        fflush(stdout);
        fflush(stderr);
    }

    bool CandidateSpace::InitRootCandidates() {
        Vertex root = dag_->GetRoot();
        Label root_label = query_->GetLabel(root);
        for (Size i = data_->GetStartOffsetByLabel(root_label);
             i < data_->GetEndOffsetByLabel(root_label); ++i) {
            Vertex cand = data_->GetVertexBySortedLabelOffset(i);
            if (data_->GetDegree(cand) < query_->GetDegree(root)) break;
            if (data_->GetCoreNum(cand) < query_->GetCoreNum(root)) continue;
            for (int l = 0; l < data_->GetNumLabels(); l++) {
                if (data_->incident_edges_[cand][l].size() < query_->incident_edges_[root][l].size()) {
                    goto nxt_candidate;
                }
            }
            candidate_set_[root].emplace_back(cand);
            BitsetCS[root][cand] = true;
            nxt_candidate:;
        }
        return !candidate_set_[root].empty();
    }

    void CandidateSpace::printCS() {
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            fprintf(stdout, "Query [%d] : \n",i);
            for (int x : candidate_set_[i]) {
                fprintf(stdout, "\t\t%d : ", x);
                for (int uc : query_->adj_list[i]) {
                    for (int y : candidate_set_[uc]) {
                        int de_id = data_->GetEdgeIndex(x, y);
                        int qe_id = query_->GetEdgeIndex(i, uc);
                        if (EdgeCandidacy(qe_id, de_id)) {
                            fprintf(stdout, "(%d, %d) ", uc, y);
                        }
                    }
                }
                printf("\n");
            }
            fprintf(stdout, "\n");
        }
    }
}