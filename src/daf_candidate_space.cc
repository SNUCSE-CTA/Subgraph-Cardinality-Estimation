#include <vector>
#include <map>
#include <random>
#include "include/daf_candidate_space.h"
#include "global/timer.h"
#define BIPARTITE_SAFETY
//#define NEIGHBOR_SAFETY
#define FOURCYCLE_SAFETY
//#define TIME_CHECK

namespace daf {
    Size global_iteration_count = 0;
    Timer steptimer_1, steptimer_2, steptimer_3, steptimer_4;
    bool is_data_sparse = true;

    BipartiteMaximumMatching BPSolver, BPTSolver;
    int temporary_left[1000];

    CandidateSpace::CandidateSpace(DataGraph *data) {
        data_ = data;
        BitsetCS = new bool*[MAX_QUERY_VERTEX];
//        BitsetCS.resize(MAX_QUERY_VERTEX);
        for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
            BitsetCS[i] = new bool[data->GetNumVertices()];
        }
        BitsetEdgeCS = new bool*[MAX_QUERY_EDGE];
        for (int i = 0; i < MAX_QUERY_EDGE; i++) {
            BitsetEdgeCS[i] = new bool[data->edge_info_.size()];
            memset(BitsetEdgeCS[i], true, data->edge_info_.size());
        }
        num_visit_cs_ = new QueryDegree[data_->GetNumVertices()];
        visited_candidates_ = new Vertex[data_->GetNumVertices() * MAX_QUERY_VERTEX];
        candidate_set_.resize(MAX_QUERY_VERTEX);
        is_data_sparse = data_->is_sparse();
    }

    CandidateSpace::~CandidateSpace() {
        for (int i = 0; i < MAX_QUERY_EDGE; i++) {
            delete[] BitsetEdgeCS[i];
        }
        delete[] BitsetEdgeCS;
        delete[] num_visit_cs_;
        delete[] visited_candidates_;
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
            memset(BitsetEdgeCS[i], true, data_->edge_info_.size());
        }
        memset(num_visit_cs_, 0, data_->GetNumVertices());
        memset(visited_candidates_, 0, data_->GetNumVertices() * query->GetNumVertices());
        num_visited_candidates = 0;
        BPSolver.global_initialize(query_->GetMaxDegree(), data_->GetMaxDegree());
        BPTSolver.global_initialize(query_->GetMaxDegree(), data_->max_num_trigs);
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            candidate_set_[i].clear();
        }
//        cs_adj_.clear();

        dag_->BuildDAG(-1);
        if (!FilterByTopDownWithInit()) return false;

        Timer csfilter_timer; csfilter_timer.Start();
        auto CandSize = [this]() {
            int totalCandidateSetSize = 0;
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                totalCandidateSetSize += GetCandidateSetSize(i);
            }
            return totalCandidateSetSize;
        };
        int cand_sz = CandSize();
//        global_iteration_count = 0;
//        fprintf(stdout, "Iter %d : Cand_size = %d\n",global_iteration_count, cand_sz);
//        while (true) {
//            global_iteration_count++;
//            bool pruned = false;
//            if (Filter(false)) pruned = true;
//            if (Filter(true)) pruned = true;
//            if (!pruned) break;
//            int new_sz = CandSize();
//            if (new_sz > 0.95 * cand_sz) break;
//            cand_sz = new_sz;
//            fprintf(stdout, "Iter %d : Cand_size = %d\n",global_iteration_count, cand_sz);
//        }

        Filter(false);

        csfilter_timer.Stop();
        fprintf(stdout, "CS prune time : %.02lf ms\n",csfilter_timer.GetTime());
        ConstructCS();
        return true;
    }


    bool CandidateSpace::FilterByTopDownWithInit() {
        bool result = true;

        uint64_t *nbr_label_bitset = new uint64_t[data_->GetNbrBitsetSize()];
        Size max_nbr_degree;
        result = InitRootCandidates();
        if (!result) return false;
        for (Size i = 1; i < query_->GetNumVertices(); ++i) {
            Vertex cur = dag_->GetVertexOrderedByBFS(i);

            Label cur_label = query_->GetLabel(cur);
            QueryDegree num_parent = 0;
            for (Size i = 0; i < dag_->GetNumParents(cur); ++i) {
                Vertex parent = dag_->GetParent(cur, i);
                int query_edge_idx = query_->GetEdgeIndex(parent, cur);

                for (Vertex parent_cand : candidate_set_[parent]) {
                    for (int data_edge_idx : data_->GetIncidentEdges(parent_cand, cur_label)) {
                        Vertex cand = data_->opposite(data_edge_idx, parent_cand);
                        if (num_visit_cs_[cand] < num_parent) continue;
                        if (data_->GetDegree(cand) < query_->GetDegree(cur)) continue;
                        if (data_->GetELabel(data_edge_idx) != query_->GetELabel(query_edge_idx)) continue;
                        if (num_visit_cs_[cand] == num_parent) {
                            num_visit_cs_[cand] += 1;
                            if (num_parent == 0) {
                                visited_candidates_[num_visited_candidates] = cand;
                                num_visited_candidates += 1;
                            }
                        }
                    }
                }
                num_parent += 1;
            }

            ComputeNbrInformation(cur, &max_nbr_degree, nbr_label_bitset);

            for (Size i = 0; i < num_visited_candidates; ++i) {
                Vertex cand = visited_candidates_[i];
                if (num_visit_cs_[cand] == num_parent &&
                    data_->GetCoreNum(cand) >= query_->GetCoreNum(cur) &&
                    data_->GetMaxNbrDegree(cand) >= max_nbr_degree &&
                    data_->CheckAllNbrLabelExist(cand, nbr_label_bitset)) {
                    candidate_set_[cur].emplace_back(cand);
                    BitsetCS[cur][cand] = true;
                }
            }
            if (candidate_set_[cur].empty()) {
                result = false;
                break;
            }

            while (num_visited_candidates > 0) {
                num_visited_candidates -= 1;
                num_visit_cs_[visited_candidates_[num_visited_candidates]] = 0;
            }
        }

        delete[] nbr_label_bitset;
        return result;
    }

    bool CandidateSpace::BipartiteSafety(Vertex cur, Vertex cand) {
#ifdef BIPARTITE_SAFETY
#ifdef TIME_CHECK
        steptimer_1.Start();
#endif
        if (query_->adj_list[cur].size() == 1) {
            int uc = query_->adj_list[cur][0];
            int query_edge_index = query_->GetEdgeIndex(cur, uc);
            for (int edge_id : data_->all_incident_edges_[cand]) {
                Vertex vc = data_->opposite(edge_id, cand);
                if (BitsetCS[uc][vc] and BitsetEdgeCS[query_edge_index][edge_id]) {
#ifdef TIME_CHECK
                    steptimer_1.Stop();
#endif
                    return true;
                }
            }
#ifdef TIME_CHECK
            steptimer_1.Stop();
#endif
            return false;
        }
        BPSolver.reset();
//        for (int i = 0; i < query_->adj_list[cur].size(); i++) {
//            Vertex uc = query_->adj_list[cur][i];
//            int query_edge_index = query_->GetEdgeIndex(cur, uc);
//            j = 0;
        int i = 0, j = 0;
        for (int query_edge_index : query_->all_incident_edges_[cur]) {
            Vertex uc = query_->opposite(query_edge_index, cur);
            j = 0;
            for (int edge_id : data_->all_incident_edges_[cand]) {
                Vertex vc = data_->opposite(edge_id, cand);
                if (BitsetCS[uc][vc] and BitsetEdgeCS[query_edge_index][edge_id]) {
                    BPSolver.add_edge(i, j);
                }
                j++;
            }
            i++;
        }
        bool ok = (BPSolver.solve() == (query_->GetDegree(cur)));
#ifdef TIME_CHECK
        steptimer_1.Stop();
#endif
        return ok;
#else
        return true;
#endif
    }


    bool CandidateSpace::BipartiteEdgeSafety(Vertex cur, Vertex cand, Vertex nxt, Vertex nxt_cand) {
        BPSolver.reset();
        int i = 0, j = 0;
        for (int query_edge_index : query_->all_incident_edges_[cur]) {
            Vertex uc = query_->opposite(query_edge_index, cur);
            j = 0;
            for (int edge_id : data_->all_incident_edges_[cand]) {
                Vertex vc = data_->opposite(edge_id, cand);
                if (uc == nxt) {
                    if (vc == nxt_cand)
                        BPSolver.add_edge(i, j);
                }
                else {
                    if (BitsetCS[uc][vc] and BitsetEdgeCS[query_edge_index][edge_id]) {
                        BPSolver.add_edge(i, j);
                    }
                }
                j++;
            }
            i++;
        }
        bool ok = (BPSolver.solve() == (query_->GetDegree(cur)));
        return ok;
    }
    Size CandidateSpace::GetDAGNextCount(Vertex cur, bool topdown) {
        return topdown? dag_->GetNumParents(cur) : dag_->GetNumChildren(cur);
    }
    Vertex CandidateSpace::GetDAGNextVertex(Vertex cur, Size idx, bool topdown) {
        return topdown? dag_->GetParent(cur, idx) : dag_->GetChild(cur, idx);
    }

    inline bool CandidateSpace::EdgeCandidacy(int query_edge_id, int data_edge_id) {
        if (query_edge_id == -1 || data_edge_id == -1) {
            return false;
        }
        return BitsetEdgeCS[query_edge_id][data_edge_id];
    }

    bool CandidateSpace::TriangleSafety(int query_edge_id, int data_edge_id) {
#ifdef TIME_CHECK
        steptimer_2.Start();
#endif
//        fflush(stdout);
        if (query_->trigvertex[query_edge_id].empty()) return true;
        auto &common_neighbors = data_->trigvertex[data_edge_id];
        if (query_->trigvertex[query_edge_id].size() > common_neighbors.size()) return false;
        if (query_->trigvertex[query_edge_id].size() == 1) {
            auto qtv = *(query_->trigvertex[query_edge_id].begin());
            for (auto tv: common_neighbors) {
                if (!BitsetCS[qtv.first][tv.first]) continue;
                if (!BitsetEdgeCS[qtv.second.first][tv.second.first]) continue;
                if (!BitsetEdgeCS[qtv.second.second][tv.second.second]) continue;
#ifdef TIME_CHECK
                steptimer_2.Stop();
#endif
                return true;
            }
#ifdef TIME_CHECK
            steptimer_2.Stop();
#endif
            return false;
        }
//        else return true;
//        fprintf(stderr, "Trianglesafety %d %d : Size %lu vs %lu MBP in\n", query_edge_id, data_edge_id, query_->trigvertex[query_edge_id].size(), common_neighbors.size());

        BPTSolver.reset();
        int triangle_index = -1;
        for (auto qtv: query_->trigvertex[query_edge_id]) {
            triangle_index++;
            int tv_idx = -1;
            bool found = false;
            for (auto tv: common_neighbors) {
                tv_idx++;
                if (!BitsetCS[qtv.first][tv.first]) continue;
                if (!BitsetEdgeCS[qtv.second.first][tv.second.first]) continue;
                if (!BitsetEdgeCS[qtv.second.second][tv.second.second]) continue;
                found = true;
                BPTSolver.add_edge(triangle_index, tv_idx);
            }
            if (!found) {
#ifdef TIME_CHECK
                steptimer_2.Stop();
#endif
                return false;
            }
        }
        bool ok = (BPTSolver.solve() == query_->triangles[query_edge_id].size());
#ifdef TIME_CHECK
        steptimer_2.Stop();
#endif
        return ok;
    }

    bool CandidateSpace::FourCycleSafetyOnline(int query_edge_id, int data_edge_id) {
        int cand, nxt_cand; std::tie(cand, nxt_cand) = data_->edge_info_[data_edge_id].vp;
        int cur, nxt; std::tie(cur, nxt) = query_->edge_info_[query_edge_id].vp;
        std::pair<int, int> edgePair = {query_edge_id, data_edge_id};
        // 4-cycle filter
        if (four_cycle_memo_old.find(edgePair) == four_cycle_memo_old.end()) {
            four_cycle_memo_old[edgePair] = std::vector<online_cycle_information>(query_->four_cycles[query_edge_id].size());
        }
        std::vector<online_cycle_information> &four_cycle_answers = four_cycle_memo_old[edgePair];
        for (int i = 0; i < query_->four_cycles[query_edge_id].size(); i++) {
            int opp_edge_idx = query_->four_cycles[query_edge_id][i].opp_edge_idx;
            auto &info = four_cycle_answers[i];
            bool validity = true;
//            steptimer_4.Start();
            if (info.d_opp_edge_idx == -1) {
                validity = false;
//                fprintf(stderr, "Solving : 4-Cycle QEdge [%d]-[%d] and dataedge [%d]-[?] for the first time\n",
//                        query_edge_id, opp_edge_idx, data_edge_id);
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
                counters[16]++;
//                steptimer_4.Stop();
                continue;
            }

//            steptimer_3.Start();
            Label third_label = query_->GetLabel(info.third);
            Label fourth_label = query_->GetLabel(info.fourth);
            std::vector <int> &third_cands  = data_->GetIncidentEdges(nxt_cand, third_label);
            std::vector <int> &fourth_cands = data_->GetIncidentEdges(cand, fourth_label);
            if (info.third_inc_idx >= third_cands.size()) return false;

            while (info.third_inc_idx < third_cands.size()) {
                counters[15]++;
                bool cand_validity = true;
                int third_cand = data_->opposite(third_cands[info.third_inc_idx], nxt_cand);
                int fourth_cand = data_->opposite(fourth_cands[info.fourth_inc_idx], cand);
                if (cand_validity) cand_validity &= ((fourth_cand != nxt_cand) && (third_cand != cand));
                if (cand_validity) cand_validity &= BitsetCS[info.third][third_cand];
                if (cand_validity) cand_validity &= EdgeCandidacy(info.q_third_edge_idx, third_cands[info.third_inc_idx]);
                if (cand_validity) cand_validity &= BitsetCS[info.fourth][fourth_cand];
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
//                steptimer_3.Stop();
                return false;
            }
            ok:
//            steptimer_3.Stop();
            continue;
        }
        return true;
    }

    bool CandidateSpace::FourCycleSafety(int query_edge_id, int data_edge_id) {
        if (data_->num_four_cycles_indexed == 0) return FourCycleSafetyOnline(query_edge_id, data_edge_id);
#ifdef TIME_CHECK
        steptimer_3.Start();
#endif
        if (query_->four_cycles[query_edge_id].size() > data_->four_cycles[data_edge_id].size()) return false;
//        BPQsolver.reset();
        for (int i = 0; i < query_->four_cycles[query_edge_id].size(); i++) {
            auto &q_info = query_->four_cycles[query_edge_id][i];
            bool found = false;
            for (int j = 0; j < data_->four_cycles[data_edge_id].size(); j++) {
                auto &d_info =  data_->four_cycles[data_edge_id][j];
                bool validity = true;
                validity &= BitsetCS[q_info.third][d_info.third];
                validity &= BitsetCS[q_info.fourth][d_info.fourth];
                if (!validity) continue;
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
//                    BPQsolver.add_edge(i, j);
                    goto nxt_cycle;
                }
            }
            if (!found) return false;
            nxt_cycle:
#ifdef TIME_CHECK
            steptimer_3.Stop();
#endif
            continue;
        }
//        bool ok = (BPQsolver.solve() == (query_->four_cycles[query_edge_id].size()));
        return true;
    }

    bool CandidateSpace::EdgeSafety(int query_edge_id, int data_edge_id) {
        // triangle filter
        if (is_data_sparse and !query_->trig_empty[query_edge_id]) {
            if (!TriangleSafety(query_edge_id, data_edge_id))
                return false;
        }
#ifdef FOURCYCLE_SAFETY
#ifdef TIME_CHECK
        steptimer_3.Start();
#endif
        if (is_data_sparse and !query_->quad_empty[query_edge_id]) {
            if (!FourCycleSafety(query_edge_id, data_edge_id)) return false;
        }
#ifdef TIME_CHECK
        steptimer_3.Stop();
#endif
#endif
        return true;
    }

    struct variable {
        double priority;
        int stage, which;
        bool operator>(const variable &o) const {
            return priority > o.priority;
        }
        bool operator<(const variable &o) const {
            return priority < o.priority;
        }
    };

    bool CandidateSpace::Filter(bool topdown) {
        auto CandSize = [this]() {
            int totalCandidateSetSize = 0;
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                totalCandidateSetSize += GetCandidateSetSize(i);
            }
            return totalCandidateSetSize;
        };
        int cand_sz = CandSize();
        fprintf(stderr, "Initial CS SIZE = %d\n",cand_sz);

        std::vector<int> local_stage(query_->GetNumVertices(), 0);
        std::vector<double> priority(query_->GetNumVertices(), 0.5);

        std::priority_queue<variable> candidate_queue;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            candidate_queue.push({priority[i], 0, i});
        }
        int queue_pop_count = 0;
        int maximum_queue_cnt = 5 * query_->GetNumVertices();
        int current_stage = 0;
        while (!candidate_queue.empty()) {
            auto [pri, stage, cur] = candidate_queue.top();
            candidate_queue.pop();
            if (stage < local_stage[cur]) continue;
            current_stage++;
            queue_pop_count++;

            int bef_cand_size = candidate_set_[cur].size();
//            fprintf(stderr, "Got vertex %d to reduce... Candidate Size of Vertex %d at "
//                            "stage %d : %d\n",cur, cur, current_stage, bef_cand_size);
            for (int i = 0; i < candidate_set_[cur].size(); i++) {
                Vertex cand = candidate_set_[cur][i];
                bool valid = true;
                int ii = 0, jj = 0;

                // 3, 4 Cycle Filter
                for (int query_edge_idx : query_->all_incident_edges_[cur]) {
                    int nxt = query_->to_[query_edge_idx];
                    int nxt_label = query_->GetLabel(nxt);
                    bool found = false;
                    for (int data_edge_idx : data_->GetIncidentEdges(cand, nxt_label)) {
                        Vertex nxt_cand = data_->to_[data_edge_idx];
                        if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                        if (!BitsetCS[nxt][nxt_cand]) continue;
                        if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
                            continue;
                        }
                        if (!EdgeSafety(query_edge_idx, data_edge_idx)) {
                            BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                            BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]] = false;
                            continue;
                        }
                        found = true;
//                        break;
                    }
                    if (!found) {
                        valid = false;
                        break;
                    }
                }
                if (!valid) goto remove_vertex_;
                if (query_->adj_list[cur].size() == 1) goto remove_vertex_;
                // Edge-Bipartite
                BPSolver.reset();
                // add all edges
                for (int query_edge_index : query_->all_incident_edges_[cur]) {
                    Vertex uc = query_->opposite(query_edge_index, cur);
//                    printf("Neighbor %d of %d is %dth neighbor\n",uc,cur,ii);
                    jj = 0;
                    for (int edge_id : data_->all_incident_edges_[cand]) {
                        Vertex vc = data_->opposite(edge_id, cand);
                        if (data_->GetDegree(vc) < query_->GetDegree(uc)) break;
//                        printf("Neighbor %d of %d is %dth neighbor\n",vc,cand,jj);
                        if (BitsetCS[uc][vc] and BitsetEdgeCS[query_edge_index][edge_id]) {
                            BPSolver.add_edge(ii, jj);
                        }
                        jj++;
                    }
                    ii++;
                }
                ii = 0;
                for (int query_edge_idx : query_->all_incident_edges_[cur]) {
                    int nxt = query_->to_[query_edge_idx];
                    bool found = false;

                    BPSolver.reset(false);
                    if (BPSolver.solve(ii) < query_->adj_list[cur].size() - 1) {
                        valid = false;
                        break;
                    }
                    jj = -1;
                    std::memcpy(temporary_left, BPSolver.left, sizeof(int) * BPSolver.left_len);
                    for (int data_edge_idx : data_->all_incident_edges_[cand]) {
                        Vertex nxt_cand = data_->to_[data_edge_idx];
                        jj++;
                        if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                        if (!BitsetCS[nxt][nxt_cand]) continue;
//                        if (!BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]]) {
//                            BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
//                        }
                        if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
                            continue;
                        }
//                        printf("Checking Bipartite Matching with %d-%d fixed from (%d, %d)\n",nxt,nxt_cand,cur,cand);
                        if (BPSolver.right[jj] == -1) {
//                            printf("Right[%d] = -1! ok\n",jj);
                            found = true;
                            continue;
                        }
                        std::memset(BPSolver.used, false, sizeof(bool) * query_->adj_list[cur].size());
                        BPSolver.used[ii] = true;
//                        printf("ii = %d(%d), jj = %d(%d), DFS[%d]...\n",
//                               ii,query_->to_[query_->all_incident_edges_[cur][ii]],
//                               jj,data_->to_[data_->all_incident_edges_[cand][jj]],BPSolver.right[jj]);
//                        BPSolver.print();
                        if (BPSolver.single_dfs(BPSolver.right[jj])) {
//                            printf("DFS[%d]...\n",BPSolver.right[jj]);
                            found = true;
                        }
                        else {
//                            printf("Edge %d-%d is dead\n",query_edge_idx, data_edge_idx);
//                            printf("%d-%d to %d-%d\n",cur, nxt, cand, nxt_cand);
                            BitsetEdgeCS[query_edge_idx][data_edge_idx] = false;
                            BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]] = false;
                        }
                    }
                    if (!found) {
//                        printf("NO valid EdgeBP! [%d, %d] dead\n",cur,cand);
                        valid = false;
                        break;
                    }
                    found = false;
                    ii++;
                }
                remove_vertex_:
                if (!valid) {
                    candidate_set_[cur][i] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    --i;
                    BitsetCS[cur][cand] = false;
                }
            }

            if (candidate_set_[cur].empty()) {
                counters[191]++;
                if (counters[191] == 1) {
                    fprintf(stderr, "FOUND INVALID %u\n",cur);
                    exit(2);
                }
                return true;
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
                priority[nxt] = 1 - (1 - out_prob) * (1 - priority[nxt]);
                if (priority[nxt] < 0.1) continue;
                local_stage[nxt] = current_stage;
                //fprintf(stderr, "    Cur=[%d] pushes Nxt=[%d] with priority %lf\n",cur,nxt,priority[nxt]);
                candidate_queue.push({priority[nxt], current_stage, nxt});
            }
        }
        return true;
    }

    void CandidateSpace::ConstructCS() {
        Timer csbuild_timer; csbuild_timer.Start();
        cs_edge_list_.clear();
        cs_edge_list_.resize(query_->GetNumVertices(), std::vector<std::vector<VertexPair>>());
        std::vector <int> CandidateIndex(data_->GetNumVertices());
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            cs_edge_list_[i].resize(GetCandidateSetSize(i));
        }
        Size CandidateSpaceSize = 0, CandidateEdges = 0;
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            Vertex u = dag_->GetVertexOrderedByBFS(i);
            Label u_label = query_->GetLabel(u);
            Size u_degree = query_->GetDegree(u);

            CandidateSpaceSize += GetCandidateSetSize(u);
            for (Size idx = 0; idx < GetCandidateSetSize(u); idx++)
                CandidateIndex[candidate_set_[u][idx]] = idx;

            for (Size i = query_->GetStartOffset(u); i < query_->GetEndOffset(u); ++i) {
                Vertex u_adj = query_->GetNeighbor(i);
                int query_edge_idx = query_->GetEdgeIndex(u_adj, u);

                for (Size v_adj_idx = 0; v_adj_idx < candidate_set_[u_adj].size(); ++v_adj_idx) {
                    Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

                    for (int data_edge_idx : data_->GetIncidentEdges(v_adj, u_label)) {
                        Vertex v = data_->opposite(data_edge_idx, v_adj);

                        if (data_->GetDegree(v) < u_degree) break;
                        if (!BitsetEdgeCS[query_edge_idx][data_edge_idx]) continue;
                        if (!BitsetEdgeCS[query_->opposite_edge[query_edge_idx]][data_->opposite_edge[data_edge_idx]]) continue;
//                        if (!EdgeSafety(query_edge_idx, data_edge_idx)) continue;
                        if (BitsetCS[u][v]){
                            CandidateEdges++;
                            Size v_idx = CandidateIndex[v];
                            VertexPair uc_pair = {u_adj, v_adj_idx};
                            cs_edge_list_[u][v_idx].emplace_back(uc_pair);
//                            cs_edge_list_[uc_pair].insert(u_pair);

                            /* Add CS_ADJ_List; */
//                            cs_adj_[{u, v}][u_adj].emplace_back(v_adj);
                        }
                    }
                }
            }
        }
        cs_edge_.clear();
        cs_edge_.resize(query_->GetNumVertices());
        for (int u = 0; u < query_->GetNumVertices(); u++) {
            cs_edge_[u].resize(GetCandidateSetSize(u));
            for (Size idx = 0; idx < GetCandidateSetSize(u); idx++) {
                auto &it = cs_edge_list_[u][idx];
                std::sort(it.begin(), it.end());
                it.erase(std::unique(it.begin(), it.end()), it.end());
                cs_edge_[u][idx].resize(query_->GetNumVertices());
                for (auto &i : it) {
                    cs_edge_[u][idx][i.first].push_back(i.second);
                }
            }
        }

        CandidateEdges /= 2;
        csbuild_timer.Stop();

        fprintf(stdout, "Counter 15 = %llu, 16 = %llu, 17 = %llu\n",counters[15], counters[16], counters[17]);
        std::cout << "GetEdgeIndex Calls : " << functionCallCounter << std::endl;
        std::cout << "StepTimer 1 (Bipartite-Safety) " << steptimer_1.GetTime() << "ms" << std::endl;
        std::cout << "StepTimer 2 (Triangle-Safety) " << steptimer_2.GetTime() << "ms" << std::endl;
        std::cout << "StepTimer 3 (4Cycle-Safety) " << steptimer_3.GetTime() << "ms" << std::endl;
        std::cout << "StepTimer 4 (4Cycle-Safety-Bypass) " << steptimer_4.GetTime() << "ms" << std::endl;
        fprintf(stdout, "Total 4-cycle checks = %llu, Passed %llu\n",counters[0], counters[1]);
        fprintf(stdout, "CS build time : %.02lf ms\n",csbuild_timer.GetTime());
        fprintf(stderr, "CS-SZ : %u %u\n", CandidateSpaceSize, CandidateEdges);
        fprintf(stdout, "#CandidateSetSize : %u\n", CandidateSpaceSize);
        fprintf(stdout, "#CandidateSetEdges : %u\n", CandidateEdges);
        double vert_cand_avg = (1.0*CandidateSpaceSize)/query_->GetNumVertices();
        fprintf(stdout, "#AVG_VERT_CAND_SIZE : %.02lf\n", vert_cand_avg);
        fprintf(stdout, "#AVG_EDGE_BTW_CANDS : %.02lf\n", (1.0*CandidateEdges)/query_->GetNumEdges()/vert_cand_avg);
////        printCS();
//        for (int i = 0; i < query_->GetNumVertices(); i++) {
//            fprintf(stdout, "Query [%d] : CS Size = %lu\n", i,candidate_set_[i].size());
////            for (int x: candidate_set_[i]) {
////                fprintf(stdout, "%d ", x);
////            }
////            fprintf(stdout, "\n");
//        }
    }

    bool CandidateSpace::InitRootCandidates() {
        Vertex root = dag_->GetRoot();
        Label root_label = query_->GetLabel(root);

        uint64_t *nbr_label_bitset = new uint64_t[data_->GetNbrBitsetSize()];
        Size max_nbr_degree;

        ComputeNbrInformation(root, &max_nbr_degree, nbr_label_bitset);

        for (Size i = data_->GetStartOffsetByLabel(root_label);
             i < data_->GetEndOffsetByLabel(root_label); ++i) {
            Vertex cand = data_->GetVertexBySortedLabelOffset(i);

            if (data_->GetDegree(cand) < query_->GetDegree(root)) continue;

            if (data_->GetCoreNum(cand) >= query_->GetCoreNum(root) &&
                data_->CheckAllNbrLabelExist(cand, nbr_label_bitset) &&
                data_->GetMaxNbrDegree(cand) >= max_nbr_degree) {
                candidate_set_[root].emplace_back(cand);
                BitsetCS[root][cand] = true;
            }
        }

        delete[] nbr_label_bitset;
        return !candidate_set_[root].empty();
    }

    void CandidateSpace::ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                               uint64_t *nbr_label_bitset) {
        *max_nbr_degree = 0;
        std::fill(nbr_label_bitset, nbr_label_bitset + data_->GetNbrBitsetSize(),
                  0ull);
        for (Size i = query_->GetStartOffset(u); i < query_->GetEndOffset(u); ++i) {
            Vertex adj = query_->GetNeighbor(i);

            nbr_label_bitset[query_->GetLabel(adj) / (sizeof(uint64_t) * CHAR_BIT)] |=
                    1ull << (query_->GetLabel(adj) % (sizeof(uint64_t) * CHAR_BIT));

            if (query_->GetDegree(adj) > *max_nbr_degree) {
                *max_nbr_degree = query_->GetDegree(adj);
            }
        }
    }

    void CandidateSpace::printCS() {
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            fprintf(stdout, "Query [%d] : ",i);
            for (int x : candidate_set_[i]) {
                fprintf(stdout, "%d ", x);
            }
            fprintf(stdout, "\n");
            for (int j = 0; j < GetCandidateSetSize(i); j++) {
                fprintf(stdout, "  CS Edges from %d\n", GetCandidate(i, j));
                for (auto it : cs_edge_list_[i][j]) {
                    fprintf(stdout, "    Edge with [%d, %d]\n",it.first, it.second);
                }
            }
        }
    }
}