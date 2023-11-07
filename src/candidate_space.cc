#include <vector>
#include <map>
#include <random>
#include "include/candidate_space.h"
#include "global/timer.h"

namespace CardEst {

    BipartiteMaximumMatching BPSolver;

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

    bool CandidateSpace::BuildCS(QueryGraph *query, OrderedQueryGraph *dag) {
        query_ = query;
        dag_ = dag;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            memset(BitsetCS[i], false, data_->GetNumVertices());
        }
        for (int i = 0; i < query_->edge_info_.size(); i++) {
            memset(BitsetEdgeCS[i], false, data_->edge_info_.size());
        }
        memset(num_visit_cs_, 0, data_->GetNumVertices());
        if (opt.neighborhood_filter >= NEIGHBOR_BIPARTITE_SAFETY) {
            BPSolver.global_initialize(query_->GetMaxDegree(), data_->GetMaxDegree());
        }
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            candidate_set_[i].clear();
        }
        dag_->BuildOrderedQueryGraph(-1);
        if (!InitializeCS()) {
            exit(1);
        }
        Timer csfilter_timer; csfilter_timer.Start();
        Filter();
        csfilter_timer.Stop();
        ConstructCS();
        return true;
    }

    bool CandidateSpace::InitializeCS() {
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

            for (Size j = 0; j < candidate_set_[cur].size(); ++j) {
                Vertex cand = candidate_set_[cur][j];

                BitsetCS[cur][cand] = false;
                if (num_visit_cs_[cand] == num_parent) {
                    BitsetCS[cur][cand] = true;
                }
                else {
                    candidate_set_[cur][j] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    j--;
                }
                num_visit_cs_[cand] = 0;
            }

            if (candidate_set_[cur].empty()) exit(1);
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
                            cs_edge++;
                        }
                    }
                }
            }
            cs_vertex += candidate_set_[i].size();
        }

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

    bool CandidateSpace::FourCycleSafety(int query_edge_id, int data_edge_id) {
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
        return true;
    }
    bool CandidateSpace::StructureSafety(int query_edge_id, int data_edge_id) {
        switch (opt.structure_filter) {
            case NO_STRUCTURE_FILTER:
                return true;
            case TRIANGLE_SAFETY:
                return TriangleSafety(query_edge_id, data_edge_id);
            case FOURCYCLE_SAFETY:
                return TriangleSafety(query_edge_id, data_edge_id) and FourCycleSafety(query_edge_id, data_edge_id);
        }
        return true;
    }
    bool CandidateSpace::StructureFilter(int cur, int cand, int direction) {
        std::vector<int>* query_neighbors = &(query_->all_incident_edges_[cur]);
        if (direction == -1) {
            query_neighbors = &(dag_->dag_child_edges_[cur]);
        }
        else if (direction == 1) {
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
                int left_unmatch = query_->all_incident_edges_[cur][i];
                int right_unmatch = data_->all_incident_edges_[cand][j];
                BitsetEdgeCS[left_unmatch][right_unmatch] = false;
                BitsetEdgeCS[query_->opposite_edge[left_unmatch]][data_->opposite_edge[right_unmatch]] = false;
            }
        }

        return true;
    }

    bool CandidateSpace::NeighborhoodFilter(int cur, int cand) {
        switch (opt.neighborhood_filter) {
            case NEIGHBOR_SAFETY:
                return CheckNeighborSafety(cur, cand);
            case NEIGHBOR_BIPARTITE_SAFETY:
                return BipartiteSafety(cur, cand);
            case EDGE_BIPARTITE_SAFETY:
                return EdgeBipartiteSafety(cur, cand);
        }
    }
    bool CandidateSpace::Filter() {
        std::vector<int> local_stage(query_->GetNumVertices(), 0);
        std::vector<double> priority(query_->GetNumVertices(), 0);
        std::priority_queue<QueryVertex> candidate_queue;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            priority[i] = 0.5;
            candidate_queue.push({priority[i], 0, i});
        }

        int queue_pop_count = 0;
        int maximum_queue_cnt = 5 * query_->GetNumVertices();
        int current_stage = 0;
        if (opt.refinement_order == DAG_DP) goto dag_dp;
        while (!candidate_queue.empty() and queue_pop_count <= maximum_queue_cnt) {
            auto [pri, stage, cur] = candidate_queue.top();
            candidate_queue.pop();
            if (stage < local_stage[cur]) continue;
            current_stage++;
            queue_pop_count+=query_->GetDegree(cur);
            int bef_cand_size = candidate_set_[cur].size();

            if (opt.neighborhood_filter == NEIGHBOR_SAFETY) {
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
                if (valid) valid = NeighborhoodFilter(cur, cand);
                if (!valid) {
                    int removed = candidate_set_[cur][i];
                    for (int query_edge_idx : query_->all_incident_edges_[cur]) {
                        int nxt = query_->to_[query_edge_idx];
                        for (int data_edge_idx : data_->GetIncidentEdges(removed, query_->GetLabel(nxt))) {
                            int nxt_cand = data_->to_[data_edge_idx];
                            if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                            if (BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
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

            if (candidate_set_[cur].empty()) exit(1);

            int aft_cand_size = candidate_set_[cur].size();
            if (aft_cand_size == bef_cand_size) {
                priority[cur] = 0;
                continue;
            }
            double out_prob = 1 - aft_cand_size * 1.0 / bef_cand_size;
            priority[cur] = 0;

            for (Vertex nxt : query_->adj_list[cur]) {
                priority[nxt] = 1 - (1 - out_prob) * (1 - priority[nxt]);
                if (priority[nxt] < opt.cutoff) continue;
                local_stage[nxt] = current_stage;

                candidate_queue.push({priority[nxt], current_stage, nxt});
            }
        }
        return true;

        dag_dp:
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
                    if (valid) valid = NeighborhoodFilter(cur, cand);
                    if (!valid) {
                        int removed = candidate_set_[cur][i];
                        for (int query_edge_idx : query_->all_incident_edges_[cur]) {
                            int nxt = query_->to_[query_edge_idx];
                            for (int data_edge_idx : data_->GetIncidentEdges(removed, query_->GetLabel(nxt))) {
                                int nxt_cand = data_->to_[data_edge_idx];
                                if (data_->GetDegree(nxt_cand) < query_->GetDegree(nxt)) break;
                                if (BitsetEdgeCS[query_edge_idx][data_edge_idx]) {
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

                if (candidate_set_[cur].empty()) exit(1);
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
        for (Vertex u = 0; u < query_->GetNumVertices(); u++) {
            Label u_label = query_->GetLabel(u);
            Size u_degree = query_->GetDegree(u);

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
                        cs_edge_[u][CandidateIndex[v]][u_adj].emplace_back(v_adj_idx);
                    }
                }
            }
        }
        csbuild_timer.Stop();
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

    Size CandidateSpace::GetDAGNextCount(Vertex cur, bool topdown) {
        return topdown? dag_->GetNumParents(cur) : dag_->GetNumChildren(cur);
    }
}