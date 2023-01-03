#include <vector>
#include <map>
#include <random>
#include "include/daf_candidate_space.h"
#include "global/timer.h"
#define BIPARTITE_SAFETY
//#define NEIGHBOR_SAFETY
#define FOURCYCLE_SAFETY

namespace daf {
    Size global_iteration_count = 0;
    Timer steptimer_1, steptimer_2, steptimer_3, steptimer_4;
    Size counters[1000];
    CandidateSpace::CandidateSpace(DataGraph &data, QueryGraph &query,
                                   DAG &dag)
            : data_(data), query_(query), dag_(dag) {

        BitsetCS.resize(query_.GetNumVertices());
        candidate_set_.resize(query_.GetNumVertices());

        for (Vertex u = 0; u < query_.GetNumVertices(); ++u) {
            BitsetCS[u].resize(data_.GetNumVertices());
        }
        num_visit_cs_ = new QueryDegree[data_.GetNumVertices()];
        visited_candidates_ =
                new Vertex[data_.GetNumVertices() * query_.GetNumVertices()];

        num_visited_candidates = 0;

        std::fill(num_visit_cs_, num_visit_cs_ + data_.GetNumVertices(), 0);
        num_cs_edges_ = 0;
        for (EdgeInfo &e : data_.edge_info_) {
            e.edge_candidacy_.clear();
            e.edge_candidacy_.resize(query_.GetNumEdges()*2, true);
        }
    }

    CandidateSpace::~CandidateSpace() {
        delete[] num_visit_cs_;
        delete[] visited_candidates_;
    }

    bool CandidateSpace::BuildCS() {

        srand(0);
        dag_.BuildDAG(-1);
        if (!FilterByTopDownWithInit()) return false;

        Timer csfilter_timer; csfilter_timer.Start();
        auto CandSize = [this]() {
            int totalCandidateSetSize = 0;
            for (int i = 0; i < query_.GetNumVertices(); i++) {
                totalCandidateSetSize += GetCandidateSetSize(i);
            }
            return totalCandidateSetSize;
        };
        int cand_sz = CandSize();
        global_iteration_count = 0;
        fprintf(stdout, "Iter %d : Cand_size = %d\n",global_iteration_count, cand_sz);
        while (true) {
            global_iteration_count++;
            bool pruned = false;
            if (Filter(false)) pruned = true;
            if (Filter(true)) pruned = true;
            if (!pruned) break;
            int new_sz = CandSize();
            if (new_sz > 0.95 * cand_sz) break;
            cand_sz = new_sz;
            fprintf(stdout, "Iter %d : Cand_size = %d\n",global_iteration_count, cand_sz);
        }
        csfilter_timer.Stop();
        fprintf(stdout, "CS prune time : %.02lf ms\n",csfilter_timer.GetTime());
        fflush(stdout);
        ConstructCS();
        return true;
    }


    bool CandidateSpace::FilterByTopDownWithInit() {
        bool result = true;

        uint64_t *nbr_label_bitset = new uint64_t[data_.GetNbrBitsetSize()];
        Size max_nbr_degree;
        result = InitRootCandidates();
        if (!result) return false;
        for (Size i = 1; i < query_.GetNumVertices(); ++i) {
            Vertex cur = dag_.GetVertexOrderedByBFS(i);

            Label cur_label = query_.GetLabel(cur);
            QueryDegree num_parent = 0;
            for (Size i = 0; i < dag_.GetNumParents(cur); ++i) {
                Vertex parent = dag_.GetParent(cur, i);
                int query_edge_idx = query_.GetEdgeIndex(parent, cur);

                for (Vertex parent_cand : candidate_set_[parent]) {
                    for (int data_edge_idx : data_.GetIncidentEdges(parent_cand, cur_label)) {
                        Vertex cand = data_.opposite(data_edge_idx, parent_cand);
                        if (data_.GetDegree(cand) < query_.GetDegree(cur)) continue;
                        if (data_.GetELabel(data_edge_idx) != query_.GetELabel(query_edge_idx)) continue;
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
                    data_.GetCoreNum(cand) >= query_.GetCoreNum(cur) &&
                    data_.GetMaxNbrDegree(cand) >= max_nbr_degree &&
                    data_.CheckAllNbrLabelExist(cand, nbr_label_bitset)) {
                    candidate_set_[cur].push_back(cand);
                    BitsetCS[cur].set(cand);
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
//        if (data_.GetDegree(cand) >= 5 * query_.GetDegree(cur)) {
//            return true;
//        }
        steptimer_1.Start();
        std::vector <int> label_edge_offset(data_.GetNumLabels()+1);
        for (int i = 0; i < data_.GetNumLabels(); i++) {
            label_edge_offset[i+1] = data_.GetIncidentEdges(cand, i).size();
            if (i > 0) label_edge_offset[i] += label_edge_offset[i-1];
        }
        BipartiteMaximumMatching BP(query_.GetDegree(cur), data_.GetDegree(cand));
        BP.clear_adj();
        for (int i = 0; i < query_.adj_list[cur].size(); i++) {
            Vertex uc = query_.adj_list[cur][i];
            int query_edge_index = query_.GetEdgeIndex(cur, uc);
            int query_next_label = query_.GetLabel(uc);
            int j = 0;
            for (int edge_id : data_.GetIncidentEdges(cand, query_next_label)) {
                Vertex vc = data_.opposite(edge_id, cand);
                if (BitsetCS[uc][vc] and EdgeCandidacy(query_edge_index, edge_id)) {
                    BP.add_edge(i, j + label_edge_offset[query_next_label]);
                }
                j++;
            }
        }
        bool ok = (BP.solve() == (query_.GetDegree(cur)));
        steptimer_1.Stop();
        return ok;
#else
        return true;
#endif
    }

    Size CandidateSpace::GetDAGNextCount(Vertex cur, bool topdown) {
        return topdown? dag_.GetNumParents(cur) : dag_.GetNumChildren(cur);
    }
    Vertex CandidateSpace::GetDAGNextVertex(Vertex cur, Size idx, bool topdown) {
        return topdown? dag_.GetParent(cur, idx) : dag_.GetChild(cur, idx);
    }

    void CandidateSpace::PrepareNeighborSafety(Vertex cur) {
#ifdef NEIGHBOR_SAFETY
        std::fill(neighbor_label_frequency.begin(), neighbor_label_frequency.end(), 0);
        std::fill(in_neighbor_cs.begin(), in_neighbor_cs.end(), 0);
        for (Size nidx = query_.GetStartOffset(cur); nidx < query_.GetEndOffset(cur); ++nidx) {
            Vertex q_neighbor = query_.GetNeighbor(nidx);
            neighbor_label_frequency[query_.GetLabel(q_neighbor)]++;
            for (Vertex d_neighbor : candidate_set_[q_neighbor]) {
                in_neighbor_cs[d_neighbor] = true;
            }
        }
#endif
    }

    bool CandidateSpace::CheckNeighborSafety(Vertex cur, Vertex cand) {
#ifdef NEIGHBOR_SAFETY
        auto tmp_label_frequency = neighbor_label_frequency;
        for (Size nidx = data_.GetStartOffset(cand); nidx < data_.GetEndOffset(cand); nidx++) {
            Vertex d_neighbor = data_.GetNeighbor(nidx);
            if (in_neighbor_cs[d_neighbor]) {
                tmp_label_frequency[data_.GetLabel(d_neighbor)]--;
            }
        }
        for (int l = 0; l < data_.GetNumLabels(); ++l) {
            if (tmp_label_frequency[l] > 0) {
                return false;
            }
        }
#endif
        return true;
    }

    inline bool CandidateSpace::EdgeCandidacy(int query_edge_id, int data_edge_id) {
        if (query_edge_id == -1 || data_edge_id == -1) return false;
        return data_.edge_info_[data_edge_id].edge_candidacy_[query_edge_id];
    }

    bool CandidateSpace::TriangleSafety(int query_edge_id, int data_edge_id) {
        steptimer_2.Start();
        int cand, nxt_cand; std::tie(cand, nxt_cand) = data_.edge_info_[data_edge_id].vp;
        int cur, nxt; std::tie(cur, nxt) = query_.edge_info_[query_edge_id].vp;
        auto &common_neighbors = data_.trigvertex[data_edge_id];
        BipartiteMaximumMatching BP(query_.triangles[query_edge_id].size(), common_neighbors.size());
        int triangle_index = -1;
        for (auto qtv: query_.trigvertex[query_edge_id]) {
            triangle_index++;
            int tv_idx = -1;
            bool found = false;
            for (auto tv: common_neighbors) {
                tv_idx++;
                if (!BitsetCS[qtv.first][tv.first]) continue;
                if (!EdgeCandidacy(qtv.second.first, tv.second.first)) continue;
                if (!EdgeCandidacy(qtv.second.second, tv.second.second)) continue;
                found = true;
                BP.add_edge(triangle_index, tv_idx);
            }
            if (!found) {
                steptimer_2.Stop();
                return false;
            }
        }
        bool ok = (BP.solve() == query_.triangles[query_edge_id].size());
        steptimer_2.Stop();
        return ok;
    }

//
//    bool CandidateSpace::FourCycleSafety(int query_edge_id, int data_edge_id) {
//        steptimer_3.Start();
//        int cand, nxt_cand; std::tie(cand, nxt_cand) = data_.edge_info_[data_edge_id].vp;
//        int cur, nxt; std::tie(cur, nxt) = query_.edge_info_[query_edge_id].vp;
//        std::pair<int, int> edgePair = {query_edge_id, data_edge_id};
//        // 4-cycle filter
//        bool first_index = false;
//
//        if (four_cycle_memo_old.find(edgePair) == four_cycle_memo_old.end()) {
//            first_index = true;
//            four_cycle_memo_old[edgePair] = std::vector<std::vector<std::tuple<int, int, int>>>(query_.four_cycles[query_edge_id].size());
//        }
//        std::vector<std::vector <std::tuple<int, int, int>>> &four_cycle_answers = four_cycle_memo_old[edgePair];
//        for (int i = 0; i < query_.four_cycles[query_edge_id].size(); i++) {
//            int opp_edge_idx = query_.four_cycles[query_edge_id][i];
//            VertexPair opp_edge = query_.edge_info_[opp_edge_idx].vp;
//            if (opp_edge.first == opp_edge.second) continue;
//            Vertex third = opp_edge.first, fourth = opp_edge.second;
//            int query_fourth_edge_ = query_.GetEdgeIndex(cur, fourth);
//            Label fourth_label = query_.GetLabel(fourth);
//            int query_third_edge_ = query_.GetEdgeIndex(nxt, third);
//            Label third_label = query_.GetLabel(third);
//
//            counters[0]++;
//            if (first_index) {
//                steptimer_3.Start();
//
//                for (auto fourth_cand_edge : data_.GetIncidentEdges(cand, fourth_label)) {
//                    if (!EdgeCandidacy(query_fourth_edge_, fourth_cand_edge)) continue;
//                    Vertex fourth_cand = data_.opposite(fourth_cand_edge, cand);
//                    if (BitsetCS[fourth][fourth_cand] and fourth_cand != nxt_cand)
//                        tmpBitset[fourth_cand] = true;
//                }
//
//                for (auto third_cand_edge : data_.GetIncidentEdges(nxt_cand, third_label)) {
//                    if (!EdgeCandidacy(query_third_edge_, third_cand_edge)) continue;
//                    Vertex third_cand = data_.opposite(third_cand_edge, nxt_cand);
//                    if (BitsetCS[third][third_cand] and third_cand != cand) {
//                        for (auto ft_edge : data_.GetIncidentEdges(third_cand, fourth_label)) {
//                            if (!EdgeCandidacy(opp_edge_idx, ft_edge)) continue;
//                            int ft = data_.opposite(ft_edge, third_cand);
//                            if (tmpBitset[ft]) {
//                                counters[17]++;
//                                int third_edge_idx = data_.GetEdgeIndex(nxt_cand, third_cand);
//                                int fourth_edge_idx = data_.GetEdgeIndex(cand, ft);
//                                four_cycle_answers[i].push_back({ft_edge, third_edge_idx, fourth_edge_idx});
//                            }
//                        }
//                    }
//                }
////                    fprintf(stderr, "Indexing Cycle %d for edge %d (with edge %d) : %d candidates found\n",i,query_edge_id,opp_edge_idx,four_cycle_answers[i].size());
////                    fprintf(stderr, "Candidates : ");
////                    for (int j = 0; j < four_cycle_answers[i].size(); j++) {
////                        fprintf(stderr, "%d ",four_cycle_answers[i][j]);
////                    }
////                    fprintf(stderr, "\n");
//                for (auto fourth_cand_edge :  data_.GetIncidentEdges(cand, fourth_label)) {
//                    Vertex fourth_cand = data_.opposite(fourth_cand_edge, cand);
//                    tmpBitset[fourth_cand] = false;
//                }
//                if (four_cycle_answers[i].empty()) {
//                    steptimer_3.Stop();
//                    return false;
//                }
//                steptimer_3.Stop();
//            }
//            else {
//                //20-19510
////                    fprintf(stderr, "Cycle %d for edge %d (with edge %d) already indexed : %d candidates found\n",i,query_edge_id,opp_edge_idx,four_cycle_answers[i].size());
//                steptimer_4.Start();
//                while (!four_cycle_answers[i].empty()) {
//                    int candidate_edge_idx, third_edge_idx, fourth_edge_idx;
//                    std::tie(candidate_edge_idx, third_edge_idx, fourth_edge_idx) = four_cycle_answers[i].back();
////                        fprintf(stderr, "    Matching edge %d with %d! ",candidate_edge_idx, opp_edge_idx);
//                    if (EdgeCandidacy(opp_edge_idx, candidate_edge_idx)) {
//                        int third_cand = -1, fourth_cand = -1;
//                        std::tie(third_cand, fourth_cand) = data_.edge_info_[candidate_edge_idx].vp;
////                            fprintf(stderr, "cand 4cycle %d-%d-%d-%d\n",cand,nxt_cand,third_cand,fourth_cand);
//                        if (EdgeCandidacy(query_third_edge_, third_edge_idx)) {
//                            if (EdgeCandidacy(query_fourth_edge_, fourth_edge_idx)) {
//                                steptimer_4.Stop();
//                                counters[15]++;
//                                goto nxt_cycle;
//                            }
//                        }
//                    }
//                    counters[16]++;
//                    four_cycle_answers[i].pop_back();
//                }
////                    fprintf(stderr, "Cycle %d for edge %d (with edge %d) already indexed : %d candidates remains\n",i,query_edge_id,opp_edge_idx,four_cycle_answers[i].size());
//                if (four_cycle_answers[i].empty()) return false;
//                nxt_cycle:
//                steptimer_4.Stop();
//            }
//        }
//        return true;
//    }

    bool CandidateSpace::FourCycleSafety(int query_edge_id, int data_edge_id) {
        steptimer_3.Start();
        int cand, nxt_cand; std::tie(cand, nxt_cand) = data_.edge_info_[data_edge_id].vp;
        int cur, nxt; std::tie(cur, nxt) = query_.edge_info_[query_edge_id].vp;
        std::pair<int, int> edgePair = {query_edge_id, data_edge_id};
        // 4-cycle filter

        if (four_cycle_memo.find(edgePair) == four_cycle_memo.end()) {
            four_cycle_memo[edgePair] = std::vector<cycle_information>(query_.four_cycles[query_edge_id].size());
        }

        std::vector<cycle_information> &four_cycle_answers = four_cycle_memo[edgePair];
        for (int i = 0; i < query_.four_cycles[query_edge_id].size(); i++) {
            int opp_edge_idx = query_.four_cycles[query_edge_id][i];
            auto &info = four_cycle_answers[i];
            bool validity = true;
            steptimer_4.Start();
            if (info.d_opp_edge_idx == -1) {
                validity = false;
//                fprintf(stderr, "Solving : 4-Cycle QEdge [%d]-[%d] and dataedge [%d]-[?] for the first time\n",
//                        query_edge_id, opp_edge_idx, data_edge_id);
                VertexPair opp_edge = query_.edge_info_[opp_edge_idx].vp;
                std::tie(info.third, info.fourth) = opp_edge;
                info.q_opp_edge_idx = opp_edge_idx;
                info.q_fourth_edge_idx = query_.GetEdgeIndex(cur, info.fourth);
                info.q_third_edge_idx = query_.GetEdgeIndex(nxt, info.third);
                info.third_inc_idx = 0;
                info.fourth_inc_idx = 0;
                if (query_.GetEdgeIndex(cur, info.third) != -1) counters[17]++;
                if (query_.GetEdgeIndex(nxt, info.fourth) != -1) counters[17]++;
            }
            else {
                if (validity) validity &= EdgeCandidacy(info.q_opp_edge_idx, info.d_opp_edge_idx);
                if (validity) validity &= EdgeCandidacy(info.q_third_edge_idx, info.d_third_edge_idx);
                if (validity) validity &= EdgeCandidacy(info.q_fourth_edge_idx, info.d_fourth_edge_idx);
            }
            if (validity) {
                counters[16]++;
                steptimer_4.Stop();
                continue;
            }

            steptimer_3.Start();
            Label third_label = query_.GetLabel(info.third);
            Label fourth_label = query_.GetLabel(info.fourth);
            std::vector <int> &third_cands  = data_.GetIncidentEdges(nxt_cand, third_label);
            std::vector <int> &fourth_cands = data_.GetIncidentEdges(cand, fourth_label);
            if (info.third_inc_idx >= third_cands.size()) return false;

            while (info.third_inc_idx < third_cands.size()) {
                counters[15]++;
                bool cand_validity = true;
                int third_cand = data_.opposite(third_cands[info.third_inc_idx], nxt_cand);
                int fourth_cand = data_.opposite(fourth_cands[info.fourth_inc_idx], cand);
                if (cand_validity) cand_validity &= ((fourth_cand != nxt_cand) && (third_cand != cand));
                if (cand_validity) cand_validity &= BitsetCS[info.third][third_cand];
                if (cand_validity) cand_validity &= EdgeCandidacy(info.q_third_edge_idx, third_cands[info.third_inc_idx]);
                if (cand_validity) cand_validity &= BitsetCS[info.fourth][fourth_cand];
                if (cand_validity) cand_validity &= EdgeCandidacy(info.q_fourth_edge_idx, fourth_cands[info.fourth_inc_idx]);
                if (cand_validity) {
                    int d_opp_idx = data_.GetEdgeIndex(third_cand, fourth_cand);
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
                steptimer_3.Stop();
                return false;
            }
            ok:
            steptimer_3.Stop();
            continue;
        }
        return true;
    }

    bool CandidateSpace::EdgeSafety(int query_edge_id, int data_edge_id) {
        // triangle filter
        if (data_.is_sparse() and !query_.triangles[query_edge_id].empty()) {
            if (!TriangleSafety(query_edge_id, data_edge_id)) return false;
        }
#ifdef FOURCYCLE_SAFETY
        if (data_.is_sparse() and !query_.four_cycles[query_edge_id].empty()) {
            if (!FourCycleSafety(query_edge_id, data_edge_id)) return false;
        }
#endif
        return true;
    }

    bool CandidateSpace::Filter(bool topdown) {
        bool result = true;
        bool pruned = false;
        neighbor_label_frequency.resize(data_.GetNumLabels(), 0);
        in_neighbor_cs.resize(data_.GetNumVertices(), false);
        tmpBitset.resize(data_.GetNumVertices(), false);

        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Size query_idx = topdown ? i : query_.GetNumVertices() - i - 1;
            Vertex cur = dag_.GetVertexOrderedByBFS(query_idx);

            if (GetDAGNextCount(cur, topdown) == 0) continue;

            Label cur_label = query_.GetLabel(cur);
            QueryDegree num_nxt = 0;

            for (Size i = 0; i < GetDAGNextCount(cur, topdown); ++i) {
                Vertex nxt = GetDAGNextVertex(cur, i, topdown);
                int query_edge_idx = query_.GetEdgeIndex(nxt, cur);
                for (Vertex nxt_cand : candidate_set_[nxt]) {
                    for (int data_edge_idx : data_.GetIncidentEdges(nxt_cand, cur_label)) {
                        Vertex cand = data_.opposite(data_edge_idx, nxt_cand);
                        if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;
                        if (data_.edge_info_[data_edge_idx].edge_candidacy_[query_edge_idx]==0) {
                            continue;
                        }
                        if (!EdgeSafety(query_edge_idx, data_edge_idx)) {
                            data_.edge_info_[data_edge_idx].edge_candidacy_[query_edge_idx] = false;
                            continue;
                        }
                        if (num_visit_cs_[cand] == num_nxt) {
                            num_visit_cs_[cand] += 1;
                            if (num_nxt == 0) {
                                visited_candidates_[num_visited_candidates] = cand;
                                num_visited_candidates += 1;
                            }
                        }
                    }
                }
                num_nxt += 1;
            }

            PrepareNeighborSafety(cur);


            for (Size i = 0; i < candidate_set_[cur].size(); ++i) {
                Vertex cand = candidate_set_[cur][i];
                bool valid = (num_visit_cs_[cand] == num_nxt);
                if (valid) valid = CheckNeighborSafety(cur, cand);
                if (valid) valid = BipartiteSafety(cur, cand);
                if (!valid) {
                    candidate_set_[cur][i] = candidate_set_[cur].back();
                    candidate_set_[cur].pop_back();
                    --i;
                    BitsetCS[cur].reset(cand);
                    pruned = true;
                }
            }

            if (candidate_set_[cur].empty()) {
                fprintf(stderr, "FOUND INVALID %u\n",cur);
                exit(2);
            }

            while (num_visited_candidates > 0) {
                num_visited_candidates -= 1;
                num_visit_cs_[visited_candidates_[num_visited_candidates]] = 0;
            }
        }

        return pruned;
    }

    void CandidateSpace::ConstructCS() {
        Timer csbuild_timer; csbuild_timer.Start();
        cs_edge_list_.clear();
        tsl::hopscotch_map<Vertex, Size> CandidateSTDSet[query_.GetNumVertices()];
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            for (Size idx = 0; idx < GetCandidateSetSize(i); idx++) {
                CandidateSTDSet[i][candidate_set_[i][idx]] = idx;
                cs_adj_[{i, candidate_set_[i][idx]}] = (tsl::robin_map<Vertex, tsl::robin_set<Vertex>>());
                for (int q_nbr : query_.adj_list[i]) {
                    cs_adj_[{i, candidate_set_[i][idx]}][q_nbr] = tsl::robin_set<Vertex>();
                }
            }
        }
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByBFS(i);
            Label u_label = query_.GetLabel(u);
            Size u_degree = query_.GetDegree(u);

            for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
                Vertex u_adj = query_.GetNeighbor(i);
                int query_edge_idx = query_.GetEdgeIndex(u_adj, u);

                for (Size v_adj_idx = 0; v_adj_idx < candidate_set_[u_adj].size(); ++v_adj_idx) {
                    Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

                    for (int data_edge_idx : data_.GetIncidentEdges(v_adj, u_label)) {
                        Vertex v = data_.opposite(data_edge_idx, v_adj);

                        if (data_.GetDegree(v) < u_degree) break;
                        if (!EdgeCandidacy(query_edge_idx, data_edge_idx)) continue;
                        if (BitsetCS[u][v]){
                            Size v_idx = CandidateSTDSet[u][v];
                            VertexPair u_pair = {u, v_idx};
                            VertexPair uc_pair = {u_adj, v_adj_idx};
                            cs_edge_list_[u_pair].insert(uc_pair);
                            cs_edge_list_[uc_pair].insert(u_pair);

                            /* Add CS_ADJ_List; */
                            cs_adj_[{u, v}][u_adj].insert(v_adj);
                            cs_adj_[{u_adj, v_adj}][u].insert(v);
                        }
                    }
                }
            }
        }
        Size CandidateSpaceSize = 0, CandidateEdges = 0;
        for (Vertex i = 0; i < query_.GetNumVertices(); ++i) {
            CandidateSpaceSize += CandidateSTDSet[i].size();
        }
        for (auto it : cs_edge_list_) {
            CandidateEdges += it.second.size();
        }
        CandidateEdges /= 2;
        csbuild_timer.Stop();

        fprintf(stdout, "Counter 15 = %d, 16 = %d, 17 = %d\n",counters[15], counters[16], counters[17]);
        std::cout << "GetEdgeIndex Calls : " << functionCallCounter << std::endl;
        std::cout << "StepTimer 1 (BP-Safety) " << steptimer_1.GetTime() << "ms" << std::endl;
        std::cout << "StepTimer 2 (Triangle-Safety) " << steptimer_2.GetTime() << "ms" << std::endl;
        std::cout << "StepTimer 3 (4Cycle-Safety) " << steptimer_3.GetTime() << "ms" << std::endl;
        std::cout << "StepTimer 4 (4Cycle-Safety-Bypass) " << steptimer_4.GetTime() << "ms" << std::endl;
        fprintf(stdout, "Total 4-cycle checks = %u, Passed %u\n",counters[0], counters[1]);
        fprintf(stdout, "CS build time : %.02lf ms\n",csbuild_timer.GetTime());
        fprintf(stderr, "CS-SZ : %u %u\n", CandidateSpaceSize, CandidateEdges);
        fprintf(stdout, "#CandidateSetSize : %u\n", CandidateSpaceSize);
        fprintf(stdout, "#CandidateSetEdges : %u\n", CandidateEdges);
        double vert_cand_avg = (1.0*CandidateSpaceSize)/query_.GetNumVertices();
        fprintf(stdout, "#AVG_VERT_CAND_SIZE : %.02lf\n", vert_cand_avg);
        fprintf(stdout, "#AVG_EDGE_BTW_CANDS : %.02lf\n", (1.0*CandidateEdges)/query_.GetNumEdges()/vert_cand_avg);
    }

    bool CandidateSpace::InitRootCandidates() {
        Vertex root = dag_.GetRoot();
        Label root_label = query_.GetLabel(root);

        uint64_t *nbr_label_bitset = new uint64_t[data_.GetNbrBitsetSize()];
        Size max_nbr_degree;

        ComputeNbrInformation(root, &max_nbr_degree, nbr_label_bitset);

        for (Size i = data_.GetStartOffsetByLabel(root_label);
             i < data_.GetEndOffsetByLabel(root_label); ++i) {
            Vertex cand = data_.GetVertexBySortedLabelOffset(i);

            if (data_.GetDegree(cand) < query_.GetDegree(root)) continue;

            if (data_.GetCoreNum(cand) >= query_.GetCoreNum(root) &&
                data_.CheckAllNbrLabelExist(cand, nbr_label_bitset) &&
                data_.GetMaxNbrDegree(cand) >= max_nbr_degree) {
                candidate_set_[root].push_back(cand);
                BitsetCS[root].set(cand);
            }
        }

        delete[] nbr_label_bitset;
        return !candidate_set_[root].empty();
    }

    void CandidateSpace::ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                               uint64_t *nbr_label_bitset) {
        *max_nbr_degree = 0;
        std::fill(nbr_label_bitset, nbr_label_bitset + data_.GetNbrBitsetSize(),
                  0ull);
        for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
            Vertex adj = query_.GetNeighbor(i);

            nbr_label_bitset[query_.GetLabel(adj) / (sizeof(uint64_t) * CHAR_BIT)] |=
                    1ull << (query_.GetLabel(adj) % (sizeof(uint64_t) * CHAR_BIT));

            if (query_.GetDegree(adj) > *max_nbr_degree) {
                *max_nbr_degree = query_.GetDegree(adj);
            }
        }
    }

    void CandidateSpace::printCS() {
        for (int i = 0; i < query_.GetNumVertices(); i++) {
            fprintf(stdout, "Query [%d] : ",i);
            for (int x : candidate_set_[i]) {
                fprintf(stdout, "%d ", x);
            }
            fprintf(stdout, "\n");
            for (int x : candidate_set_[i]) {
                fprintf(stdout, "  CS Edges from %d\n",x);
                for (auto it : cs_adj_[{i, x}]) {
                    fprintf(stdout, "    Edge with [%d] : ",it.first);
                    for (auto itt: it.second) {
                        fprintf(stdout, "%d ", itt);
                    }
                    fprintf(stdout, "\n");
                }
            }
        }
    }
}
