#include <vector>
#include <map>
#include <random>
#include "include/daf_candidate_space.h"
//#define BIPARTITE_SAFETY
#define NEIGHBOR_SAFETY
//#define FOURCYCLE_SAFETY
#ifdef NEIGHBOR_SAFETY
    #undef BIPARTITE_SAFETY
#endif
#ifdef BIPARTITE_SAFETY
    #undef NEIGHBOR_SAFETY
#endif


namespace daf {
    CandidateSpace::CandidateSpace(DataGraph &data, QueryGraph &query,
                                   DAG &dag)
            : data_(data), query_(query), dag_(dag) {
        candidate_set_size_ = new Size[query_.GetNumVertices()];
        candidate_set_ = new Vertex *[query_.GetNumVertices()];
        candidate_offsets_ = new Size **[query_.GetNumVertices()];
        BitsetCS.resize(query_.GetNumVertices());

        for (Vertex u = 0; u < query_.GetNumVertices(); ++u) {
            BitsetCS[u].resize(data_.GetNumVertices());
            if (query_.IsInNEC(u) && !query_.IsNECRepresentation(u)) {
                candidate_set_[u] = nullptr;
                candidate_offsets_[u] = nullptr;
            }
            else {
                candidate_set_[u] = new Vertex[dag_.GetInitCandSize(u)];
                candidate_offsets_[u] = new Size *[query_.GetDegree(u)];
            }
        }

        linear_cs_adj_list_ = nullptr;

        num_visit_cs_ = new QueryDegree[data_.GetNumVertices()];
        visited_candidates_ =
                new Vertex[data_.GetNumVertices() * query_.GetNumVertices()];
        cand_to_cs_idx_ = new Size[data_.GetNumVertices()];

        num_visitied_candidates_ = 0;

        std::fill(num_visit_cs_, num_visit_cs_ + data_.GetNumVertices(), 0);
        std::fill(candidate_set_size_, candidate_set_size_ + query_.GetNumVertices(),
                  0);
        std::fill(cand_to_cs_idx_, cand_to_cs_idx_ + data_.GetNumVertices(),
                  INVALID_SZ);
        num_cs_edges_ = 0;
    }

    CandidateSpace::~CandidateSpace() {
        delete[] candidate_set_size_;

        for (Vertex u = 0; u < query_.GetNumVertices(); ++u) {
            if (candidate_set_[u] != nullptr) {
                delete[] candidate_set_[u];
                for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
                    Vertex u_adj = query_.GetNeighbor(i);
                    if (candidate_offsets_[u][i - query_.GetStartOffset(u)] != nullptr) {
//                        fprintf(stderr, "Deleting candidate_offsets_[%d][%d]\n", i, i - query_.GetStartOffset(u));
                        delete[] candidate_offsets_[u][i - query_.GetStartOffset(u)];
                    }

                }
                delete[] candidate_offsets_[u];
            }
        }
        delete[] candidate_set_;
        delete[] candidate_offsets_;

        if (linear_cs_adj_list_ != nullptr) delete[] linear_cs_adj_list_;
        delete[] num_visit_cs_;
        delete[] visited_candidates_;
        delete[] cand_to_cs_idx_;
    }


    bool CandidateSpace::BuildCS() {
        srand(0);
        dag_.BuildDAG(-1);
        if (!FilterByTopDownWithInit()) return false;
        Size stable_iter = 0;
        while (stable_iter < 3) {
            stable_iter++;
            bool pruned = false;
            if (Filter(false)) pruned = true;
            if (Filter(true)) pruned = true;
            if (!pruned) break;
        }
        ConstructCS();
        return true;
    }


    bool CandidateSpace::FilterByTopDownWithInit() {
        bool result = true;

        uint64_t *nbr_label_bitset = new uint64_t[data_.GetNbrBitsetSize()];
        Size max_nbr_degree;

        if (!InitRootCandidates()) {
            result = false;
        }
        else {
            for (Size i = 1; i < query_.GetNumVertices(); ++i) {
                Vertex cur = dag_.GetVertexOrderedByBFS(i);

                if (query_.IsInNEC(cur) && !query_.IsNECRepresentation(cur)) continue;

                Label cur_label = query_.GetLabel(cur);
                QueryDegree num_parent = 0;
                for (Size i = 0; i < dag_.GetNumParents(cur); ++i) {
                    Vertex parent = dag_.GetParent(cur, i);

                    for (Size i = 0; i < candidate_set_size_[parent]; ++i) {
                        Vertex parent_cand = candidate_set_[parent][i];

                        for (Size i = data_.GetStartOffset(parent_cand, cur_label);
                             i < data_.GetEndOffset(parent_cand, cur_label); ++i) {
                            Vertex cand = data_.GetNeighbor(i);

                            if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;

                            if (num_visit_cs_[cand] == num_parent) {
                                num_visit_cs_[cand] += 1;
                                if (num_parent == 0) {
                                    visited_candidates_[num_visitied_candidates_] = cand;
                                    num_visitied_candidates_ += 1;
                                }
                            }
                        }
                    }
                    num_parent += 1;
                }

                ComputeNbrInformation(cur, &max_nbr_degree, nbr_label_bitset);

                for (Size i = 0; i < num_visitied_candidates_; ++i) {
                    Vertex cand = visited_candidates_[i];
                    if (num_visit_cs_[cand] == num_parent &&
                        data_.GetCoreNum(cand) >= query_.GetCoreNum(cur) &&
                        data_.GetMaxNbrDegree(cand) >= max_nbr_degree &&
                        data_.CheckAllNbrLabelExist(cand, nbr_label_bitset)) {
                        candidate_set_[cur][candidate_set_size_[cur]] = cand;
                        candidate_set_size_[cur] += 1;
                        BitsetCS[cur].set(cand);
                    }
                }

                if (candidate_set_size_[cur] == 0) {
                    result = false;
                    break;
                }

                while (num_visitied_candidates_ > 0) {
                    num_visitied_candidates_ -= 1;
                    num_visit_cs_[visited_candidates_[num_visitied_candidates_]] = 0;
                }
            }
        }

        delete[] nbr_label_bitset;
        return result;
    }

    bool CandidateSpace::BipartiteSafety(Vertex cur, Vertex cand) {
#ifdef BIPARTITE_SAFETY
        BipartiteMaximumMatching BP(query_.GetDegree(cur)+1, data_.GetDegree(cand)+1);
        BP.clear_adj();
        BP.add_edge(query_.GetDegree(cur), data_.GetDegree(cand));
        for (int i = 0; i < query_.adj_list[cur].size(); i++) {
            Vertex uc = query_.adj_list[cur][i];
            for (int j = 0; j < data_.adj_list[cand].size(); j++) {
                Vertex vc = data_.adj_list[cand][j];
                if (BitsetCS[uc][vc])
                    BP.add_edge(i, j);
            }
        }
        return (BP.solve() == (query_.GetDegree(cur)+1));
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
            for (Size ncsidx = 0; ncsidx < candidate_set_size_[q_neighbor]; ++ncsidx) {
                Vertex d_neighbor = candidate_set_[q_neighbor][ncsidx];
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

    bool CandidateSpace::EdgeSafety(Vertex cur, Vertex cand, Vertex nxt, Vertex nxt_cand) {
        if (query_.GetEdgeLabel(cur, nxt) != data_.GetEdgeLabel(cand, nxt_cand)) return false;
        std::vector<Vertex> datagraph_cand_intersection;
        std::set_intersection(data_.adj_list[cand].begin(), data_.adj_list[cand].end(),
                              data_.adj_list[nxt_cand].begin(), data_.adj_list[nxt_cand].end(),
                              std::back_inserter(datagraph_cand_intersection));
        if (query_.triangles[cur][nxt].size() > datagraph_cand_intersection.size()) {
            return false;
        }
#ifdef BIPARTITE_SAFETY
        BipartiteMaximumMatching BP(query_.triangles[cur][nxt].size(), datagraph_cand_intersection.size());
        BP.clear_adj();
        for (int i = 0; i < query_.triangles[cur][nxt].size(); i++) {
            Vertex uu = query_.triangles[cur][nxt][i];
            for (int j = 0; j < datagraph_cand_intersection.size(); j++) {
                Vertex vv = datagraph_cand_intersection[j];
                if (BitsetCS[uu][vv]) {
                    BP.add_edge(i, j);
                }
            }
        }
        return (BP.solve() == query_.triangles[cur][nxt].size());
#else
        // triangle filter
        for (auto trig_vertex_ : query_.triangles[cur][nxt]) {
            bool found = false;
            for (Vertex &tcand : datagraph_cand_intersection) {
                if (BitsetCS[trig_vertex_][tcand]) {
                    found = true;
                    break;
                }
            }
            if (!found) return false;
        }
#endif
#ifdef FOURCYCLE_SAFETY
        // 4-cycle filter
        VertexPair opp_edge = query_.four_cycles[cur][nxt];
        if (opp_edge.first != opp_edge.second) {
            Vertex third = opp_edge.first, fourth = opp_edge.second;
            for (auto fourth_cand : data_.adj_list[cand]) {
                if (BitsetCS[fourth][fourth_cand] and fourth_cand != nxt_cand)
                    tmpBitset[fourth_cand] = true;
            }
            std::vector <Vertex> third_cands;
            for (auto third_cand : data_.adj_list[nxt_cand]) {
                if (BitsetCS[third][third_cand] and third_cand != cand) {
                    third_cands.push_back(third_cand);
                }
            }
            bool found_valid_cycle = false;
            for (auto t : third_cands) {
                for (auto ft : data_.adj_list[t]) {
                    if (tmpBitset[ft]) {
                        found_valid_cycle = true;
                        goto checked;
                    }
                }
            }
            checked:
            for (auto fourth_cand : data_.adj_list[cand]) {
                tmpBitset[fourth_cand] = false;
            }
            if (!found_valid_cycle) {
                return false;
            }
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
                for (Size i = 0; i < candidate_set_size_[nxt]; ++i) {
                    Vertex nxt_cand = candidate_set_[nxt][i];
                    for (Size i = data_.GetStartOffset(nxt_cand, cur_label);
                         i < data_.GetEndOffset(nxt_cand, cur_label); ++i) {
                        Vertex cand = data_.GetNeighbor(i);
                        if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;
                        if (!EdgeSafety(cur, cand, nxt, nxt_cand)) continue;
                        if (num_visit_cs_[cand] == num_nxt) {
                            num_visit_cs_[cand] += 1;
                            if (num_nxt == 0) {
                                visited_candidates_[num_visitied_candidates_] = cand;
                                num_visitied_candidates_ += 1;
                            }
                        }
                    }
                }
                num_nxt += 1;
            }

            PrepareNeighborSafety(cur);


            for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
                Vertex cand = candidate_set_[cur][i];
                bool valid = (num_visit_cs_[cand] == num_nxt);
                if (valid) valid = CheckNeighborSafety(cur, cand);
                if (valid) valid = BipartiteSafety(cur, cand);
                if (!valid) {
                    candidate_set_[cur][i] =
                            candidate_set_[cur][candidate_set_size_[cur] - 1];
                    candidate_set_size_[cur] -= 1;
                    --i;
                    BitsetCS[cur].reset(cand);
                    pruned = true;
                }
                else {
                    num_cs_edges_ += cand_to_cs_idx_[cand];
                }
            }

            if (candidate_set_size_[cur] == 0) {
                fprintf(stderr, "FOUND INVALID %u\n",cur);
                exit(2);
            }

            while (num_visitied_candidates_ > 0) {
                num_visitied_candidates_ -= 1;
                num_visit_cs_[visited_candidates_[num_visitied_candidates_]] = 0;
            }
        }

        return pruned;
    }

    void CandidateSpace::ConstructCS() {
        std::map<Vertex, Size> CandidateSTDSet[query_.GetNumVertices()];
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            for (Size idx = 0; idx < GetCandidateSetSize(i); idx++) {
                CandidateSTDSet[i][candidate_set_[i][idx]] = idx;
            }
        }
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByBFS(i);
            Label u_label = query_.GetLabel(u);
            Size u_degree = query_.GetDegree(u);

            for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
                Vertex u_adj = query_.GetNeighbor(i);
                Size u_adj_idx = i - query_.GetStartOffset(u);

                for (Size v_adj_idx = 0; v_adj_idx < candidate_set_size_[u_adj]; ++v_adj_idx) {
                    Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

                    for (Size i = data_.GetStartOffset(v_adj, u_label);
                         i < data_.GetEndOffset(v_adj, u_label); ++i) {
                        Vertex v = data_.GetNeighbor(i);

                        if (data_.GetDegree(v) < u_degree) break;
                        if (!EdgeSafety(u, v, u_adj, v_adj)) continue;
                        if (CandidateSTDSet[u].find(v) != CandidateSTDSet[u].end()){
                            Size v_idx = CandidateSTDSet[u][v];
                            VertexPair u_pair = {u, v_idx};
                            VertexPair uc_pair = {u_adj, v_adj_idx};
                            cs_edge_list_[u_pair].insert(uc_pair);
                            cs_edge_list_[uc_pair].insert(u_pair);
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
        Size chk = 0;
        for (auto &it : BitsetCS) {
            chk += it.count();
        }
        int outcsDataVertex = 0;
        for (int i = 0; i < data_.GetNumVertices(); i++) {
            bool found = false;
            for (int j = 0; j < query_.GetNumVertices(); j++) {
                if (BitsetCS[j][i]) {
                    found = true;
                    break;
                }
            }
            if (!found) {
                outcsDataVertex++;
            }
        }
        fprintf(stdout, "#BitCSSize : %u\n", chk);
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

            if (data_.GetDegree(cand) < query_.GetDegree(root)) break;

            if (data_.GetCoreNum(cand) >= query_.GetCoreNum(root) &&
                data_.CheckAllNbrLabelExist(cand, nbr_label_bitset) &&
                data_.GetMaxNbrDegree(cand) >= max_nbr_degree) {
                candidate_set_[root][candidate_set_size_[root]] = cand;
                candidate_set_size_[root] += 1;
                BitsetCS[root].set(cand);
            }
        }

        delete[] nbr_label_bitset;
        return candidate_set_size_[root] > 0;
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
}
