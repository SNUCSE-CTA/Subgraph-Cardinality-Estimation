#include <vector>
#include <map>
#include <random>
#include "include/daf_candidate_space.h"
//#define BIPARTITE_SAFETY

namespace daf {
    CandidateSpace::CandidateSpace(const DataGraph &data, const QueryGraph &query,
                                   DAG &dag)
            : data_(data), query_(query), dag_(dag) {
        candidate_set_size_ = new Size[query_.GetNumVertices()];
        candidate_set_ = new Vertex *[query_.GetNumVertices()];
        candidate_offsets_ = new Size **[query_.GetNumVertices()];

        for (Vertex u = 0; u < query_.GetNumVertices(); ++u) {
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
        auto CandidateSpaceSize = [this]() {
            Size cand_size_ = 0;
            for (Vertex i = 0; i < query_.GetNumVertices(); ++i) {
                cand_size_ += candidate_set_size_[i];
            }
            return cand_size_;
        };
        auto CandidateEdges = [this]() {
            Size num_edges_ = 0;
            for (auto &it : cs_edge_list_) {
                num_edges_ += it.second.size();
            }
            return num_edges_;
        };
        srand(0);
        dag_.BuildDAG(-1);
        if (!FilterByTopDownWithInit()) return false;

//        fprintf(stderr, "candidate_set_size = %u\n", CandidateSpaceSize());

        Size stable_iter = 0;
        while (true) {
            stable_iter++;
            bool pruned = false;
            if (FilterByBottomUp()) pruned = true;
            if (FilterByTopDown()) pruned = true;
//            fprintf(stderr, "Iteration [%u] candidate_set_size = %u\n",stable_iter, CandidateSpaceSize());
            if (!pruned) break;
        }
        ConstructCS();
        fprintf(stdout, "#CandidateSetSize : %u\n", CandidateSpaceSize());
        fprintf(stdout, "#CandidateSetEdges : %u\n", CandidateEdges());

//        fprintf(stderr, "    Candidate_set_size = %u, Edge %u\n", CandidateSpaceSize(),CandidateEdges());
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

//#define BIPARTITE_SAFETY
    bool CandidateSpace::BipartiteSafety(Vertex cur, Vertex cand) {
#ifdef BIPARTITE_SAFETY
        BipartiteMaximumMatching BP(query_.GetNumVertices(), data_.GetNumVertices());
        BP.clear_adj();
        BP.add_edge(cur, cand);
        for (Vertex ui = 0; ui < query_.GetNumVertices(); ++ui) {
            if (ui == cur) continue;
            bool is_nb = query_.CheckEdgeExist(ui, cur);
            for (Size lvcandidx = 0; lvcandidx < candidate_set_size_[ui]; ++lvcandidx) {
                Vertex lv_cand = candidate_set_[ui][lvcandidx];
                if (lv_cand == cand) continue;
                if (is_nb && !data_.CheckEdgeExist(lv_cand, cand)) continue;
                BP.add_edge(ui, lv_cand);
            }
        }
        return (BP.solve() == query_.GetNumVertices());
#else
        return true;
#endif
    }

    bool CandidateSpace::FilterByBottomUp() {
        bool result = true;
        bool pruned = false;
        int *neighbor_label_frequency = new int[data_.GetNumLabels()];
        int *tmp_label_frequency = new int[data_.GetNumLabels()];
        bool *in_neighbor_cs = new bool[data_.GetNumVertices()];
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex cur = dag_.GetVertexOrderedByBFS(query_.GetNumVertices() - i - 1);

            if (dag_.GetNumChildren(cur) == 0) continue;

            Label cur_label = query_.GetLabel(cur);

            QueryDegree num_child = 0;
            for (Size i = 0; i < dag_.GetNumChildren(cur); ++i) {
                Vertex child = dag_.GetChild(cur, i);

                if (query_.IsInNEC(child) && !query_.IsNECRepresentation(child)) continue;

                for (Size i = 0; i < candidate_set_size_[child]; ++i) {
                    Vertex child_cand = candidate_set_[child][i];

                    for (Size i = data_.GetStartOffset(child_cand, cur_label);
                         i < data_.GetEndOffset(child_cand, cur_label); ++i) {
                        Vertex cand = data_.GetNeighbor(i);
                        if (data_.GetDegree(cand) < query_.GetDegree(cur)) break;

                        if (num_visit_cs_[cand] == num_child) {
                            num_visit_cs_[cand] += 1;
                            if (num_child == 0) {
                                visited_candidates_[num_visitied_candidates_] = cand;
                                num_visitied_candidates_ += 1;
                            }
                        }
                    }
                }
                num_child += 1;
            }

            memset(neighbor_label_frequency, 0, sizeof(Size) * data_.num_label_);
            memset(in_neighbor_cs, 0, sizeof(bool) * data_.GetNumVertices());
            for (Size nidx = query_.GetStartOffset(cur); nidx < query_.GetEndOffset(cur); ++nidx) {
                Vertex q_neighbor = query_.GetNeighbor(nidx);
                neighbor_label_frequency[query_.GetLabel(q_neighbor)]++;
                for (Size ncsidx = 0; ncsidx < candidate_set_size_[q_neighbor]; ++ncsidx) {
                    Vertex d_neighbor = candidate_set_[q_neighbor][ncsidx];
                    in_neighbor_cs[d_neighbor] = true;
                }
            }

            for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
                Vertex cand = candidate_set_[cur][i];
                bool valid = (num_visit_cs_[cand] == num_child);
                if (valid) {
                    memcpy(tmp_label_frequency, neighbor_label_frequency, sizeof(int) * data_.GetNumLabels());
                    // Check neighbor safety
                    for (Size nidx = data_.GetStartOffset(cand); nidx < data_.GetEndOffset(cand); nidx++) {
                        Vertex d_neighbor = data_.GetNeighbor(nidx);
                        if (in_neighbor_cs[d_neighbor]) {
                            tmp_label_frequency[data_.GetLabel(d_neighbor)]--;
                        }
                    }
                    for (int l = 0; l < data_.GetNumLabels(); ++l) {
                        if (tmp_label_frequency[l] > 0) {
                            valid = false;
                            break;
                        }
                    }
                }
                if (valid) {
                    valid = BipartiteSafety(cur, cand);
                }

                if (!valid) {
                    candidate_set_[cur][i] =
                            candidate_set_[cur][candidate_set_size_[cur] - 1];
                    candidate_set_size_[cur] -= 1;
                    --i;
                    pruned = true;
                }
                else {
                    num_cs_edges_ += cand_to_cs_idx_[cand];
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
        delete[] neighbor_label_frequency;
        delete[] tmp_label_frequency;
        delete[] in_neighbor_cs;

        return pruned;
    }

    bool CandidateSpace::FilterByTopDown() {
        bool result = true;
        bool pruned = false;
        int *neighbor_label_frequency = new int[data_.GetNumLabels()];
        int *tmp_label_frequency = new int[data_.GetNumLabels()];
        bool *in_neighbor_cs = new bool[data_.GetNumVertices()];

        for (Size i = 1; i < query_.GetNumVertices(); ++i) {
            Vertex cur = dag_.GetVertexOrderedByBFS(i);
            Label cur_label = query_.GetLabel(cur);

            // ComputeCounter()
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
                                cand_to_cs_idx_[cand] = 0;
                            }
                            cand_to_cs_idx_[cand] += 1;
                        }
                        else if (num_visit_cs_[cand] > num_parent) {
                            cand_to_cs_idx_[cand] += 1;
                        }
                    }
                }
                num_parent += 1;
            }


            memset(neighbor_label_frequency, 0, sizeof(Size) * data_.num_label_);
            memset(in_neighbor_cs, 0, sizeof(bool) * data_.GetNumVertices());
            for (Size nidx = query_.GetStartOffset(cur); nidx < query_.GetEndOffset(cur); ++nidx) {
                Vertex q_neighbor = query_.GetNeighbor(nidx);
                neighbor_label_frequency[query_.GetLabel(q_neighbor)]++;
                for (Size ncsidx = 0; ncsidx < candidate_set_size_[q_neighbor]; ++ncsidx) {
                    Vertex d_neighbor = candidate_set_[q_neighbor][ncsidx];
                    in_neighbor_cs[d_neighbor] = true;
                }
            }

            std::vector<Vertex> sameLabelQueryVertices = query_.verticesbyLabel[query_.GetLabel(cur)];
            for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
                Vertex cand = candidate_set_[cur][i];

                bool valid = (num_visit_cs_[cand] == num_parent);
                if (valid) {
                    memcpy(tmp_label_frequency, neighbor_label_frequency, sizeof(int) * data_.GetNumLabels());
                    for (Size nidx = data_.GetStartOffset(cand); nidx < data_.GetEndOffset(cand); nidx++) {
                        Vertex d_neighbor = data_.GetNeighbor(nidx);
                        if (in_neighbor_cs[d_neighbor]) {
                            tmp_label_frequency[data_.GetLabel(d_neighbor)]--;
                        }
                    }
                    for (int l = 0; l < data_.GetNumLabels(); ++l) {
                        if (tmp_label_frequency[l] > 0) {
                            valid = false;
                            break;
                        }
                    }
                }
                if (valid) valid = BipartiteSafety(cur, cand);

                if (!valid) {
                    candidate_set_[cur][i] =
                            candidate_set_[cur][candidate_set_size_[cur] - 1];
                    candidate_set_size_[cur] -= 1;
                    --i;
                    pruned = true;
                }
                else {
                    num_cs_edges_ += cand_to_cs_idx_[cand];
                }
            }

            if (candidate_set_size_[cur] == 0) {
                result = false;
                break;
            }

            while (num_visitied_candidates_ > 0) {
                num_visitied_candidates_ -= 1;
                num_visit_cs_[visited_candidates_[num_visitied_candidates_]] = 0;
                cand_to_cs_idx_[visited_candidates_[num_visitied_candidates_]] = INVALID_SZ;
            }
        }
        delete[] neighbor_label_frequency;
        delete[] tmp_label_frequency;
        delete[] in_neighbor_cs;
        return pruned;
    }

    void CandidateSpace::ConstructCS() {
        linear_cs_adj_list_ = new Vertex[num_cs_edges_ * 2];

        Size cur_cand_offset = 0;

        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByBFS(i);

            if (query_.IsInNEC(u) && !query_.IsNECRepresentation(u)) continue;

            Label u_label = query_.GetLabel(u);
            Size u_degree = query_.GetDegree(u);

            for (Size i = 0; i < candidate_set_size_[u]; ++i) {
                cand_to_cs_idx_[candidate_set_[u][i]] = i;
            }

            for (Size i = query_.GetStartOffset(u); i < query_.GetEndOffset(u); ++i) {
                Vertex u_adj = query_.GetNeighbor(i);

                if (query_.IsInNEC(u_adj) && !query_.IsNECRepresentation(u_adj)) continue;

                Size u_adj_idx = i - query_.GetStartOffset(u);

                candidate_offsets_[u][u_adj_idx] = new Size[GetCandidateSetSize(u) + 1];

                std::fill(candidate_offsets_[u][u_adj_idx],
                          candidate_offsets_[u][u_adj_idx] + GetCandidateSetSize(u) + 1,
                          cur_cand_offset);

                for (Size v_adj_idx = 0; v_adj_idx < candidate_set_size_[u_adj]; ++v_adj_idx) {
                    Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

                    for (Size i = data_.GetStartOffset(v_adj, u_label);
                         i < data_.GetEndOffset(v_adj, u_label); ++i) {
                        Vertex v = data_.GetNeighbor(i);

                        if (data_.GetDegree(v) < u_degree) break;

                        if (cand_to_cs_idx_[v] != INVALID_SZ) {
                            candidate_offsets_[u][u_adj_idx][cand_to_cs_idx_[v] + 1] += 1;
                        }
                    }
                }

                for (Size i = 2; i < GetCandidateSetSize(u) + 1; ++i) {
                    candidate_offsets_[u][u_adj_idx][i] +=
                            candidate_offsets_[u][u_adj_idx][i - 1] - cur_cand_offset;
                }

                for (Size v_adj_idx = 0; v_adj_idx < candidate_set_size_[u_adj];
                     ++v_adj_idx) {
                    Vertex v_adj = candidate_set_[u_adj][v_adj_idx];

                    for (Size i = data_.GetStartOffset(v_adj, u_label);
                         i < data_.GetEndOffset(v_adj, u_label); ++i) {
                        Vertex v = data_.GetNeighbor(i);

                        if (data_.GetDegree(v) < u_degree) break;

                        Size v_idx = cand_to_cs_idx_[v];
                        if (v_idx != INVALID_SZ) {
                            Vertex u_cand = GetCandidate(u, v_idx);
                            Vertex uc_cand = GetCandidate(u_adj, v_adj_idx);
                            VertexPair u_pair = {u, v_idx};
                            VertexPair uc_pair = {u_adj, v_adj_idx};
//                            fprintf(stderr, "ADDEDGE (%u, %u) <-> (%u, %u)\n",u_pair.first,u_pair.second,uc_pair.first,uc_pair.second);
                            cs_edge_list_[u_pair].insert(uc_pair);
                            cs_edge_list_[uc_pair].insert(u_pair);
                            linear_cs_adj_list_[candidate_offsets_[u][u_adj_idx][v_idx]] =
                                    v_adj_idx;
                            candidate_offsets_[u][u_adj_idx][v_idx] += 1;
                        }
                    }
                }

                for (Size i = GetCandidateSetSize(u) - 1; i--;) {
                    candidate_offsets_[u][u_adj_idx][i + 1] =
                            candidate_offsets_[u][u_adj_idx][i];
                }

                candidate_offsets_[u][u_adj_idx][0] = cur_cand_offset;

                cur_cand_offset =
                        candidate_offsets_[u][u_adj_idx][GetCandidateSetSize(u)];
            }

            for (Size i = 0; i < candidate_set_size_[u]; ++i) {
                cand_to_cs_idx_[candidate_set_[u][i]] = INVALID_SZ;
            }
        }
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
