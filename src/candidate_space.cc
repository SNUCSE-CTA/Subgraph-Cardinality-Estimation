#include <vector>
#include <map>
#include <random>
#include "include/candidate_space.h"

namespace daf {
    CandidateSpace::CandidateSpace(const DataGraph &data, const QueryGraph &query,
                                   DAG &dag)
            : data_(data), query_(query), dag_(dag) {
        candidate_set_size_ = new Size[query_.GetNumVertices()];
        candidate_set_ = new Vertex *[query_.GetNumVertices()];
        num_trees_ = new double *[query_.GetNumVertices()];
        candidate_offsets_ = new Size **[query_.GetNumVertices()];

        for (Vertex u = 0; u < query_.GetNumVertices(); ++u) {
            if (query_.IsInNEC(u) && !query_.IsNECRepresentation(u)) {
                candidate_set_[u] = nullptr;
                num_trees_[u] = nullptr;
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
        delete[] candidate_weights_;
        delete[] num_visit_cs_;
        delete[] visited_candidates_;
        delete[] cand_to_cs_idx_;
    }

    bool CandidateSpace::BuildCS() {
        srand(0);
        dag_.BuildDAG(-1);
        if (!FilterByTopDownWithInit()) return false;
        if (!FilterByBottomUp()) return false;
        if (!FilterByTopDown()) return false;
        for (Size iteration = 0; iteration < 5; iteration++) {
            dag_.BuildDAG(rand()%query_.GetNumVertices());
            if (!FilterByTopDown()) return false;
            if (!FilterByBottomUp()) return false;
            if (!FilterByTopDown()) return false;
        }
        dag_.BuildDAG(-1);
        if (!FilterByTopDown()) return false;
        if (!FilterByBottomUp()) return false;
        if (!FilterByTopDown()) return false;

        ConstructCS();
        BuildQueryTree();
        ConstructTreeDP();

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

    bool CandidateSpace::FilterByBottomUp() {
        bool result = true;
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

            for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
                Vertex cand = candidate_set_[cur][i];
                if (num_visit_cs_[cand] != num_child) {
                    candidate_set_[cur][i] =
                            candidate_set_[cur][candidate_set_size_[cur] - 1];
                    candidate_set_size_[cur] -= 1;
                    --i;
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

        return result;
    }

    bool CandidateSpace::FilterByTopDown() {
        bool result = true;
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

            for (Size i = 0; i < candidate_set_size_[cur]; ++i) {
                Vertex cand = candidate_set_[cur][i];

                memcpy(tmp_label_frequency, neighbor_label_frequency, sizeof(int) * data_.GetNumLabels());
                // Check neighbor safety
                bool is_neighbor_safe = true;
                for (Size nidx = data_.GetStartOffset(cand); nidx < data_.GetEndOffset(cand); nidx++) {
                    Vertex d_neighbor = data_.GetNeighbor(nidx);
                    if (in_neighbor_cs[d_neighbor]) {
                        tmp_label_frequency[data_.GetLabel(d_neighbor)]--;
                    }
                }
                for (int l = 0; l < data_.GetNumLabels(); ++l) {
                    if (tmp_label_frequency[l] > 0) {
                        is_neighbor_safe = false;
                        break;
                    }
                }
                if ((num_visit_cs_[cand] != num_parent) || (is_neighbor_safe == false)) {
                    candidate_set_[cur][i] =
                            candidate_set_[cur][candidate_set_size_[cur] - 1];
                    candidate_set_size_[cur] -= 1;
                    --i;
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
        return result;
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

    void CandidateSpace::BuildQueryTree() {

        std::vector<std::pair<double, std::pair<Vertex, Vertex>>> edges;
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            for (Size nidx = query_.GetStartOffset(i); nidx < query_.GetEndOffset(i); ++nidx) {
                Vertex q_neighbor = query_.GetNeighbor(nidx);
                if (i > q_neighbor) continue;
                Size ij_cs_edge = 0;
                for (Size cs_idx = 0; cs_idx < candidate_set_size_[i]; ++cs_idx) {
                    Size num_cs_neighbor =
                            GetCandidateEndOffset(i, nidx-query_.GetStartOffset(i), cs_idx) -
                            GetCandidateStartOffset(i, nidx-query_.GetStartOffset(i), cs_idx);
                    ij_cs_edge += num_cs_neighbor;
                }
                ij_cs_edge /= ((1 + candidate_set_size_[i]) * (1 + candidate_set_size_[q_neighbor]));
                edges.push_back({ij_cs_edge * 1.0,{i, q_neighbor}});
            }
        }
        std::sort(edges.begin(), edges.end());
        UnionFind uf(query_.GetNumVertices());
        for (auto e : edges) {
            auto [u, v] = e.second;
            if (uf.unite(u, v)) {
//                fprintf(stderr, "Tree %d-%d %.02lf\n",u,v,e.first);
                dag_.AddTreeEdge(u, v);
            }
        }
        dag_.BuildTree();
    }




    void CandidateSpace::ConstructTreeDP() {
//        fprintf(stderr, "DAGROOT = %d\n",dag_.GetRoot());
        double *num_tree_child = new double[data_.GetNumVertices()];
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByTree(query_.GetNumVertices() - i - 1);
            Label u_label = query_.GetLabel(u);
            num_trees_[u] = new double[GetCandidateSetSize(u)];
            std::fill(num_trees_[u], num_trees_[u] + GetCandidateSetSize(u), 1.0);
            for (Size child_idx = 0; child_idx < dag_.GetNumTreeChildren(u); ++child_idx) {
                std::fill(num_tree_child, num_tree_child + data_.GetNumVertices(), 0.0);
                Vertex child = dag_.GetTreeChild(u, child_idx);
//                fprintf(stderr, "Consider child %d of %d\n",child,u);
                for (Size child_cand_idx = 0; child_cand_idx < candidate_set_size_[child]; ++child_cand_idx) {
                    Vertex child_cand = candidate_set_[child][child_cand_idx];
                    for (Size cc_nbr_idx = data_.GetStartOffset(child_cand, u_label);
                         cc_nbr_idx < data_.GetEndOffset(child_cand, u_label); ++cc_nbr_idx) {
                        Vertex cand = data_.GetNeighbor(cc_nbr_idx);
                        num_tree_child[cand] += num_trees_[child][child_cand_idx];
                    }
                }
                for (Size c_idx = 0; c_idx < GetCandidateSetSize(u); ++c_idx) {
                    num_trees_[u][c_idx] *= num_tree_child[candidate_set_[u][c_idx]];
//                    fprintf(stderr,"NUM_TREE[%d][%d] = %.02lf\n", u, c_idx, num_trees_[u][c_idx]);
                }
            }
        }
        sample_dist.resize(query_.GetNumVertices());
        candidate_weights_ = new double***[query_.GetNumVertices()];
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByTree(query_.GetNumVertices() - i - 1);
            candidate_weights_[u] = new double**[dag_.GetNumTreeChildren(u)];
            sample_dist[u].resize(dag_.GetNumTreeChildren(u));
            for (Size ic = 0; ic < dag_.GetNumTreeChildren(u); ++ic) {
                Vertex uc = dag_.GetTreeChild(u, ic);
                candidate_weights_[u][ic] = new double*[GetCandidateSetSize(u)];
                sample_dist[u][ic].resize(GetCandidateSetSize(u));
                for (Size iv = 0; iv < GetCandidateSetSize(u); ++iv) {
                    Vertex v = GetCandidate(u, iv);
                    candidate_weights_[u][ic][iv] = new double[GetCandidateSetSize(uc)];
                    for (Size ivc = 0; ivc < GetCandidateSetSize(uc); ++ivc) {
                        Vertex vc = GetCandidate(uc, ivc);
                        candidate_weights_[u][ic][iv][ivc] = num_trees_[uc][ivc] * (data_.CheckEdgeExist(v, vc) ? 1 : 0) ;
                        //if (ivc > 0) candidate_weights_[u][ic][iv][ivc] += candidate_weights_[u][ic][iv][ivc-1];
                    }
                    sample_dist[u][ic][iv] = std::discrete_distribution<int>(candidate_weights_[u][ic][iv], candidate_weights_[u][ic][iv] + GetCandidateSetSize(uc));
                }
            }
        }
        Vertex root = dag_.GetRoot();
        std::vector <double> root_weight(GetCandidateSetSize(root), 0.0);
        for (int root_candidate_idx = 0; root_candidate_idx < GetCandidateSetSize(root); ++root_candidate_idx) {
            root_weight[root_candidate_idx] = num_trees_[root][root_candidate_idx];
        }
        root_weights_ = std::discrete_distribution<int>(root_weight.begin(), root_weight.end());
        total_trees_ = 0.0;
        for (Size i = 0; i < GetCandidateSetSize(root); ++i) {
            total_trees_ += num_trees_[root][i];
        }
        delete[] num_tree_child;
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    inline Size sample(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }

    void CandidateSpace::SampleCSTree(Vertex *tree_sample) {
        memset(tree_sample, -1, sizeof(Vertex) * query_.GetNumVertices());
        tree_sample[dag_.GetRoot()] = sample(root_weights_);
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByTree(i);
            for (Size ic = 0; ic < dag_.GetNumTreeChildren(u); ++ic) {
                Vertex uc = dag_.GetTreeChild(u, ic);
                tree_sample[uc] = sample(sample_dist[u][ic][tree_sample[u]]);
            }
        }
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            tree_sample[i] = GetCandidate(i, tree_sample[i]);
        }
    }

    double CandidateSpace::EstimateEmbeddings(Size num_samples) {
        bool *seen = new bool[data_.GetNumVertices()];
        memset(seen, 0, sizeof(bool) * data_.GetNumVertices());
        Vertex *tree_sample = new Vertex[query_.GetNumVertices()];
        Size success = 0, t;
        for (t = 1; t <= num_samples; ++t) {
            SampleCSTree(tree_sample);
            for (Size i = 0; i < query_.GetNumVertices(); ++i) {
                Vertex cand = tree_sample[i];
                if (seen[cand]) {
                    goto next_sample;
                }
                seen[cand] = true;
            }
            for (auto &[i, j] : query_.edge_exists) {
                if (!data_.CheckEdgeExist(tree_sample[i], tree_sample[j])) {
                    goto next_sample;
                }

            }
            success++;
            next_sample:
            for (Size i = 0; i < query_.GetNumVertices(); ++i) {
                seen[tree_sample[i]] = false;
            }
            double rhohat = (success * 1.0 / t);
            if (t >= 1000.0 / rhohat) break;
        }
//        fprintf(stderr,"TOTAL TREES %.04lf, SAMPLE %u/%u -> Estimate %.04lf\n",
//                total_trees_, success, num_samples, total_trees_ * (success * 1.0 / t));
        delete[] seen;
        delete[] tree_sample;
        return total_trees_ * (success * 1.0 / t);
    }
}
