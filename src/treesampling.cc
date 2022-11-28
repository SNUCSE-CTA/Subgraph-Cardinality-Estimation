#include "include/treesampling.h"
namespace daf {
    std::random_device rd;
    std::mt19937 gen(rd());
    inline Size sample(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }
    inline Size sample_no_weight(std::vector<Vertex> &candidate_indices) {
        return candidate_indices[rand() % candidate_indices.size()];
    }

    TreeSampling::TreeSampling(DataGraph &data, QueryGraph &query,
                               DAG &dag)
            : data_(data), query_(query), dag_(dag), CS(data, query, dag) {
        CS.BuildCS();
        BuildQueryTree();
        ConstructTreeDP();
        seen_.resize(data_.GetNumVertices(), false);
    }

    TreeSampling::~TreeSampling() {}

    void TreeSampling::BuildQueryTree() {
//        for (auto &e : CS.cs_edge_list_) {
//            VertexPair u_pair = e.first;
//            fprintf(stderr, "[%u, %u] -> ",u_pair.first,u_pair.second);
//            for (auto &ee : e.second) {
//                fprintf(stderr, "[%u, %u] ",ee.first,ee.second);
//            }
//            fprintf(stderr, "\n");
//        }
        std::vector<std::pair<double, std::pair<Vertex, Vertex>>> edges;
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            for (Size nidx = query_.GetStartOffset(i); nidx < query_.GetEndOffset(i); ++nidx) {
                Vertex q_neighbor = query_.GetNeighbor(nidx);
                if (i > q_neighbor) continue;
                double ij_cs_edge = 0;
                for (Size cs_idx = 0; cs_idx < CS.GetCandidateSetSize(i); ++cs_idx) {
                    Size num_cs_neighbor = 0;
                    for (VertexPair vp : CS.cs_edge_list_[{i, cs_idx}]) {
                        if (vp.first == q_neighbor)
                            num_cs_neighbor++;
                    }
                    ij_cs_edge += num_cs_neighbor;
                }
                ij_cs_edge /= ((1 + CS.GetCandidateSetSize(i)) * (1 + CS.GetCandidateSetSize(q_neighbor)));
                if (ij_cs_edge > 0) {
                    edges.push_back({ij_cs_edge * 1.0,{i, q_neighbor}});
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        UnionFind uf(query_.GetNumVertices());
        for (auto e : edges) {
            auto [u, v] = e.second;
            if (uf.unite(u, v)) {
//                fprintf(stderr, "%u %u\n",u,v);
                dag_.AddTreeEdge(u, v);
            }
        }
        dag_.BuildTree();
    }

    void TreeSampling::ConstructTreeDP() {
//        fprintf(stderr, "ConstructTreeDP\n");
        num_trees_.resize(query_.GetNumVertices());
        sample_candidates_.resize(query_.GetNumVertices());
        sample_candidate_weights_.resize(query_.GetNumVertices());
        sample_dist_.resize(query_.GetNumVertices());
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByTree(query_.GetNumVertices() - i - 1);
//            fprintf(stderr, "u = %u, parent = %u\n",u, dag_.GetTreeParent(u, 0));
            Size num_cands = CS.GetCandidateSetSize(u);
            Size num_children = dag_.GetNumTreeChildren(u);
            sample_candidates_[u].resize(num_cands);
            sample_candidate_weights_[u].resize(num_cands);
            sample_dist_[u].resize(num_cands);
            std::vector<double> tmp_num_child(num_children, 0.0);
            for (Size cs_idx = 0; cs_idx < num_cands; ++cs_idx) {
                sample_candidates_[u][cs_idx].resize(num_children);
                sample_candidate_weights_[u][cs_idx].resize(num_children);
                sample_dist_[u][cs_idx].resize(num_children);
                double num_ = 1.0;
                Vertex v = CS.GetCandidate(u, cs_idx);
                VertexPair u_pair = {u, cs_idx};
                std::fill(tmp_num_child.begin(), tmp_num_child.end(), 0.0);
                for (auto &[uc, vc_idx] : CS.cs_edge_list_[u_pair]){
                    if (dag_.GetTreeParent(uc, 0) != u) continue;
                    Vertex vc = CS.GetCandidate(uc, vc_idx);
                    Size uc_idx = dag_.GetChildIndex(uc);
//                    fprintf(stderr, "Consider [%u %u(%u)] -> [%u %u(%u)]\n",u,cs_idx,v, uc, vc_idx, vc);
                    tmp_num_child[uc_idx] += num_trees_[uc][vc];
                    sample_candidates_[u][cs_idx][uc_idx].push_back(vc_idx);
                    sample_candidate_weights_[u][cs_idx][uc_idx].push_back(num_trees_[uc][vc]);
                }
                for (Size j = 0; j < num_children; ++j) {
                    num_ *= tmp_num_child[j];
                    sample_dist_[u][cs_idx][j] = std::discrete_distribution<int>(
                            sample_candidate_weights_[u][cs_idx][j].begin(),
                            sample_candidate_weights_[u][cs_idx][j].end());
//                    sample_dist_[u][uc] = std::discrete_distribution<int>(
//                                        sample_candidate_weights_[u_pair][uc].begin(),
//                                        sample_candidate_weights_[u_pair][uc].end());
                }
                num_trees_[u][v] = num_;
            }
//            fprintf(stderr, "Numtree %u : ",u);
//            for (Size cand_idx = 0; cand_idx < CS.GetCandidateSetSize(u); ++cand_idx) {
//                Vertex v = CS.GetCandidate(u, cand_idx);
//                fprintf(stderr, "%u[%.04lf]\t",v,num_trees_[u][v]);
//            }
//            fprintf(stderr, "\n");
        }

        total_trees_ = 0.0;
        Vertex root = dag_.GetRoot();
        Size root_candidate_size = CS.GetCandidateSetSize(root);
        std::vector <double> root_weight;
        for (int root_candidate_idx = 0; root_candidate_idx < root_candidate_size; ++root_candidate_idx) {
            Vertex root_candidate = CS.GetCandidate(root, root_candidate_idx);
            total_trees_ += num_trees_[root][root_candidate];
            if (num_trees_[root][root_candidate] > 0) {
                root_candidates_.push_back(root_candidate_idx);
                root_weight.push_back(num_trees_[root][root_candidate]);
            }
        }
//        fprintf(stderr, "ROOTDISTR-ROOT %u : ",root);
        for (Size cand_idx = 0; cand_idx < CS.GetCandidateSetSize(root); ++cand_idx) {
            Vertex v = CS.GetCandidate(root, cand_idx);
//            fprintf(stderr, "%u[%.04lf]\t",v,num_trees_[root][v]);
        }
//        fprintf(stderr, "\n");
        sample_root_dist_ = std::discrete_distribution<int>(root_weight.begin(), root_weight.end());
    }

    double TreeSampling::SampleCSTree(std::vector<int> &tree_sample) {
        tree_sample[dag_.GetRoot()] = root_candidates_[sample(sample_root_dist_)];

//        fprintf(stderr, "DAGRoot = %u, sampled %uth\n", dag_.GetRoot(), tree_sample[dag_.GetRoot()]);
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByTree(i);
            Vertex v_idx = tree_sample[u];
            for (Size uc_idx = 0; uc_idx < dag_.GetNumTreeChildren(u); ++uc_idx) {
                Vertex uc = dag_.GetTreeChild(u, uc_idx);
                Size vc_idx = sample(sample_dist_[u][v_idx][uc_idx]);
//                fprintf(stderr, "vertex %u-%u, child %u => Sample %u...",u,v,uc,vc);
                tree_sample[uc] = sample_candidates_[u][v_idx][uc_idx][vc_idx];
            }
        }
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            tree_sample[i] = CS.GetCandidate(i, tree_sample[i]);
        }
//        fprintf(stderr, "SAMPLE : ");
//        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
//            fprintf(stderr, " %u(%u)", tree_sample[i], i);
//        }
//        fprintf(stderr, "\n");
        return 1.0;
    }

    int TreeSampling::CheckSample(std::vector<int> &sample) {
        int result = 1;
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            if (sample[i] >= 0) {
                Vertex cand = sample[i];
                if (seen_[cand]) {
                    result = 0;
                    goto done;
                }
                seen_[cand] = true;
            }
            else result = -1;
        }

        for (auto &[i, j] : query_.edge_exists) {
            if (sample[i] < 0 || sample[j] < 0) continue;
            if (!data_.CheckEdgeExist(sample[i], sample[j])) {
                result = 0;
                goto done;
            }
        }
        done:
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            if (sample[i] >= 0)
                seen_[sample[i]] = false;
        }
        return result;
    }

    double TreeSampling::EstimateEmbeddings(Size num_samples) {
        std::vector<int> tree_sample(query_.GetNumVertices(), -1);
        Size success = 0, t = 0;
        double ht_est = 0.0;
        for (t = 1; t <= num_samples; ++t) {
            std::fill(tree_sample.begin(), tree_sample.end(), -1);
            double ht_prob = SampleCSTree(tree_sample);
            if (CheckSample(tree_sample) == 1) {
                ht_est += ht_prob;
                success++;
            }
            double rhohat = (success * 1.0 / t);
            if (t >= 50000.0 / rhohat) break;
        }
        return total_trees_ * (success * 1.0 / t);
    }

}