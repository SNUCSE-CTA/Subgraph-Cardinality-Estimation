#include "include/sampling.h"
namespace daf {
    TreeSampling::TreeSampling(const DataGraph &data, const QueryGraph &query,
                               DAG &dag)
            : data_(data), query_(query), dag_(dag), CS(data, query, dag) {
        CS.BuildCS();
        BuildQueryTree();
//        fprintf(stderr, "QueryTree:\t");
//        for (Vertex i = 0; i < query_.GetNumVertices(); i++) {
//            int num_children = dag_.GetNumTreeChildren(i);
//            fprintf(stderr, "Tree_adj[%u] : ", i);
//            for (Vertex j = 0; j < num_children; j++) {
//                fprintf(stderr, "%u[%u]\t", dag_.GetTreeChild(i, j), dag_.GetChildIndex(dag_.GetTreeChild(i, j)));
//            }
//            fprintf(stderr, "\n");
//        }
        ConstructTreeDP();
    }
    TreeSampling::~TreeSampling() {

    }
    void TreeSampling::BuildQueryTree() {

        std::vector<std::pair<double, std::pair<Vertex, Vertex>>> edges;
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            for (Size nidx = query_.GetStartOffset(i); nidx < query_.GetEndOffset(i); ++nidx) {
                Vertex q_neighbor = query_.GetNeighbor(nidx);
                if (i > q_neighbor) continue;
                Size ij_cs_edge = 0;
                for (Size cs_idx = 0; cs_idx < CS.GetCandidateSetSize(i); ++cs_idx) {
                    Size num_cs_neighbor =
                            CS.GetCandidateEndOffset(i, nidx-query_.GetStartOffset(i), cs_idx) -
                            CS.GetCandidateStartOffset(i, nidx-query_.GetStartOffset(i), cs_idx);
                    ij_cs_edge += num_cs_neighbor;
                }
                ij_cs_edge /= ((1 + CS.GetCandidateSetSize(i)) * (1 + CS.GetCandidateSetSize(q_neighbor)));
                edges.push_back({ij_cs_edge * 1.0,{i, q_neighbor}});
            }
        }
        std::sort(edges.begin(), edges.end());
        UnionFind uf(query_.GetNumVertices());
        for (auto e : edges) {
            auto [u, v] = e.second;
            if (uf.unite(u, v)) {
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
        sample_root_dist_ = std::discrete_distribution<int>(root_weight.begin(), root_weight.end());
    }

    std::random_device rd;
    std::mt19937 gen(rd());
    inline Size sample(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }
    inline Size sample_no_weight(std::vector<Vertex> &candidate_indices) {
        return candidate_indices[rand() % candidate_indices.size()];
    }

    double TreeSampling::SampleCSTree(Vertex *tree_sample, bool HTSampling) {
        memset(tree_sample, -1, sizeof(Vertex) * query_.GetNumVertices());
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

    double TreeSampling::EstimateEmbeddings(Size num_samples, bool HTSampling) {
        bool *seen = new bool[data_.GetNumVertices()];
        memset(seen, 0, sizeof(bool) * data_.GetNumVertices());
        Vertex *tree_sample = new Vertex[query_.GetNumVertices()];
        Size success = 0, t;
        double ht_est = 0.0;
        for (t = 1; t <= num_samples; ++t) {
            double ht_prob = SampleCSTree(tree_sample, HTSampling);
            for (Size i = 0; i < query_.GetNumVertices(); ++i) {
                Vertex cand = tree_sample[i];
                if (seen[cand]) {
                    goto next_sample;
                }
                seen[cand] = true;
            }
            for (auto &[i, j] : query_.edge_exists) {
                if (!data_.CheckEdgeExist(tree_sample[i], tree_sample[j])) {
//                    fprintf(stderr, "Sample : ");
//                    for (int ii = 0; ii < query_.GetNumVertices(); ++ii) {
//                        fprintf(stderr, "%u(%u) ", ii,tree_sample[ii]);
//                    }
//                    fprintf(stderr, "\n");
//                    exit(1);
                    goto next_sample;
                }
            }
            ht_est += ht_prob;
            success++;
            next_sample:
            for (Size i = 0; i < query_.GetNumVertices(); ++i) {
                seen[tree_sample[i]] = false;
            }
            if (!HTSampling) {
                double rhohat = (success * 1.0 / t);
                if (t >= 1000.0 / rhohat) break;
            }
        }
        delete[] seen;
        delete[] tree_sample;
        if (HTSampling) return ht_est / t;
        else return total_trees_ * (success * 1.0 / t);
    }

}