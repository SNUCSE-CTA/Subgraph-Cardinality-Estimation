#include "include/treesampling.h"
#include "global/timer.h"

namespace daf {
    inline Size sample(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }

    TreeSampling::TreeSampling(DataGraph &data, QueryGraph &query,
                               DAG &dag)
            : data_(data), query_(query), dag_(dag), CS(data, query, dag) {
        Timer treeSamplingTimer; treeSamplingTimer.Start();
        CS.BuildCS();
        std::cout << "CSTIME " << treeSamplingTimer.Peek() << " ms\n";
        Timer querytree_timer; querytree_timer.Start();
        BuildQueryTree();
        querytree_timer.Stop();
        std::cout << "Query Tree Building Time: " << treeSamplingTimer.Peek() << " ms\n";
        RWI_ = RWI();
        RWI_.init(data, query, dag, CS);
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
                    int v = CS.GetCandidate(i, cs_idx);
                    Size num_cs_neighbor = CS.cs_adj_[{i, v}][q_neighbor].size();
//                    Size num_cs_neighbor = 0;
//                    for (VertexPair vp : CS.cs_edge_list_[{i, cs_idx}]) {
//                        if (vp.first == q_neighbor)
//                            num_cs_neighbor++;
//                    }
                    ij_cs_edge += num_cs_neighbor;
                }
                ij_cs_edge /= ((1 + CS.GetCandidateSetSize(i)) * (1 + CS.GetCandidateSetSize(q_neighbor)));
                if (ij_cs_edge > 0) {
                    edges.push_back({ij_cs_edge * 1.0,{i, q_neighbor}});
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        std::vector<std::pair<double, std::pair<Vertex, Vertex>>> cur_best_config;
        double cur_best_trees = std::numeric_limits<double>::max();
        for (int build_tree_iteration = 0; build_tree_iteration < 2; build_tree_iteration++) {
            UnionFind uf(query_.GetNumVertices());
            for (auto &t_neighbor : dag_.tree_neighbors_) {
                t_neighbor.clear();
            }
            for (auto e : edges) {
                auto [u, v] = e.second;
                if (uf.unite(u, v)) {
//                fprintf(stderr, "%u %u\n",u,v);
                    dag_.AddTreeEdge(u, v);
                }
            }
            dag_.BuildTree();
            double num_tree_homo = ConstructTreeDP();
            if (num_tree_homo < cur_best_trees) {
                cur_best_trees = num_tree_homo;
                cur_best_config = edges;
            }
//            fprintf(stderr, "ITER %d NUM_TREE = %.4le BEST %.04le\n",build_tree_iteration, num_tree_homo, cur_best_trees);
            if (build_tree_iteration == 1) {
                for (auto &e : edges) {
                    int i, q_neighbor;
                    std::tie(i, q_neighbor) = e.second;
                    e.first *= ((1 + CS.GetCandidateSetSize(i)) * (1 + CS.GetCandidateSetSize(q_neighbor)));
                }
                std::sort(edges.begin(), edges.end());
            }
            else{
                std::shuffle(edges.begin(), edges.end(), gen);
            }
        }
        edges = cur_best_config;
        for (auto &t_neighbor : dag_.tree_neighbors_) {
            t_neighbor.clear();
        }
        UnionFind uf(query_.GetNumVertices());
        for (auto e : edges) {
            auto [u, v] = e.second;
            if (uf.unite(u, v)) {
//                fprintf(stderr, "%u %u\n",u,v);
                dag_.AddTreeEdge(u, v);
            }
        }
        dag_.BuildTree();
        double num_tree_homo = ConstructTreeDP();
        fprintf(stderr, "NUM_TREE = %.04lE\n",num_tree_homo);
    }

    double TreeSampling::ConstructTreeDP() {
//        fprintf(stderr, "ConstructTreeDP\n");
        num_trees_.clear();
        num_trees_.resize(query_.GetNumVertices());
        for (auto &it : num_trees_) it.clear();
        sample_candidates_.clear();
        sample_candidates_.resize(query_.GetNumVertices());
        sample_candidate_weights_.clear();
        sample_candidate_weights_.resize(query_.GetNumVertices());
        sample_dist_.clear();
        sample_dist_.resize(query_.GetNumVertices());
        for (Size i = 0; i < query_.GetNumVertices(); ++i) {
            Vertex u = dag_.GetVertexOrderedByTree(query_.GetNumVertices() - i - 1);
//            fprintf(stderr, "u = %u, parent = %u\n",u, dag_.GetTreeParent(u, 0));
            Size num_cands = CS.GetCandidateSetSize(u);
            Size num_children = dag_.GetNumTreeChildren(u);

            sample_candidates_[u].clear();
            sample_candidates_[u].resize(num_cands);
            sample_candidate_weights_[u].clear();
            sample_candidate_weights_[u].resize(num_cands);
            sample_dist_[u].clear();
            sample_dist_[u].resize(num_cands);

            std::vector<double> tmp_num_child(num_children);
            for (Size cs_idx = 0; cs_idx < num_cands; ++cs_idx) {

                sample_candidates_[u][cs_idx].clear();
                sample_candidates_[u][cs_idx].resize(num_children);
                sample_candidate_weights_[u][cs_idx].clear();
                sample_candidate_weights_[u][cs_idx].resize(num_children);
                sample_dist_[u][cs_idx].clear();
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
                }
                num_trees_[u][v] = num_;
            }
        }

        total_trees_ = 0.0;
        root_candidates_.clear();
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
//        fprintf(stderr, "\n");
        sample_root_dist_ = std::discrete_distribution<int>(root_weight.begin(), root_weight.end());
        return total_trees_;
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
//                    fprintf(stderr, "Dead by Noninjectiveness : %d is seen twice\n",cand);
                    result = -1;
                    goto done;
                }
                seen_[cand] = true;
            }
        }
        if (result == -1) goto done;

        for (auto &[i, j] : query_.all_edges) {
            if (sample[i] < 0 || sample[j] < 0) continue;
            if (!data_.CheckEdgeExist(sample[i], sample[j])) {
                result = 0;
//                fprintf(stderr, "There should be an edge between %u and %u but there isnt in %u %u\n",
//                        i,j,sample[i],sample[j]);
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

    std::pair<double, int> TreeSampling::UniformSamplingEstimate(Size num_samples) {
        std::vector<int> tree_sample(query_.GetNumVertices(), -1);
        Size success = 0, t = 0;
        double ht_est = 0.0;
        int reject_homo = 0, reject_nontree = 0;
        while (num_samples--) {
            std::fill(tree_sample.begin(), tree_sample.end(), -1);
            double ht_prob = SampleCSTree(tree_sample);
            t++;
            int sample_record = CheckSample(tree_sample);
            if (sample_record == 1) {
                ht_est += ht_prob;
                success++;
            }
            else if (sample_record == -1) {
                reject_homo++;
            }
            else if (sample_record == 0) {
                reject_nontree++;
            }
            double rhohat = (success * 1.0 / t);
            if (success >= 1000) break;
        }
        fprintf(stdout, "#NUM_SAMPLES : %u\n", t);
        fprintf(stdout, "#NUM_SUCCESS : %u\n", success);
        fprintf(stderr,"Rejected by Injective : %d, Rejected by Edge : %d\n", reject_homo, reject_nontree);
        return {total_trees_ * (success * 1.0 / t), success};
    }


    void inplace_set_intersection(std::set<int> &st1, std::set<int> &st2) {
        std::set<int>::iterator it1 = st1.begin();
        std::set<int>::iterator it2 = st2.begin();
        while ( (it1 != st1.end()) && (it2 != st2.end()) ) {
            if (*it1 < *it2) {
                st1.erase(it1++);
            }
            else if (*it2 < *it1) {
                ++it2;
            }
            else {
                ++it1;
                ++it2;
            }
        }
        st1.erase(it1, st1.end());
    }

    double TreeSampling::EstimateEmbeddings(Size num_samples) {
//        Timer sampletimer_uni, sampletimer_inter;
//        sampletimer_uni.Start();
//        std::pair<double, int> uniformResult = UniformSamplingEstimate(num_samples);
//        sampletimer_uni.Stop();
//        std::cout << "Uniform Sampling time: " << std::fixed << sampletimer_uni.GetTime() << " ms\n";
//        if (uniformResult.second == 0) {
//            sampletimer_inter.Start();
//            double intersectionResult = IntersectionSamplingEstimate(num_samples);
//            sampletimer_inter.Stop();
//            std::cout << "Intersection Sampling time: " << std::fixed << sampletimer_inter.GetTime() << " ms\n";
//            return intersectionResult;
//        }
//        return uniformResult.first;
//        CS.printCS();
        double intersectionResult = RWI_.IntersectionSamplingEstimate(1000000);
        return intersectionResult;
    }

}