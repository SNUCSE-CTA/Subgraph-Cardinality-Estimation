#include "include/treesampling.h"
#include "global/timer.h"

namespace daf {
    inline Size sample(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }

    TreeSampling::TreeSampling(DataGraph *data) {
        data_ = data;
        CS = new CandidateSpace(data);
        seen_.resize(data_->GetNumVertices());
        RWI_ = new RWI(data);
        
    }

    void TreeSampling::RegisterQuery(QueryGraph *query, DAG *dag) {
        query_ = query;
        dag_ = dag;
        Timer CSTimer; CSTimer.Start();
        CS->BuildCS(query, dag);
        std::cout << "CSTIME " << CSTimer.Peek() << " ms\n";
        std::fill(seen_.begin(), seen_.end(), 0);
        sample_dist_.clear();
        num_trees_.clear();
        sample_candidates_.clear();
        sample_candidate_weights_.clear();
        sample_dist_.clear();
        root_candidates_.clear();
        Timer querytree_timer; querytree_timer.Start();
        BuildQueryTree();
        querytree_timer.Stop();
        std::cout << "Query Tree Building Time: " << querytree_timer.GetTime() << " ms\n";
        RWI_->init(query, CS);
    }

    TreeSampling::~TreeSampling() {}

    void TreeSampling::BuildQueryTree() {
        std::vector<std::pair<double, std::pair<Vertex, Vertex>>> edges;
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            for (Size nidx = query_->GetStartOffset(i); nidx < query_->GetEndOffset(i); ++nidx) {
                Vertex q_neighbor = query_->GetNeighbor(nidx);
                if (i > q_neighbor) continue;
                double ij_cs_edge = 0;
                for (Size cs_idx = 0; cs_idx < CS->GetCandidateSetSize(i); ++cs_idx) {
                    Size num_cs_neighbor = CS->cs_edge_[i][cs_idx][q_neighbor].size();
                    //ij_cs_edge += num_cs_neighbor * num_cs_neighbor;
                    ij_cs_edge += num_cs_neighbor;
                }
                
                ij_cs_edge /= (CS->GetCandidateSetSize(i) * CS->GetCandidateSetSize(q_neighbor));
                if (ij_cs_edge > 0) {
                    edges.push_back({ij_cs_edge * 1.0,{i, q_neighbor}});
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        UnionFind uf(query_->GetNumVertices());
        for (auto e : edges) {
            auto [u, v] = e.second;
            if (uf.unite(u, v)) {
                dag_->AddTreeEdge(u, v);
            }
        }
        dag_->BuildTree();
        double num_tree_homo = ConstructTreeDP();
        fprintf(stderr, "NUM_TREE = %.04lE\n",num_tree_homo);
    }

    double TreeSampling::ConstructTreeDP() {
        num_trees_.resize(query_->GetNumVertices());
        for (auto &it : num_trees_) it.clear();
        sample_candidates_.resize(query_->GetNumVertices());
        sample_candidate_weights_.resize(query_->GetNumVertices());
        sample_dist_.resize(query_->GetNumVertices());
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            Vertex u = dag_->GetVertexOrderedByTree(query_->GetNumVertices() - i - 1);
//            fprintf(stderr, "u = %u, parent = %u\n",u, dag_->GetTreeParent(u, 0));
            Size num_cands = CS->GetCandidateSetSize(u);
            Size num_children = dag_->GetNumTreeChildren(u);

            sample_candidates_[u].resize(num_cands);
            sample_candidate_weights_[u].resize(num_cands);
            sample_dist_[u].resize(num_cands);

            std::vector<double> tmp_num_child(num_children);
            for (Size cs_idx = 0; cs_idx < num_cands; ++cs_idx) {
                sample_candidates_[u][cs_idx].resize(num_children);
                sample_candidate_weights_[u][cs_idx].resize(num_children);
                sample_dist_[u][cs_idx].resize(num_children);

                double num_ = 1.0;
                Vertex v = CS->GetCandidate(u, cs_idx);
                std::fill(tmp_num_child.begin(), tmp_num_child.end(), 0.0);
                for (auto &uc : query_->adj_list[u]){
                    if (dag_->GetTreeParent(uc, 0) != u) continue;
                    for (auto &vc_idx : CS->cs_edge_[u][cs_idx][uc]) {
                        Vertex vc = CS->GetCandidate(uc, vc_idx);
                        Size uc_idx = dag_->GetChildIndex(uc);
//                    fprintf(stderr, "Consider [%u %u(%u)] -> [%u %u(%u)]\n",u,cs_idx,v, uc, vc_idx, vc);
                        tmp_num_child[uc_idx] += num_trees_[uc][vc];
                        sample_candidates_[u][cs_idx][uc_idx].emplace_back(vc_idx);
                        sample_candidate_weights_[u][cs_idx][uc_idx].emplace_back(num_trees_[uc][vc]);
                    }
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
        Vertex root = dag_->GetRoot();
        Size root_candidate_size = CS->GetCandidateSetSize(root);
        std::vector <double> root_weight;
        for (int root_candidate_idx = 0; root_candidate_idx < root_candidate_size; ++root_candidate_idx) {
            Vertex root_candidate = CS->GetCandidate(root, root_candidate_idx);
            total_trees_ += num_trees_[root][root_candidate];
            if (num_trees_[root][root_candidate] > 0) {
                root_candidates_.emplace_back(root_candidate_idx);
                root_weight.emplace_back(num_trees_[root][root_candidate]);
            }
        }
//        fprintf(stderr, "\n");
        sample_root_dist_ = std::discrete_distribution<int>(root_weight.begin(), root_weight.end());
        return total_trees_;
    }

    double TreeSampling::SampleCSTree(std::vector<int> &tree_sample) {
        tree_sample[dag_->GetRoot()] = root_candidates_[sample(sample_root_dist_)];
        seen_[CS->GetCandidate(dag_->GetRoot(), tree_sample[dag_->GetRoot()])] = true;
        bool valid = true;
//        fprintf(stderr, "DAGRoot = %u, sampled %uth\n", dag_->GetRoot(), tree_sample[dag_->GetRoot()]);
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            Vertex u = dag_->GetVertexOrderedByTree(i);
            Vertex v_idx = tree_sample[u];
            for (Size uc_idx = 0; uc_idx < dag_->GetNumTreeChildren(u); ++uc_idx) {
                Vertex uc = dag_->GetTreeChild(u, uc_idx);
                Size vc_idx = sample(sample_dist_[u][v_idx][uc_idx]);
//                fprintf(stderr, "vertex %u-%u, child %u => Sample %u...",u,v,uc,vc);
                tree_sample[uc] = sample_candidates_[u][v_idx][uc_idx][vc_idx];
                int cand = CS->GetCandidate(uc, tree_sample[uc]);
                if (seen_[cand]) {
                    valid = false;
                    goto INJECTIVE_VIOLATED;
                }
                seen_[cand] = true;
            }
        }
        INJECTIVE_VIOLATED:
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            if (tree_sample[i] >= 0) {
                int cand = CS->GetCandidate(i, tree_sample[i]);
                seen_[cand] = false;
            }
        }
        if (!valid) return 0.0;
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            tree_sample[i] = CS->GetCandidate(i, tree_sample[i]);
        }
        return 1.0;
    }

    int TreeSampling::CheckSample(std::vector<int> &sample) {
        for (auto &[i, j] : query_->all_edges) {
            if (sample[i] < 0 || sample[j] < 0) continue;
            if (!data_->CheckEdgeExist(sample[i], sample[j])) {
                return false;
            }
        }
        return true;
    }

    boost::math::normal dist(0., 1.);
    long double z = quantile(dist, 0.975);
    std::pair<double, int> TreeSampling::UniformSamplingEstimate() {
        std::vector<int> tree_sample(query_->GetNumVertices(), -1);
        long long success = 0, t = 0;
        int reject_homo = 0, reject_nontree = 0;
        while (true) {
            std::fill(tree_sample.begin(), tree_sample.end(), -1);
            double sample_result = SampleCSTree(tree_sample);
            t++;
            if (sample_result > 0.5) {
                int sample_record = CheckSample(tree_sample);
                if (sample_record == 1) {
                    success++;
                }
                else if (sample_record == 0) {
                    reject_nontree++;
                }
            }
            else reject_homo++;
            long double rhohat = (success * 1.0 / t);
            if(t>=1000 && t%100==0){
                if(t==50000 && success*5000<=t){
                    fprintf(stderr, "NUM_HOMO  %d\tNUM_CFLCT %d\n",reject_homo, reject_nontree);
                    fprintf(stderr, "#NUM_SAMPLES : %lld, #NUM_SUCCESS : %lld\n", t, success);
                    fprintf(stdout, "#NUM_SAMPLES : %lld\n", t);
                    fprintf(stdout, "#NUM_SUCCESS : %lld\n", success);
                    return {-1, success};
                }
//                if(success >= 100) break;
                // Wilson Confidence Interval
                long double wminus = (success + z*z/2 - z*sqrt(success*(t-success)/(long double)t + z*z/4))/(t+z*z);
                long double wplus  = (success + z*z/2 + z*sqrt(success*(t-success)/(long double)t + z*z/4))/(t+z*z);

                if(rhohat * 0.8 < wminus && wplus < rhohat * 1.25){
                    break;
                }
/*
                //Agresti-Coull
                boost::math::normal dist(0., 1.);
                long double z = quantile(dist, 0.975);
                long double ntilde = t + z*z;
                long double rhotilde = (success+z*z/2)/ntilde;
                long double wminus = rhotilde - z*sqrt(rhotilde*(1-rhotilde)/ntilde);
                long double wplus = rhotilde + z*sqrt(rhotilde*(1-rhotilde)/ntilde);
                //if(rhohat * 0.8 < wminus && wplus < rhohat * 1.25) break;
                if(rhohat * (1-t/1500000.0) < wminus && wplus < rhohat * (1+t/1500000.0)) break;
*/
            }
        }
        fprintf(stderr, "NUM_HOMO  %d\tNUM_CFLCT %d\n",reject_homo, reject_nontree);
        fprintf(stderr, "#NUM_SAMPLES : %lld, #NUM_SUCCESS : %lld\n", t, success);

        fprintf(stdout, "#NUM_SAMPLES : %lld\n", t);
        fprintf(stdout, "#NUM_SUCCESS : %lld\n", success);
//        fprintf(stderr,"Rejected by Injective : %d, Rejected by Edge : %d\n", reject_homo, reject_nontree);
        return {total_trees_ * (success * 1.0 / t), success};
    }

    double TreeSampling::EstimateEmbeddings(Size num_samples) {
        Timer sampletimer_uni, sampletimer_inter;
        sampletimer_uni.Start();
        std::pair<double, int> uniformResult = UniformSamplingEstimate();
        sampletimer_uni.Stop();
        std::cout << "Uniform Sampling time: " << std::fixed << sampletimer_uni.GetTime() << " ms\n";
        if (uniformResult.first < 0) {
            sampletimer_inter.Start();
            double intersectionResult = RWI_->IntersectionSamplingEstimate(1000000 / (uniformResult.second + 1));
            sampletimer_inter.Stop();
            std::cout << "Intersection Sampling time: " << std::fixed << sampletimer_inter.GetTime() << " ms\n";
            return intersectionResult;
        }
        return uniformResult.first;
//        CS->printCS();
//        double intersectionResult = RWI_.IntersectionSamplingEstimate(1000000);
//        return intersectionResult;
    }

}