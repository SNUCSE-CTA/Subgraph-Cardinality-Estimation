#include <stack>
#include <queue>
#include "include/FaSTest.h"
#include "global/global.h"
#include "global/timer.h"

namespace CardEst {
    inline Size sample_from_distribution(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }

    FaSTest::FaSTest(DataGraph *data, Option opt_) {
        data_ = data;
        opt = opt_;
        CS = new CandidateSpace(data, opt);
        seen_.resize(data_->GetNumVertices());
        GSSolver = new GraphSampling(data, opt);
        num_trees_ = new double*[MAX_QUERY_VERTEX];
        for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
            num_trees_[i] = new double[data_->GetNumVertices()];
        }
    }

    void FaSTest::PrepareQuery(QueryGraph *query, OrderedQueryGraph *dag) {
        query_ = query;
        dag_ = dag;
        CS->BuildCS(query, dag);
        std::fill(seen_.begin(), seen_.end(), 0);
        sample_dist_.clear();
        sample_candidates_.clear();
        sample_candidate_weights_.clear();
        sample_dist_.clear();
        root_candidates_.clear();
        Timer querytree_timer; querytree_timer.Start();
        BuildQueryTree();
        querytree_timer.Stop();
        GSSolver->init(query, CS);
    }

    FaSTest::~FaSTest() {
        for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
            delete[] num_trees_[i];
        }
        delete[] num_trees_;
    }

    std::vector<std::pair<Vertex, Vertex>> BuildTreeByPrim(std::vector<std::pair<double, std::pair<Vertex, Vertex>>> &edges) {
        std::vector<std::pair<Vertex, Vertex>> T;
        std::vector<std::vector<std::pair<double, int>>> G(50);
        int qV = 0;
        for (auto &[w, e] : edges) {
            auto &[u, v] = e;
            G[u].emplace_back(w, v);
            G[v].emplace_back(w, u);
            qV = std::max(qV, std::max(u, v)+1);
        }
        for (int i = 0; i < qV; i++) {
            std::sort(G[i].begin(), G[i].end());
        }
        std::vector<int> deg(qV, 0);
        UnionFind uf(qV);
        while (T.size() + 1 < qV) {
            double minw = 1e9;
            std::pair<Vertex, Vertex> me;
            for (auto &[w, e] : edges) {
                auto [u, v] = e;
                if (uf.find(u) == uf.find(v)) continue;
                if (minw > w + deg[u] * 1e-7) {
                    minw = w + deg[u] * 1e-7;
                    me = e;
                }
            }
            uf.unite(me.first, me.second);
            T.push_back(me);
            deg[me.first]++;
            deg[me.second]++;
        }
        return T;
    }

    void FaSTest::BuildQueryTree() {
        std::vector<std::pair<double, std::pair<Vertex, Vertex>>> edges;
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            for (Size nidx = query_->GetStartOffset(i); nidx < query_->GetEndOffset(i); ++nidx) {
                Vertex q_neighbor = query_->GetNeighbor(nidx);
                if (i > q_neighbor) continue;
                double ij_cs_edge = 0;
                for (Size cs_idx = 0; cs_idx < CS->GetCandidateSetSize(i); ++cs_idx) {
                    Size num_cs_neighbor = CS->cs_edge_[i][cs_idx][q_neighbor].size();
                    ij_cs_edge += num_cs_neighbor;
                }
                
                ij_cs_edge /= (CS->GetCandidateSetSize(i) * CS->GetCandidateSetSize(q_neighbor));
                if (ij_cs_edge > 0) {
                    edges.push_back({ij_cs_edge * 1.0,{i, q_neighbor}});
                }
            }
        }
        std::sort(edges.begin(), edges.end());
        auto T = BuildTreeByPrim(edges);
        for (auto e : T) {
            auto [u, v] = e;
            dag_->AddTreeEdge(u, v);
        }
        dag_->BuildTree();
        double num_tree_homo = CountCandidateTrees();
    }

    double FaSTest::CountCandidateTrees() {
        sample_candidates_.resize(query_->GetNumVertices());
        sample_candidate_weights_.resize(query_->GetNumVertices());
        sample_dist_.resize(query_->GetNumVertices());
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            Vertex u = dag_->GetVertexOrderedByTree(query_->GetNumVertices() - i - 1);
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
                        tmp_num_child[uc_idx] += num_trees_[uc][vc_idx];
                        sample_candidates_[u][cs_idx][uc_idx].emplace_back(vc_idx);
                        sample_candidate_weights_[u][cs_idx][uc_idx].emplace_back(num_trees_[uc][vc_idx]);
                    }
                }
                for (Size j = 0; j < num_children; ++j) {
                    num_ *= tmp_num_child[j];
                    sample_dist_[u][cs_idx][j] = std::discrete_distribution<int>(
                            sample_candidate_weights_[u][cs_idx][j].begin(),
                            sample_candidate_weights_[u][cs_idx][j].end());
                }
                num_trees_[u][cs_idx] = num_;
            }
        }

        total_trees_ = 0.0;
        root_candidates_.clear();
        Vertex root = dag_->GetRoot();
        Size root_candidate_size = CS->GetCandidateSetSize(root);
        std::vector <double> root_weight;
        for (int root_candidate_idx = 0; root_candidate_idx < root_candidate_size; ++root_candidate_idx) {
            total_trees_ += num_trees_[root][root_candidate_idx];
            if (num_trees_[root][root_candidate_idx] > 0) {
                root_candidates_.emplace_back(root_candidate_idx);
                root_weight.emplace_back(num_trees_[root][root_candidate_idx]);
            }
        }
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            memset(num_trees_[i], 0, sizeof(double) * CS->GetCandidateSetSize(i));
        }
        sample_root_dist_ = std::discrete_distribution<int>(root_weight.begin(), root_weight.end());
        return total_trees_;
    }

    double FaSTest::GetSampleTree(std::vector<int> &sample) {
        sample[dag_->GetRoot()] = root_candidates_[sample_from_distribution(sample_root_dist_)];
        seen_[CS->GetCandidate(dag_->GetRoot(), sample[dag_->GetRoot()])] = true;
        bool valid = true;
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            Vertex u = dag_->GetVertexOrderedByTree(i);
            Vertex v_idx = sample[u];
            for (Size uc_idx = 0; uc_idx < dag_->GetNumTreeChildren(u); ++uc_idx) {
                Vertex uc = dag_->GetTreeChild(u, uc_idx);
                Size vc_idx = sample_from_distribution(sample_dist_[u][v_idx][uc_idx]);
                sample[uc] = sample_candidates_[u][v_idx][uc_idx][vc_idx];
                int cand = CS->GetCandidate(uc, sample[uc]);
                if (seen_[cand]) {
                    valid = false;
                    goto INJECTIVE_VIOLATED;
                }
                seen_[cand] = true;
            }
        }
        INJECTIVE_VIOLATED:
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            if (sample[i] >= 0) {
                int cand = CS->GetCandidate(i, sample[i]);
                seen_[cand] = false;
            }
        }
        if (!valid) return 0.0;
        for (Size i = 0; i < query_->GetNumVertices(); ++i) {
            sample[i] = CS->GetCandidate(i, sample[i]);
        }
        return 1.0;
    }

    int FaSTest::CheckSampleTree(std::vector<int> &sample) {
        for (auto &[i, j] : query_->all_edges) {
            if (sample[i] < 0 || sample[j] < 0) continue;
            if (!data_->CheckEdgeExist(sample[i], sample[j])) {
                return false;
            }
        }
        return true;
    }

    std::pair<double, int> FaSTest::UniformSamplingEstimate() {
        std::vector<int> tree_sample(query_->GetNumVertices(), -1);
        long long success = 0, t = 0;
        while (true) {
            std::fill(tree_sample.begin(), tree_sample.end(), -1);
            t++;
            if (GetSampleTree(tree_sample) == 1.0 and CheckSampleTree(tree_sample) == 1) {
                success++;
            }
            long double rhohat = (success * 1.0 / t);
            if (t >= 1000 && t % 100 == 0) {
                if (t == 50000 && success <= 10) {
                    return {-1, success};
                }
                long double wplus = boost::math::binomial_distribution<>::find_upper_bound_on_p(t, success, 0.05/2);
                long double wminus = boost::math::binomial_distribution<>::find_lower_bound_on_p(t, success, 0.05/2);
                if (rhohat * 0.8 < wminus && wplus < rhohat * 1.25) {
                    break;
                }
            }
        }
        return {total_trees_ * (success * 1.0 / t), success};
    }

    double FaSTest::EstimateEmbeddings() {
        Timer ts_timer, gs_timer;
        ts_timer.Start();
        std::pair<double, int> ts_result = UniformSamplingEstimate();
        ts_timer.Stop();
        if (ts_result.first < 0) {
            gs_timer.Start();
            double gs_result = GSSolver->IntersectionSamplingEstimate(ceil(opt.sample_size_K * query_->GetNumVertices() / sqrt(ts_result.second + 1)));
            gs_timer.Stop();
            return gs_result;
        }
        return ts_result.first;
    }

}