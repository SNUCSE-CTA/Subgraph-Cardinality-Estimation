#include <stack>
#include <queue>
#include "include/treesampling.h"
#include "global/global.h"
#include "global/timer.h"

struct BridgeFinder {
    int n;
    std::vector<std::pair<Vertex, Vertex>> bridges;
    std::vector<std::vector<std::pair<double, int>>> *G;
    int low[100], order[100];
    bool vst[100];
    int t = 0, rt = 0;
    void dfs(int r, int p) {
        vst[r] = true;
        order[r] = t++;
        low[r] = t;
        for (auto &[w, nxt] : (*G)[r]) {
            if (nxt == p) continue;
            if (!vst[nxt]) {
                dfs(nxt, r);
                if (low[nxt] > order[r])
                    bridges.emplace_back(std::min(r, nxt), std::max(r, nxt));
                low[r] = std::min(low[nxt], low[r]);
            }
            else low[r] = std::min(low[r], order[nxt]);
        }
    }
    std::vector<std::pair<Vertex, Vertex>> FindBridges(std::vector<std::vector<std::pair<double, int>>> *g, int _n) {
        G = g;
        n = _n;
        memset(low, 0, sizeof(low));
        memset(order, 0, sizeof(order));
        memset(vst, 0, sizeof(vst));
        bridges.clear();
        t = rt = 0;
        for (int i = 0; i < n; i++)
            if (!vst[i])
                dfs(i, -1);
        sort(bridges.begin(), bridges.end());
        bridges.erase(unique(bridges.begin(), bridges.end()), bridges.end());
        return bridges;
    }
} BF;

namespace daf {
    inline Size sample(std::discrete_distribution<int> &weighted_distr) {
        return weighted_distr(gen);
    }

    TreeSampling::TreeSampling(DataGraph *data, FilterOption opt) {
        data_ = data;
        if (!data_->is_sparse()) {
            opt.structure_filter = NO_STRUCTURE_FILTER;
        }
        opt.print();
        CS = new CandidateSpace(data, opt);
        seen_.resize(data_->GetNumVertices());
        RWI_ = new RWI(data);
        num_trees_ = new double*[MAX_QUERY_VERTEX];
        for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
            num_trees_[i] = new double[data_->GetNumVertices()];
        }
    }

    void TreeSampling::RegisterQuery(QueryGraph *query, DAG *dag) {
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
        std::cout << "Query Tree Building Time: " << querytree_timer.GetTime() << " ms\n";
        RWI_->init(query, CS);
    }

    TreeSampling::~TreeSampling() {
        for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
            delete[] num_trees_[i];
        }
        delete[] num_trees_;
    }

    static double pairwise_corr = 0.0;
    static int num_pairs = 0;

    std::vector<std::pair<Vertex, Vertex>> BuildTreeByBFS(std::vector<std::pair<double, std::pair<Vertex, Vertex>>> &edges) {
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
        std::vector<int> visited(qV, 0);
        int root = 0;
        std::queue<int> q;
        q.push(root);
        visited[root] = true;
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (auto &[w, v] : G[u]) {
                if (!visited[v]) {
//                    fprintf(stderr, "%d %d with weight = %.08lf\n", u, v, w);
                    T.emplace_back(u, v);
                    visited[v] = true;
                    q.push(v);
                }
            }
        }
        return T;
    }
    std::vector<std::pair<Vertex, Vertex>> BuildTreeByDFS(std::vector<std::pair<double, std::pair<Vertex, Vertex>>> &edges) {
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
        std::vector<int> visited(qV, 0);
        int root = 0;
        std::stack<std::pair<int, int>> stk;
        stk.emplace(root, 0);
        visited[root] = true;
        while (!stk.empty()) {
            auto [u, s] = stk.top();
            stk.pop();
            for (int i = s; i < G[u].size(); i++) {
                auto [w, v] = G[u][i];
                if (!visited[v]) {
                    fprintf(stderr, "%d %d with weight = %.08lf\n", u, v, w);
                    T.emplace_back(u, v);
                    visited[v] = true;
                    stk.emplace(u, i+1);
                    stk.emplace(v, 0);
                    break;
                }
            }
        }
        return T;
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
//        auto bridges = BF.FindBridges(&G, qV);
        std::vector<int> deg(qV, 0);
        UnionFind uf(qV);
//        for (auto &[u, v] : bridges) {
////            fprintf(stderr, "Add %d %d is bridge\n", u, v);
//            T.emplace_back(u, v);
//            deg[u]++;
//            deg[v]++;
//            uf.unite(u, v);
//        }
//        int root = 0;
//        visited[root] = true;
        while (T.size() + 1 < qV) {
            double minw = 1e9;
            std::pair<Vertex, Vertex> me;
            for (auto &[w, e] : edges) {
                auto [u, v] = e;
                if (uf.find(u) == uf.find(v)) continue;
//                fprintf(stderr, "    Considering %d %d with weight = %.08lf\n", u, v, w);
                if (minw > w + deg[u] * 1e-7) {
                    minw = w + deg[u] * 1e-7;
                    me = e;
                }
            }
//            fprintf(stderr, "Add %d %d with weight = %.08lf\n", me.first, me.second, minw);
            uf.unite(me.first, me.second);
            T.push_back(me);
            deg[me.first]++;
            deg[me.second]++;
        }
        return T;
    }

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
        // Tree building via MST
//        UnionFind uf(query_->GetNumVertices());
//        for (auto e : edges) {
//            auto [u, v] = e.second;
//
//            fprintf(stderr, "%d %d with weight = %.08lf\n", u, v, e.first);
//            if (uf.unite(u, v)) {
//                dag_->AddTreeEdge(u, v);
//            }
//        }
//        dag_->BuildTree();

//         Tree building via BFS/DFS Tree
        auto T = BuildTreeByPrim(edges);
        for (auto e : T) {
            auto [u, v] = e;
            dag_->AddTreeEdge(u, v);
        }
        dag_->BuildTree();

        double num_tree_homo = ConstructTreeDP();
        fprintf(stderr, "NUM_TREE = %.04lE\n",num_tree_homo);
        // Compute Pairwise Correlation Statistics
        /*
        for (int u = 0; u < query_->GetNumVertices(); u++) {
            if (query_->adj_list[u].size() <= 1) continue;
            if (CS->GetCandidateSetSize(u) <= 10) continue;
            std::vector<std::vector<double>> degree_sequences(query_->adj_list[u].size());
            std::vector<double> means, stdevs;
            int uc_idx = 0;
            for (int u1 : query_->adj_list[u]) {
                double mean = 0.0;
                for (int v_idx = 0; v_idx < CS->GetCandidateSetSize(u); v_idx++) {
                    degree_sequences[uc_idx].push_back(CS->cs_edge_[u][v_idx][u1].size() * 1.0);
                    mean += degree_sequences[uc_idx].back();
                }
                mean /= CS->GetCandidateSetSize(u);
                double sample_stdev = 0.0;
                for (double &x : degree_sequences[uc_idx]) {
                    sample_stdev += (x - mean) * (x - mean);
                }
                sample_stdev /= (1.0 * (CS->GetCandidateSetSize(u) - 1));
                sample_stdev = sqrt(sample_stdev);
                means.push_back(mean);
                stdevs.push_back(sample_stdev);
                if (sample_stdev > 1e-6) {
                    for (int v_idx = 0; v_idx < CS->GetCandidateSetSize(u); v_idx++) {
                        degree_sequences[uc_idx][v_idx] = (degree_sequences[uc_idx][v_idx] - mean) / sample_stdev;
                    }
                }
//                fprintf(stderr, "u = %d, Neighbor %d : ", u, u1);
//                for (int v_idx = 0; v_idx < CS->GetCandidateSetSize(u); v_idx++) {
//                    fprintf(stderr, " %02lf", degree_sequences[uc_idx][v_idx]);
//                }
//                fprintf(stderr, "\n");
                uc_idx++;
            }
            for (int i = 0; i < uc_idx; i++) {
                for (int j = i + 1; j < uc_idx; j++) {
                    std::vector<double> &D1 = degree_sequences[i], &D2 = degree_sequences[j];
                    double c = 0.0;
                    if (stdevs[i] <= 1e-6 || stdevs[j] <= 1e-6) {
                        c = 0.0;
                    }
                    else {
                        for (int idx = 0; idx < D1.size(); idx++) {
                            c += D1[idx] * D2[idx];
                        }
                        c /= (D1.size() - 1);
//                        if (c >= 0.999) {
//                            for (int v_idx = 0; v_idx < CS->GetCandidateSetSize(u); v_idx++) {
//                                fprintf(stderr, " %.02lf", D1[v_idx]);
//                            }
//                            fprintf(stderr, "\n");
//                            for (int v_idx = 0; v_idx < CS->GetCandidateSetSize(u); v_idx++) {
//                                fprintf(stderr, " %.02lf", D2[v_idx]);
//                            }
//                            fprintf(stderr, "\n");
//                        }
                    }
                    fprintf(stderr, "u = %d, corr[%d, %d] = %.04lf\n",u,query_->adj_list[u][i],query_->adj_list[u][j],c);
                    pairwise_corr += abs(c);
                    num_pairs += 1;
                }
            }
        }
        fprintf(stderr, "Pairwise absolute pearson correlation = %.04lf among %d pairs\n", pairwise_corr / (1.0 *num_pairs), num_pairs);
        */
    }

    double TreeSampling::ConstructTreeDP() {
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
            if (t >= 1000 && t % 100 == 0) {
                if (t == 50000 && success * 5000 <= t) {
                    fprintf(stderr, "NUM_HOMO  %d\tNUM_CFLCT %d\n", reject_homo, reject_nontree);
                    fprintf(stderr, "#NUM_SAMPLES : %lld, #NUM_SUCCESS : %lld\n", t, success);
                    fprintf(stdout, "#NUM_SAMPLES : %lld\n", t);
                    fprintf(stdout, "#NUM_SUCCESS : %lld\n", success);
                    return {-1, success};
                }
                /*
                // Wilson Confidence Interval
                long double wminus = (success + z * z / 2 - z * sqrt(success * (t - success) / (long double) t + z * z / 4)) / (t + z * z);
                long double wplus = (success + z * z / 2 + z * sqrt(success * (t - success) / (long double) t + z * z / 4)) / (t + z * z);
                */

                // Clopper-Pearson bounds
                long double wplus = boost::math::binomial_distribution<>::find_upper_bound_on_p(t, success, 0.05/2);
                long double wminus = boost::math::binomial_distribution<>::find_lower_bound_on_p(t, success, 0.05/2);
                
                if (rhohat * 0.8 < wminus && wplus < rhohat * 1.25) {
                    break;
                }
                /*
                // Jeffreys bounds
                long double wplus = boost::math::binomial_distribution<>::find_upper_bound_on_p(t, success, 0.05/2, boost::math::binomial_distribution<>::jeffreys_prior_interval);
                long double wminus = boost::math::binomial_distribution<>::find_lower_bound_on_p(t, success, 0.05/2, boost::math::binomial_distribution<>::jeffreys_prior_interval);
                assert(wplus >= wminus);
                
                if (rhohat * 0.8 < wminus && wplus < rhohat * 1.25) {
                    break;
                }
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
//        return 0.0;
        Timer sampletimer_uni, sampletimer_inter;
        sampletimer_uni.Start();
        std::pair<double, int> uniformResult = UniformSamplingEstimate();
        sampletimer_uni.Stop();
        std::cout << "Uniform Sampling time: " << std::fixed << sampletimer_uni.GetTime() << " ms\n";
        if (uniformResult.first < 0) {
            sampletimer_inter.Start();
            double intersectionResult = RWI_->IntersectionSamplingEstimate(ceil(50000 * query_->GetNumVertices() / sqrt(uniformResult.second + 1)));
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