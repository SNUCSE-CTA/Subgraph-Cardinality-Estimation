#include "include/RWI.h"
#include "global/timer.h"

using namespace daf;
namespace daf {

    int rwi_sample_count = 1000000;
    long long local_cand_cnt = 0, local_cand_sum = 0;

    std::map<int, int> mp;
    void RWI::multivector_intersection(int index, bool debug = false) {
        if (local_candidate_size[index] > 0) return;
        int num_vectors = iterators.size();
        if (num_vectors == 1) {
            while (iterators[0].first != iterators[0].second) {
                local_candidates[index][local_candidate_size[index]++] = (*iterators[0].first);
                ++iterators[0].first;
            }
            return;
        }
        while (iterators[0].first != iterators[0].second) {
            int target = *iterators[0].first;
            for (int i = 1; i < num_vectors; i++) {
                while (iterators[i].first != iterators[i].second) {
                    if (*iterators[i].first < target) {
//                        fprintf(stderr, "Increment vector %d from %d\n",i,*iterators[i].first);
                        ++iterators[i].first;
                    }
                    else if (*iterators[i].first > target) {
//                        fprintf(stderr, "Dead by vector %d : go to nxt target\n",i);
                        goto nxt_target;
                    }
                    else break;
                }
                if (iterators[i].first == iterators[i].second) return;
            }
//            fprintf(stderr, "Found %d in all vectors\n",target);
            local_candidates[index][local_candidate_size[index]++] = target;
            nxt_target:
            ++iterators[0].first;
        }
    }

    int RWI::ChooseExtendableVertex(int vertex_id) {
        int u = -1;
        if (opt.sampling_order == OPENNEIGHBORS) {
            int max_open_neighbors = 0;
            int min_nbr_cnt = 1e9;
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (dag_sample[i] != -1) continue;
                int nbr_cnt = 1e9;
                int open_neighbors = 0;
                for (int q_nbr : query_->adj_list[i]) {
                    if (dag_sample[q_nbr] != -1) {
                        open_neighbors++;
                        int num_nbr = CS->cs_edge_[q_nbr][dag_sample[q_nbr]][i].size();
                        if (num_nbr < nbr_cnt) {
                            nbr_cnt = num_nbr;
                        }
                    }
                }
                if (open_neighbors > max_open_neighbors) {
                    max_open_neighbors = open_neighbors;
                    min_nbr_cnt = nbr_cnt;
                    u = i;
                }
                else if (open_neighbors == max_open_neighbors) {
                    if (nbr_cnt < min_nbr_cnt) {
                        min_nbr_cnt = nbr_cnt;
                        u = i;
                    }
                }
            }
        }
        else if (opt.sampling_order == APPROXEXTCAND) {
            // Vertex with minimum number of (1-edge) candidate
            double min_value = 1e9;
            for (int i = 0; i < query_->GetNumVertices(); i++) {
                if (dag_sample[i] != -1) continue;
                double nbr_cnt = 1e9;
                int nbr_seen = 0;
                for (int q_nbr : query_->adj_list[i]) {
                    if (dag_sample[q_nbr] != -1) {
                        nbr_seen += 1;
                        int num_nbr = CS->cs_edge_[q_nbr][dag_sample[q_nbr]][i].size();
                        if (num_nbr < nbr_cnt) {
                            nbr_cnt = num_nbr;
                        }
                    }
                }
                if (nbr_seen == query_->adj_list[i].size()) {
                    if (local_candidate_size[i] == 0) {
                        BuildExtendableCandidates(i);
                    }
                    nbr_cnt = local_candidate_size[i];
                }
                if (nbr_cnt < min_value) {
                    min_value = nbr_cnt;
                    u = i;
                }
            }
        }
        return u;
    }

    void RWI::BuildExtendableCandidates(int u) {
        local_candidate_size[u] = 0;
        iterators.clear();
        for (int q_nbr : query_->adj_list[u]) {
            if (dag_sample[q_nbr] == -1) continue;
            iterators.emplace_back(CS->cs_edge_[q_nbr][dag_sample[q_nbr]][u].begin(), CS->cs_edge_[q_nbr][dag_sample[q_nbr]][u].end());
        }
        std::sort(iterators.begin(), iterators.end(), [](auto &a, auto &b) -> bool {
            return a.second - a.first < b.second - b.first;
        });
        multivector_intersection(u);
    }
    std::pair<double, int> RWI::SampleDAGVertex(int vertex_id, int num_samples, double w) {

        int u = ChooseExtendableVertex(vertex_id);

        // GetExtendableCandidates(u);
        BuildExtendableCandidates(u);
        if (local_candidate_size[u] == 0) {
            return {0, 1};
        }
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            if (dag_sample[i] == -1) continue;
            seen[CS->GetCandidate(i, dag_sample[i])] = true;
        }

        for (int i = 0; i < local_candidate_size[u]; ++i) {
//            fprintf(stderr, "i = %d, local_candidates[%d][i] = %d\n",i,u, local_candidates[u][i]);
            if (seen[CS->GetCandidate(u, local_candidates[u][i])]) {
                local_candidates[u][i] = local_candidates[u][local_candidate_size[u]-1];
                local_candidate_size[u]--;
                i--;
            }
        }

        for (int i = 0; i < query_->GetNumVertices(); i++) {
            if (dag_sample[i] == -1) continue;
            seen[CS->GetCandidate(i, dag_sample[i])] = false;
        }

        // Extendable Candidates OK
        if (local_candidate_size[u] == 0) {
            local_candidate_size[u] = 0;
            dag_sample[u] = -1;
            return {0, 1};
        }
        if (vertex_id == query_->GetNumVertices()-1) {
            double return_value = local_candidate_size[u] * 1.0;
            local_candidate_size[u] = 0;
            dag_sample[u] = -1;
            return {w * return_value, 1};
        }
        // Branching sample
        local_cand_sum += local_candidate_size[u];
        local_cand_cnt += 1;
        mp[vertex_id]++;

        int sample_space_size = local_candidate_size[u];
//        int num_branches = 1 + std::max(sample_space_size >> 3, std::min(sample_space_size, 4));
        int num_branches = std::min(std::max((int)sqrt(sample_space_size), 4), num_samples);
        num_branches = std::min(num_branches, sample_space_size);
//        if (num_seen[u] == query_->adj_list[u].size()) num_branches = 1;
//        if (vertex_id == 1) {
//            fprintf(stderr, "sample_space_size = %d, num_branches = %d\n", sample_space_size, num_branches);
//        }
        int i = 0;
        int num_used = 0, last_used = 0;
        double est = 0.0;
        while (num_used < num_samples and local_candidate_size[u] > 0) {
            int idx = gen()%local_candidate_size[u];
            dag_sample[u] = local_candidates[u][idx];
            int num_next_samples = (num_samples - num_used) / std::max(num_branches-i, 1);
            if (num_next_samples == 0) num_next_samples = (num_samples - num_used);
//            if (num_next_samples < last_used / 10) break;
            double est_; int num_used_;
//            printf("Try Recursion to %d with %d samples\n", vertex_id+1, num_next_samples);
            std::tie(est_, num_used_) = SampleDAGVertex(vertex_id+1, num_next_samples, w * sample_space_size);
            est += est_;
            num_used += num_used_;
            local_candidates[u][idx] = local_candidates[u][local_candidate_size[u]-1];
            local_candidate_size[u]--;
            dag_sample[u] = -1;
            i++;
            if (i == num_branches) break;
        }
        dag_sample[u] = -1;
        local_candidate_size[u] = 0;
        return {est / i, num_used};
    }

    double RWI::IntersectionSamplingEstimate(Size num_samples) {
        local_cand_sum = local_cand_cnt = 0;
        // Choose the vertex with minimum C(u) as root
        std::vector <int> num_cands(query_->GetNumVertices());
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            num_cands[i] = CS->GetCandidateSetSize(i);
        }
        root = std::min_element(num_cands.begin(), num_cands.end()) - num_cands.begin();
        root_candidates_ = CS->candidate_set_[root];
        std::shuffle(root_candidates_.begin(), root_candidates_.end(), gen);
        double ht_est = 0.0;
        int ht_count = 0;
        rwi_sample_count = num_samples;
        fprintf(stderr, "num_samples = %d\n", num_samples);
        int num_root_samples = (root_candidates_.size()), last_used = 0;
//        num_root_samples = std::min(num_root_samples, 512);
//        num_root_samples = std::max(num_root_samples, 32);
        num_root_samples = std::min(num_root_samples, (int)root_candidates_.size());
        fprintf(stderr, "num_root_samples = %d (root_cand_size = %d)\n", num_root_samples, (int)root_candidates_.size());
        int used_samples = 0;
        while (used_samples < rwi_sample_count) {
            std::fill(dag_sample.begin(), dag_sample.end(), -1);
            memset(local_candidate_size, 0, query_->GetNumVertices());
            dag_sample[root] = (ht_count % root_candidates_.size());
            int num_sample_use = (rwi_sample_count - used_samples) / (std::max(num_root_samples-ht_count, 1));
//            if (num_sample_use < last_used / 10) break;
//            printf("Try to use %d samples...\n", num_sample_use);
            auto recursion_result = SampleDAGVertex(1, num_sample_use, 1.0 * root_candidates_.size());
//            printf("Recursion result: %lf %d\n", recursion_result.first, recursion_result.second);
            ht_est += recursion_result.first;
            used_samples += recursion_result.second;
            last_used = num_sample_use;
            ht_count++;
            if (ht_count == num_root_samples) break;
        }
        fprintf(stdout, "#Candset : %.02lf\n",local_cand_sum*1.0/local_cand_cnt);
//        ht_est *= root_candidates_.size();
        ht_est /= ht_count;
//        fprintf(stdout, "Total %d nodes explored\n",local_cand_cnt);
//        for (auto &[u, v] : mp) {
//            fprintf(stderr, "Instances of sample space %d = %d\n",u,v);
//        }
        return ht_est;
    }
}