#include "include/graphsampling.h"
#include "global/timer.h"

using namespace CardEst;

int rwi_sample_count = -1;
long long local_cand_cnt = 0, local_cand_sum = 0;

std::map<int, int> mp;
void GraphSampling::multivector_intersection(int index, bool debug = false) {
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
                    ++iterators[i].first;
                }
                else if (*iterators[i].first > target) {
                    goto nxt_target;
                }
                else break;
            }
            if (iterators[i].first == iterators[i].second) return;
        }
        local_candidates[index][local_candidate_size[index]++] = target;
        nxt_target:
        ++iterators[0].first;
    }
}

int GraphSampling::ChooseExtendableVertex(int vertex_id) {
    int u = -1;
    int max_open_neighbors = 0;
    int min_nbr_cnt = 1e9;
    for (int i = 0; i < query_->GetNumVertices(); i++) {
        if (sample[i] != -1) continue;
        int nbr_cnt = 1e9;
        int open_neighbors = 0;
        for (int q_nbr : query_->adj_list[i]) {
            if (sample[q_nbr] != -1) {
                open_neighbors++;
                int num_nbr = CS->cs_edge_[q_nbr][sample[q_nbr]][i].size();
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
    return u;
}

void GraphSampling::BuildExtendableCandidates(int u) {
    local_candidate_size[u] = 0;
    iterators.clear();
    for (int q_nbr : query_->adj_list[u]) {
        if (sample[q_nbr] == -1) continue;
        iterators.emplace_back(CS->cs_edge_[q_nbr][sample[q_nbr]][u].begin(), CS->cs_edge_[q_nbr][sample[q_nbr]][u].end());
    }
    std::sort(iterators.begin(), iterators.end(), [](auto &a, auto &b) -> bool {
        return a.second - a.first < b.second - b.first;
    });
    multivector_intersection(u);
}


std::pair<double, int> GraphSampling::EstimateWM(int vertex_id, int num_samples, double w) {
    int u = ChooseExtendableVertex(vertex_id);
    BuildExtendableCandidates(u);
    if (local_candidate_size[u] == 0) {
        return {0, 1};
    }
    for (int i = 0; i < query_->GetNumVertices(); i++) {
        if (sample[i] == -1) continue;
        seen[CS->GetCandidate(i, sample[i])] = true;
    }

    for (int i = 0; i < local_candidate_size[u]; ++i) {

        if (seen[CS->GetCandidate(u, local_candidates[u][i])]) {
            local_candidates[u][i] = local_candidates[u][local_candidate_size[u]-1];
            local_candidate_size[u]--;
            i--;
        }
    }

    for (int i = 0; i < query_->GetNumVertices(); i++) {
        if (sample[i] == -1) continue;
        seen[CS->GetCandidate(i, sample[i])] = false;
    }

    if (local_candidate_size[u] == 0) {
        local_candidate_size[u] = 0;
        sample[u] = -1;
        return {0, 1};
    }
    if (vertex_id == query_->GetNumVertices()-1) {
        double return_value = local_candidate_size[u] * 1.0;
        local_candidate_size[u] = 0;
        sample[u] = -1;
        return {w * return_value, 1};
    }
    
    local_cand_sum += local_candidate_size[u];
    local_cand_cnt += 1;
    mp[vertex_id]++;

    int sample_space_size = local_candidate_size[u];

    int num_branches = std::min(std::max(sample_space_size/4, 4), num_samples);
    num_branches = std::min(num_branches, sample_space_size);

    int i = 0;
    int num_used = 0;
    double est = 0.0;
    while (num_used < num_samples and local_candidate_size[u] > 0) {
        int idx = gen()%local_candidate_size[u];
        sample[u] = local_candidates[u][idx];
        int num_next_samples = (num_samples - num_used) / std::max(num_branches-i, 1);
        if (num_next_samples == 0) num_next_samples = (num_samples - num_used);
        double est_; int num_used_;
        std::tie(est_, num_used_) = EstimateWM(vertex_id + 1, num_next_samples, w * sample_space_size);
        est += est_;
        num_used += num_used_;
        local_candidates[u][idx] = local_candidates[u][local_candidate_size[u]-1];
        local_candidate_size[u]--;
        sample[u] = -1;
        i++;
        if (i == num_branches) break;
    }
    sample[u] = -1;
    local_candidate_size[u] = 0;
    return {est / i, num_used};
}

double GraphSampling::IntersectionSamplingEstimate(Size num_samples) {
    local_cand_sum = local_cand_cnt = 0;
    
    std::vector <int> num_cands(query_->GetNumVertices());
    for (int i = 0; i < query_->GetNumVertices(); i++) {
        num_cands[i] = CS->GetCandidateSetSize(i);
    }

    root = std::min_element(num_cands.begin(), num_cands.end()) - num_cands.begin();
    root_candidates_ = CS->candidate_set_[root];
    std::shuffle(root_candidates_.begin(), root_candidates_.end(), gen);
    double est = 0.0;
    int num_roots = 0;
    rwi_sample_count = num_samples;
    int num_root_samples = std::min((int)root_candidates_.size(), (int)rwi_sample_count);
    num_root_samples = std::min(num_root_samples, (int)root_candidates_.size());
    int used_samples = 0;
    while (used_samples < rwi_sample_count and num_roots < num_root_samples) {
        std::fill(sample.begin(), sample.end(), -1);
        memset(local_candidate_size, 0, query_->GetNumVertices());
        sample[root] = (num_roots % root_candidates_.size());
        int num_sample_use = (rwi_sample_count - used_samples) / (std::max(num_root_samples - num_roots, 1));
        auto recursion_result = EstimateWM(1, num_sample_use, 1.0 * root_candidates_.size());
        est += recursion_result.first;
        used_samples += recursion_result.second;
        num_roots++;
    }
    est /= num_roots;
    return est;
}