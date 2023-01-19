#include "include/RWI.h"
#include "global/timer.h"

using namespace daf;
namespace daf {

    std::vector<int> population;

    inline Size sample_without_weight(std::vector<int> &cand) {
        int idx = gen() % cand.size();
        return cand[idx];
    }

    inline Size sample_without_weight(std::set<int> &cand) {
        int idx = gen() % cand.size();
        std::set<int>::const_iterator it(cand.begin());
        advance(it, idx);
        return *it;
    }

    int rwi_sample_count = 1000000;
    int local_cand_cnt = 0, local_cand_sum = 0;

    void RWI::multivector_intersection(int index, bool debug = false) {
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

    double RWI::SampleDAGVertex(std::vector<int> &dag_sample, int vertex_id) {
        // Vertex with minimum number of (1-edge) candidate
        std::fill(num_seen.begin(), num_seen.end(), 0);
        int u = -1;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            if (dag_sample[i] != -1) continue;
            int nbr_cnt = 1e9;
            for (int q_nbr : query_->adj_list[i]) {
                if (dag_sample[q_nbr] != -1) {
                    num_seen[i]++;
                    int num_nbr = CS->cs_edge_[q_nbr][dag_sample[q_nbr]][i].size();
                    if (num_nbr < nbr_cnt) {
                        nbr_cnt = num_nbr;
                    }
                }
            }
        }
        u = std::max_element(num_seen.begin(), num_seen.end()) - num_seen.begin();
        local_candidate_size[u] = 0;
        iterators.clear();
        for (int q_nbr : query_->adj_list[u]) {
            if (dag_sample[q_nbr] == -1) continue;
            iterators.push_back({CS->cs_edge_[q_nbr][dag_sample[q_nbr]][u].begin(), CS->cs_edge_[q_nbr][dag_sample[q_nbr]][u].end()});
        }
        std::sort(iterators.begin(), iterators.end(), [](auto &a, auto &b) -> bool {
            return a.second - a.first < b.second - b.first;
        });
//        std::cerr << "INTERSECTION" << std::endl;
        multivector_intersection(u);
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

        if (local_candidate_size[u] == 0) {
            rwi_sample_count--;
            return 0;
        }
        if (vertex_id == query_->GetNumVertices()-1) {
            rwi_sample_count--;
            return local_candidate_size[u] * 1.0;
        }
        // Branching sample
        double ht_s = 0.0;
        local_cand_sum += local_candidate_size[u];
        local_cand_cnt += 1;

        int sample_space_size = local_candidate_size[u];
        int num_branches = 1 + std::max(sample_space_size >> 5, std::min(sample_space_size, 4));
        num_branches = std::min(num_branches, sample_space_size);
        if (num_seen[u] == query_->adj_list[u].size()) num_branches = 1;

        int skipped = 0;
        if (num_branches == 1) {
            dag_sample[u] = local_candidates[u][gen()%local_candidate_size[u]];
            auto nxt_result = SampleDAGVertex(dag_sample, vertex_id+1);
            ht_s += nxt_result;
        }
        else {
            if (num_branches < local_candidate_size[u]) {
                std::shuffle(local_candidates[u], local_candidates[u] + local_candidate_size[u], gen);
            }
            for (int b = 0; b < num_branches; b++) {
                if (rwi_sample_count <= 0) {
                    skipped = num_branches - b;
                    break;
                }
                dag_sample[u] = local_candidates[u][b];
                auto nxt_result = SampleDAGVertex(dag_sample, vertex_id+1);
                ht_s += nxt_result;
            }
        }
        dag_sample[u] = -1;
        return (sample_space_size * ht_s / (num_branches - skipped));

    }




    double RWI::IntersectionSamplingEstimate(Size num_samples) {
        population.resize(data_->GetNumVertices());
        for (int i = 0; i < data_->GetNumVertices(); i++) population[i] = i;
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
        std::vector<int> dag_sample(query_->GetNumVertices(), -1);
        rwi_sample_count = num_samples;
        while (rwi_sample_count > 0) {
            memset(local_candidate_size, 0, query_->GetNumVertices());
//            dag_sample[root] = root_candidates_[ht_count % root_candidates_.size()];
            dag_sample[root] = (ht_count % root_candidates_.size());
            auto recursion_result = SampleDAGVertex(dag_sample, 1);
            ht_count++;
            ht_est += recursion_result;
//            break;
        }
        fprintf(stderr, "AVG Local Candset %.02lf\n",local_cand_sum*1.0/local_cand_cnt);
        ht_est *= root_candidates_.size();
        ht_est /= ht_count;
        return ht_est;
    }
}
