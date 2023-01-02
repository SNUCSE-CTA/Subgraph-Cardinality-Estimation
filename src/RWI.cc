#include "include/tsl/hopscotch_set.h"
#include "include/tsl/hopscotch_map.h"
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
    inline std::vector<int> sample_indices_without_replacement(int num_cands, int num_samples) {
        std::vector<int> indices;
        std::sample(population.begin(), population.begin()+num_cands,
                    std::back_inserter(indices),
                    num_samples, gen);
        std::sort(indices.begin(), indices.end());
        return indices;
    }
    inline std::vector<int> sample_without_replacement(tsl::robin_set<int> &cand, int num_samples) {
        std::vector<int> indices = sample_indices_without_replacement(cand.size(), num_samples);
        std::vector<int> result;
        auto it(cand.begin());
        for (int i = 0; i < indices.size(); i++) {
            if (i == 0) std::advance(it, indices[i]);
            else std::advance(it, indices[i] - indices[i-1]);
            result.push_back(*it);
        }
        return result;
    }

    int rwi_sample_count = 1000000;
    int local_cand_cnt = 0, local_cand_sum = 0;
    double RWI::SampleDAGVertex(std::vector<int> &dag_sample, int vertex_id) {
        // Vertex with minimum number of (1-edge) candidate
        std::vector<int> num_seen(query_->GetNumVertices(), 0);
        int u = -1;
        for (int i = 0; i < query_->GetNumVertices(); i++) {
            if (dag_sample[i] != -1) continue;
            for (int q_nbr : query_->adj_list[i]) {
                if (dag_sample[q_nbr] != -1) {
                    num_seen[i]++;
                }
            }
        }
        u = std::max_element(num_seen.begin(), num_seen.end()) - num_seen.begin();

        local_candidate_set_[u].clear();
        for (int q_nbr : query_->adj_list[u]) {
            if (dag_sample[q_nbr] == -1) continue;
            VertexPair uPair = {q_nbr, dag_sample[q_nbr]};
            if (local_candidate_set_[u].empty()) {
                local_candidate_set_[u] = CS->cs_adj_[uPair][u];
            }
            else {
                auto it1 = local_candidate_set_[u].begin();
                while (it1 != local_candidate_set_[u].end()) {
                    if (CS->cs_adj_[uPair][u].find(*it1) == CS->cs_adj_[uPair][u].end()) {
                        it1 = local_candidate_set_[u].erase(it1);
                    }
                    else ++it1;
                }
            }
        }
        for (int seen_candidate : dag_sample) {
            local_candidate_set_[u].erase(seen_candidate);
        }

        if (local_candidate_set_[u].empty()) {
            rwi_sample_count--;
            return 0;
        }
        if (vertex_id == query_->GetNumVertices()-1) {
            rwi_sample_count--;
            return local_candidate_set_[u].size() * 1.0;
        }
        // Branching sample
        double ht_s = 0.0;
        local_cand_sum += local_candidate_set_[u].size();
        local_cand_cnt += 1;

        int sample_space_size = local_candidate_set_[u].size();
        int num_branches = 1 + std::max(sample_space_size >> 5, std::min(sample_space_size, 4));
        num_branches = std::min(num_branches, sample_space_size);
        if (num_seen[u] == query_->adj_list[u].size()) num_branches = 1;
//        local_cand_sum += num_branches;
//        int num_branches = 1;

//        fprintf(stderr, "local_cand_size of vertexid %d (%d) = %lu\n", vertex_id, u, local_candidate_set_[u].size());
//        local_cand_sum += local_candidate_set_[u].size();
        local_cand_cnt++;
        std::vector <int> branches = sample_without_replacement(local_candidate_set_[u],num_branches);
        int skipped = 0;
        for (int b = 0; b < num_branches; b++) {
            if (rwi_sample_count <= 0) {
                skipped = num_branches - b;
                break;
            }
            dag_sample[u] = branches[b];
            auto nxt_result = SampleDAGVertex(dag_sample, vertex_id+1);
            ht_s += nxt_result;
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
//        std::cerr << "root : " << root << std::endl;
//        CS->printCS();
//        for (int x : root_candidates_) std::cerr << x << ' ';
//        std::cerr << std::endl;
        double ht_est = 0.0;
        int ht_count = 0;
        std::vector<int> dag_sample(query_->GetNumVertices(), -1);
        rwi_sample_count = num_samples;
        while (rwi_sample_count > 0) {
            for (int i = 0; i < query_->GetNumVertices(); i++)
                local_candidate_set_[i].clear();
            dag_sample[root] = root_candidates_[ht_count % root_candidates_.size()];
//            fprintf(stderr, "Matching %d-%d\n",dag_->GetRoot(),dag_sample[dag_->GetRoot()]);
//            fprintf(stderr, "%d samples left->",rwi_sample_count);
            auto recursion_result = SampleDAGVertex(dag_sample, 1);
//            fprintf(stderr, "%d samples left, weight = %.02lf\n",rwi_sample_count, recursion_result * root_candidates_.size());
            ht_count++;
            ht_est += recursion_result;
//            break;
        }
//        fprintf(stderr, "AVG Local Candset %.02lf\n",local_cand_sum*1.0/local_cand_cnt);
        ht_est *= root_candidates_.size();
        ht_est /= ht_count;
        return ht_est;
    }
}
