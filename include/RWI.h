#ifndef SUBGRAPHSAMPLER_RWI_H
#define SUBGRAPHSAMPLER_RWI_H
#include "include/daf_candidate_space.h"
#pragma once
namespace daf {
    class RWI {
    public:
        CandidateSpace *CS;
        DataGraph *data_;
        QueryGraph *query_;
        Vertex root;
        bool* seen;
        int **local_candidates, *local_candidate_size;
        RWI(DataGraph *data){
            data_ = data;
            seen = new bool[data->GetNumVertices()];
            local_candidates = new int*[MAX_QUERY_VERTEX];
            for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
                local_candidates[i] = new int[data->GetNumVertices()];
            }
            local_candidate_size = new int[MAX_QUERY_VERTEX];
        };
        ~RWI(){
            for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
                delete[] local_candidates[i];
            }
            delete[] local_candidates;
            delete[] local_candidate_size;
            delete[] seen;
        };
        void init(QueryGraph *query, CandidateSpace *cs) {
            CS = cs;
            query_ = query;
            num_seen.resize(query_->GetNumVertices(), 0);
            memset(seen, 0, data_->GetNumVertices());
            memset(local_candidate_size, 0, query->GetNumVertices());
        }
        std::vector<int> num_seen;
        std::vector<std::pair<std::vector<Vertex>::iterator, std::vector<Vertex>::iterator>> iterators;
        std::vector<int> root_candidates_;
        double IntersectionSamplingEstimate(Size num_samples);
        double SampleDAGVertex(std::vector<int> &dag_sample, int vertex_id);

        void multivector_intersection(int index, bool debug);
    };

}
#endif //SUBGRAPHSAMPLER_RWI_H
