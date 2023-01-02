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
        DAG *dag_;
        Vertex root;

        RWI(){};
        ~RWI(){};
        void init (DataGraph &data, QueryGraph &query, DAG &dag, CandidateSpace &cs) {
            CS = &cs;
            data_ = &data;
            query_ = &query;
            dag_ = &dag;
            local_candidate_set_.resize(query_->GetNumVertices());
        }
        std::vector<tsl::robin_set<Vertex>> local_candidate_set_;
        std::vector<int> root_candidates_;
        double SampleCSDAGWithIntersection(std::vector<int> &sample);
        double IntersectionSamplingEstimate(Size num_samples);
        double SampleDAGVertex(std::vector<int> &dag_sample, int vertex_id);
    };

}
#endif //SUBGRAPHSAMPLER_RWI_H
