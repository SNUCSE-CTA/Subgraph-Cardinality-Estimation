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
        boost::dynamic_bitset<uint64_t> seen;

        RWI(){};
        ~RWI(){};
        void init (DataGraph *data, QueryGraph *query, DAG *dag, CandidateSpace *cs) {
            CS = cs;
            data_ = data;
            query_ = query;
            dag_ = dag;
            local_candidates_.resize(query_->GetNumVertices());
            local_candidate_set_.resize(query_->GetNumVertices());
            num_seen.resize(query_->GetNumVertices(), 0);
            best_neighbor.resize(query_->GetNumVertices(), -1);
            seen.resize(data_->GetNumVertices(), false);
        }
        std::vector<int> num_seen, best_neighbor;
        std::vector<tsl::robin_set<Vertex>> local_candidate_set_;
        std::vector<std::vector<Vertex>> local_candidates_;
        std::vector<std::pair<std::vector<Vertex>::iterator, std::vector<Vertex>::iterator>> iterators;
        std::vector<int> root_candidates_;
        double IntersectionSamplingEstimate(Size num_samples);
        double SampleDAGVertex(std::vector<int> &dag_sample, int vertex_id);

        void multivector_intersection(int index, bool debug);
    };

}
#endif //SUBGRAPHSAMPLER_RWI_H
