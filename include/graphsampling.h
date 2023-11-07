#ifndef GRAPHSAMPLING_H
#define GRAPHSAMPLING_H
#include "include/candidate_space.h"
#pragma once
namespace CardEst {
    class GraphSampling {
    public:
        Option opt;
        CandidateSpace *CS;
        DataGraph *data_;
        QueryGraph *query_;
        Vertex root;
        bool* seen;
        int **local_candidates, *local_candidate_size;
        GraphSampling(DataGraph *data, Option opt_){
            data_ = data;
            seen = new bool[data->GetNumVertices()];
            local_candidates = new int*[MAX_QUERY_VERTEX];
            for (int i = 0; i < MAX_QUERY_VERTEX; i++) {
                local_candidates[i] = new int[data->GetNumVertices()];
            }
            local_candidate_size = new int[MAX_QUERY_VERTEX];
            opt = opt_;
        };
        ~GraphSampling(){
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
            sample.resize(query_->GetNumVertices(), -1);
            memset(seen, 0, data_->GetNumVertices());
            memset(local_candidate_size, 0, query->GetNumVertices());
        }
        std::vector<int> num_seen, sample;
        std::vector<std::pair<std::vector<Vertex>::iterator, std::vector<Vertex>::iterator>> iterators;
        std::vector<int> root_candidates_;
        double IntersectionSamplingEstimate(Size num_samples);
        std::pair<double, int> EstimateWM(int vertex_id, int num_samples, double w);

        void multivector_intersection(int index, bool debug);

        int ChooseExtendableVertex(int vertex_id);

        void BuildExtendableCandidates(int u);
    };

}
#endif 