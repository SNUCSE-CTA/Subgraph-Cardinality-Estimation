#ifndef TREESAMPLING_H
#define TREESAMPLING_H
#include <utility>
#include <vector>
#include <random>

#include "global/global.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"
#include "include/candidate_space.h"
#include "include/graphsampling.h"

namespace CardEst {
class FaSTest {
public:
    FaSTest(DataGraph *data, Option opt);
    ~FaSTest();
    FaSTest &operator=(const FaSTest &) = delete;
    FaSTest(const FaSTest &) = delete;

    std::pair<double, int> UniformSamplingEstimate();
    double EstimateEmbeddings();

    double total_trees_ = 0;
    CandidateSpace *CS;
    DataGraph *data_;
    QueryGraph *query_;
    OrderedQueryGraph *dag_;
    GraphSampling *GSSolver;
    Option opt;

    std::vector<bool> seen_, query_seen_;
    double **num_trees_;
    std::vector<std::vector<std::vector<std::vector<Vertex>>>> sample_candidates_;
    std::vector<std::vector<std::vector<std::vector<double>>>> sample_candidate_weights_;
    std::vector<std::vector<std::vector<std::discrete_distribution<int>>>> sample_dist_;

    std::vector<Vertex> root_candidates_;
    std::discrete_distribution<int> sample_root_dist_;
    double CountCandidateTrees();
    double GetSampleTree(std::vector<int> &sample);
    void BuildQueryTree();
    int CheckSampleTree(std::vector<int> &sample);
    void PrepareQuery(QueryGraph *query, OrderedQueryGraph *dag);

};
}
#endif 
