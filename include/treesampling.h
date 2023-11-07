#ifndef SUBGRAPHSAMPLER_SAMPLING_H
#define SUBGRAPHSAMPLER_SAMPLING_H
#include <utility>
#include <vector>
#include <random>

#include "global/global.h"
#include "include/daf_dag.h"
#include "include/daf_data_graph.h"
#include "include/daf_query_graph.h"
#include "include/daf_candidate_space.h"
#include "include/graphsampling.h"
namespace daf {

    class TreeSampling {
public:
    TreeSampling(DataGraph *data, Option opt);
    ~TreeSampling();
    TreeSampling &operator=(const TreeSampling &) = delete;
    TreeSampling(const TreeSampling &) = delete;

    std::pair<double, int> UniformSamplingEstimate();
    double EstimateEmbeddings();

    double total_trees_ = 0;
    CandidateSpace *CS;
    DataGraph *data_;
    QueryGraph *query_;
    DAG *dag_;
    GraphSampling *RWI_;
    Option opt;

    std::vector<bool> seen_, query_seen_;
//    std::vector<std::unordered_map<Vertex, double>> num_trees_;
    double **num_trees_;
    std::vector<std::vector<std::vector<std::vector<Vertex>>>> sample_candidates_;
    std::vector<std::vector<std::vector<std::vector<double>>>> sample_candidate_weights_;
    std::vector<std::vector<std::vector<std::discrete_distribution<int>>>> sample_dist_;


    std::vector<Vertex> root_candidates_;
    std::discrete_distribution<int> sample_root_dist_;
    double ConstructTreeDP();
    double SampleCSTree(std::vector<int> &sample);
    void BuildQueryTree();

    int CheckSample(std::vector<int> &sample);
    void RegisterQuery(QueryGraph *query, DAG *dag);

};


}
#endif //SUBGRAPHSAMPLER_SAMPLING_H
