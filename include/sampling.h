#ifndef SUBGRAPHSAMPLER_SAMPLING_H
#define SUBGRAPHSAMPLER_SAMPLING_H
#include <utility>
#include <vector>
#include <random>

#include "global/global.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"
#include "include/candidate_space.h"
namespace daf {
class TreeSampling {
public:
    TreeSampling(const DataGraph &data, const QueryGraph &query, DAG &dag);
    ~TreeSampling();
    TreeSampling &operator=(const TreeSampling &) = delete;
    TreeSampling(const TreeSampling &) = delete;

    double EstimateEmbeddings(Size num_samples, bool HTSampling=false);
    double total_trees_ = 0;
    CandidateSpace CS;
    const DataGraph &data_;
    const QueryGraph &query_;
    DAG &dag_;

    std::vector<std::unordered_map<Vertex, double>> num_trees_;
    std::vector<std::vector<std::vector<std::vector<Vertex>>>> sample_candidates_;
    std::vector<std::vector<std::vector<std::vector<double>>>> sample_candidate_weights_;
    std::vector<std::vector<std::vector<std::discrete_distribution<int>>>> sample_dist_;
//    std::unordered_map<VertexPair, std::unordered_map<Vertex, std::vector<Vertex>>> sample_candidates_;
//    std::unordered_map<VertexPair, std::unordered_map<Vertex, std::vector<double>>> sample_candidate_weights_;
//    std::unordered_map<VertexPair, std::unordered_map<Vertex, std::discrete_distribution<int>>> sample_dist_;
    std::vector<Vertex> root_candidates_;
    std::discrete_distribution<int> sample_root_dist_;
    void ConstructTreeDP();
    double SampleCSTree(std::vector<Vertex> &sample, bool HTSampling=false);
    void BuildQueryTree();
};


}
#endif //SUBGRAPHSAMPLER_SAMPLING_H
