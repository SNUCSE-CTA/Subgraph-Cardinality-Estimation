#include <algorithm>
#include <iostream>
#include <string>
#include <iomanip>
#include <vector>

#include "global/timer.h"
#include "include/candidate_space.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"
#include "include/sampling.h"

const int UNIFORMRANDOM = 1;
const int HORVITZTHOMPSON = 2;

Timer total_timer, sample_timer;

std::string data_name = "../../dataset/wordnet/data_graph/wordnet.graph";
std::string query_name = "../../dataset/wordnet/query_graph/query_dense_4_1.graph";
uint32_t estimator = UNIFORMRANDOM;
uint64_t limit = std::numeric_limits<uint64_t>::max();
int num_samples = 1000000;

void run_treesample (DataGraph &data, QueryGraph &query) {
    total_timer.Start();
    query.LoadAndProcessGraph(data);
    total_timer.Stop();

    daf::DAG dag(data, query);

    total_timer.Start();
    dag.BuildDAG();
    total_timer.Stop();

    total_timer.Start();
    daf::TreeSampling treesampling(data, query, dag);
    total_timer.Stop();

    for (auto u = 0; u < query.GetNumVertices(); ++u) {
        if (treesampling.CS.GetCandidateSetSize(u) == 0) {
            std::cout << "Total time: " << total_timer.GetTime() << " ms\n";
            exit(1);
        }
    }
    sample_timer.Start();
    double est = treesampling.EstimateEmbeddings(num_samples, (estimator == HORVITZTHOMPSON));
    sample_timer.Stop();

    total_timer.Add(sample_timer);

    std::cout << "Total time: " << total_timer.GetTime() << " ms\n";
    std::cout << "Search time: " << sample_timer.GetTime() << " ms\n";
    std::cout << "#Matches(Approx) : " << std::fixed << std::setprecision(6) << est << std::endl;
    std::cout << "#Tree : " << std::fixed << std::setprecision(6) << treesampling.total_trees_ << std::endl;
}


int get_estimator(char *est_name) {
    if (strcmp(est_name, "uniform-random") == 0) return UNIFORMRANDOM;
    if (strcmp(est_name, "horvitz-thompson") == 0) return HORVITZTHOMPSON;
    return UNIFORMRANDOM;
}
int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'd':
                    data_name = argv[i + 1];
                    break;
                case 'q':
                    query_name = argv[i + 1];
                    break;
                case 'm':
                    limit = std::atoi(argv[i + 1]);
                case 's':
                    num_samples = std::atoi(argv[i + 1]);
                case 't':
                    estimator = get_estimator(argv[i + 1]);
                    break;

            }
        }
    }

    std::cout << "Loading data graph...\n";
    DataGraph data(data_name);
    data.LoadAndProcessGraph();

    std::cout << "Loading query graph...\n";
    QueryGraph query(query_name);
    run_treesample(data, query);
}
