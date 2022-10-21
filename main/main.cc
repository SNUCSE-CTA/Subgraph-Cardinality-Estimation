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

int main(int argc, char *argv[]) {
    daf::Timer total_timer, sample_timer;

    std::string data_name = "../../dataset/yeast/data_graph/yeast.graph";
    std::string query_name = "../../dataset/yeast/query_graph/query_dense_4_1.graph";
    uint64_t limit = std::numeric_limits<uint64_t>::max();
    int num_samples = 1000000;
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
            }
        }
    }

    std::cout << "Loading data graph...\n";
    daf::DataGraph data(data_name);
    data.LoadAndProcessGraph();

    std::cout << "Loading query graph...\n";
    daf::QueryGraph query(query_name);

    total_timer.Start();
    query.LoadAndProcessGraph(data);
    total_timer.Stop();

    daf::DAG dag(data, query);

    total_timer.Start();
    dag.BuildDAG();
    total_timer.Stop();

    daf::CandidateSpace cs(data, query, dag);

    total_timer.Start();
    bool cs_constructed = cs.BuildCS();
    total_timer.Stop();

    for (auto u = 0; u < query.GetNumVertices(); ++u) {
        if (cs.GetCandidateSetSize(u) == 0) {
            std::cout << "Total time: " << total_timer.GetTime() << " ms\n";
            return 1;
        }
    }
    sample_timer.Start();
    double est = cs.EstimateEmbeddings(num_samples);
    sample_timer.Stop();

    total_timer.Add(sample_timer);

    std::cout << "Total time: " << total_timer.GetTime() << " ms\n";
    std::cout << "Search time: " << sample_timer.GetTime() << " ms\n";
    std::cout << "#Matches(Approx) : " << std::fixed << std::setprecision(6) << est << std::endl;
    std::cout << "#Tree : " << std::fixed << std::setprecision(6) << cs.total_trees_ << std::endl;
}
