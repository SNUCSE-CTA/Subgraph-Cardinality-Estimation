#include <algorithm>
#include <iostream>
#include <string>
#include <deque>

#include "global/timer.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"
#include "include/FaSTest.h"

using namespace CardEst;

std::string dataset = "yeast", ans_file_name, data_root = "../dataset/";
std::string data_name, query_name = "../dataset/yeast/query_graph/query_sparse_16_110.graph";
std::deque<std::string> query_names;

bool solve_all = false;
int q_cnt = 0;
std::unordered_map<std::string, long long> true_cnt;
FaSTest *Estimator;
Option OPTION;

void estimate(DataGraph &data, QueryGraph &query) {
    GlobalTimer.Start();
    Timer total_timer;
    query.LoadAndProcessGraph(data);

    total_timer.Start();
    CardEst::OrderedQueryGraph g(data, query); g.BuildOrderedQueryGraph();
    Estimator->PrepareQuery(&query, &g);
    double est = Estimator->EstimateEmbeddings();
    total_timer.Stop();
    q_cnt++;
    long long tcnt = true_cnt.find(query_name) == true_cnt.end() ? -1 : true_cnt[query_name];
    fprintf(stdout, "[%d / %d] Query : %s\n",q_cnt, (int)query_names.size(),query_name.c_str());
    fprintf(stdout, "  Est:  %20lld\n  True: %20lld\n  Time: %20.02lf ms\n\n",
         llrint(est), tcnt, total_timer.GetTime());
}

void PrepareAnswers() {
    ans_file_name = data_root+dataset+"/"+dataset+"_ans.txt";
    data_name = data_root+dataset+"/data_graph/"+dataset+".graph";
    if (solve_all) {
        query_names.clear();
    }
    std::ifstream ans_in(ans_file_name);
    while (!ans_in.eof()) {
        std::string name, c;
        ans_in >> name >> c;
        if (name.empty() || c.empty()) continue;
        name = data_root+dataset+"/query_graph/"+name;
        true_cnt[name] = atoll(c.c_str());
        if (solve_all) {
            query_names.push_back(name);
        }
    }
}


void read_filter_option(const std::string& opt, const std::string &filter) {
    if (opt.substr(2) == "NEIGHBORHOOD") {
        if (filter == "NS")
            OPTION.neighborhood_filter = CardEst::NEIGHBOR_SAFETY;
        else if (filter == "NB")
            OPTION.neighborhood_filter = CardEst::NEIGHBOR_BIPARTITE_SAFETY;
        else
            OPTION.neighborhood_filter = CardEst::EDGE_BIPARTITE_SAFETY;
    }
    else if (opt.substr(2) == "STRUCTURE") {
        if (filter == "X")
            OPTION.structure_filter = CardEst::NO_STRUCTURE_FILTER;
        else if (filter == "3")
            OPTION.structure_filter = CardEst::TRIANGLE_SAFETY;
        else if (filter == "4")
            OPTION.structure_filter = CardEst::FOURCYCLE_SAFETY;
    }
    else if (opt.substr(2) == "CUTOFF") {
        OPTION.cutoff = atof(opt.c_str());
    }
}


int main(int argc, char *argv[]) {
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'D':
                    dataset = argv[i + 1];
                    break;
                case 'd':
                    data_name = argv[i + 1];
                    break;
                case 'q':
                    query_name = argv[i + 1];
                    query_names = {query_name};
                    break;
                case 's':
                    OPTION.sample_size_K = std::atoi(argv[i + 1]);
                    break;
                case 'A':
                    solve_all = true;
                    break;
                case '-':
                    read_filter_option(std::string(argv[i]), std::string(argv[i+1]));
                    break;
            }
        }
    }
    PrepareAnswers();

    std::cout << "Loading data graph..." << std::endl;
    DataGraph data(data_name);
    data.LoadAndProcessGraph();

    Estimator = new FaSTest(&data, OPTION);

    for (std::string &qname : query_names) {
        query_name = qname;
        QueryGraph query(qname);
        estimate(data, query);
    }
}
