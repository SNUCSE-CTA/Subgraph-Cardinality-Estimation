#ifndef QUERY_GRAPH_H_
#define QUERY_GRAPH_H_

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "include/data_graph.h"
#include "include/graph.h"

namespace CardEst {

class QueryGraph : public Graph {
public:
    explicit QueryGraph(const std::string &filename);

    ~QueryGraph();

    QueryGraph &operator=(const QueryGraph &) = delete;

    QueryGraph(const QueryGraph &) = delete;

    bool LoadAndProcessGraph(const DataGraph &data);
    inline Label GetMaxLabel() const;

    std::vector<std::vector<Vertex>> verticesbyLabel;

    std::vector<std::vector<Vertex>> triangles;

    std::vector<Vertex> GetVerticesByLabel(Size l) const {
        return verticesbyLabel[l];
    }

    bool ProcessLabeledGraph(const DataGraph &data);
    std::string myname;


private:
    Label max_label_;
};


inline Label QueryGraph::GetMaxLabel() const { return max_label_; }
}
#endif  // QUERY_GRAPH_H_
