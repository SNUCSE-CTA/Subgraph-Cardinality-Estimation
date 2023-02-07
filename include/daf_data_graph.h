#ifndef DATA_GRAPH_H_
#define DATA_GRAPH_H_

#include <algorithm>
#include <climits>
#include <string>
#include <set>
#include <unordered_map>
#include <utility>
#include <vector>
#include "global/global.h"
#include "include/daf_graph.h"

struct triginfo {
    int vertex;
    int fst_edge_idx, snd_edge_idx;
};

namespace std {
    template <>
    struct hash<triginfo> {
        auto operator()(const triginfo &x) const -> size_t {
            return x.vertex;
        }
    };
}

namespace daf {
class DataGraph : public Graph {
public:
    explicit DataGraph(const std::string &filename);

    ~DataGraph();

    DataGraph &operator=(const DataGraph &) = delete;

    DataGraph(const DataGraph &) = delete;

    void LoadAndProcessGraph();
    void ProceesLabeledGraph();

    using Graph::GetEndOffset;
    using Graph::GetStartOffset;

    inline Label GetTransferredLabel(Label l) const;

    inline Size GetStartOffsetByLabel(Label l) const;

    inline Size GetEndOffsetByLabel(Label l) const;

    inline Size GetVertexBySortedLabelOffset(Size i) const;

    inline Size GetInitCandSize(Label l, Size d) const;

    bool is_sparse();
    int max_num_trigs;

private:
    Label *transferred_label_;
    std::pair<Size, Size> *adj_offs_by_label_;

    Size *offs_by_label_;
    Vertex *vertices_sorted_;

    uint64_t *linear_nbr_bitset_;
    Size *max_nbr_degree_;

    Size nbr_bitset_size_;
    Size max_label_frequency_;
};


inline Label DataGraph::GetTransferredLabel(Label l) const {
    return transferred_label_[l];
}

inline Size DataGraph::GetStartOffsetByLabel(Label l) const {
    return offs_by_label_[l];
}

inline Size DataGraph::GetEndOffsetByLabel(Label l) const {
    return offs_by_label_[l + 1];
}

inline Size DataGraph::GetVertexBySortedLabelOffset(Size i) const {
    return vertices_sorted_[i];
}

inline Size DataGraph::GetInitCandSize(Label l, Size d) const {
    Size s = GetStartOffsetByLabel(l);
    Size e = GetEndOffsetByLabel(l);
    auto pos = std::lower_bound(
            vertices_sorted_ + s, vertices_sorted_ + e, d,
            [this](Vertex v, Size d) -> bool { return GetDegree(v) >= d; });
    return pos - (vertices_sorted_ + s);
}

inline bool DataGraph::is_sparse() {
    return GetNumEdges() / GetNumVertices() < 10;
}

}

#endif  // DATA_GRAPH_H_
