#ifndef GRAPH_H_
#define GRAPH_H_

#include <map>
#include <set>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

#include "global/global.h"
#include <boost/dynamic_bitset.hpp>

namespace daf {
    struct EdgeInfo {
        int index, edge_label, vsum;
        VertexPair vp;
        boost::dynamic_bitset<uint64_t> edge_candidacy_;
    };
    class Graph {
    public:
        explicit Graph(const std::string &filename);

        ~Graph();

        Graph &operator=(const Graph &) = delete;

        Graph(const Graph &) = delete;

        inline Size GetNumLabels() const;

        inline Size GetNumVertices() const;

        inline Size GetNumEdges() const;

        inline Size GetMaxDegree() const;

        inline Label GetLabel(Vertex v) const;

        inline Size GetStartOffset(Vertex v) const;

        inline Size GetEndOffset(Vertex v) const;

        inline Size GetDegree(Vertex v) const;

        inline double GetAvgDegree() const;

        inline Size GetCoreNum(Vertex v) const;

        inline Label GetLabelFrequency(Label l) const;

        inline Vertex GetNeighbor(Size i) const;

        inline bool CheckEdgeExist(Vertex u, Vertex v) const;
        inline Label GetELabel(int edge_idx) const {
            return edge_info_[edge_idx].edge_label;
        }

        tsl::robin_map<std::pair<Vertex, Vertex>, Label> edge_labels_;
        std::vector<tsl::robin_set<Vertex>> edge_exists;
        std::vector<VertexPair> all_edges;
        Size num_label_;


        inline std::vector<int>& GetIncidentEdges(Vertex v, Label l);

        std::vector<EdgeInfo> edge_info_;
        std::vector<std::vector<std::vector<int>>> incident_edges_;
        tsl::robin_map<VertexPair, Size> edge_index_map_;

    protected:
        std::unordered_map<Label, Label> transferred_label_map;
        Label *true_label_;

        void LoadRoughGraph(std::vector<std::vector<Vertex>> *graph);

        void computeCoreNum();

        Size num_vertex_;
        Size num_edge_;

        Size max_degree_;

        Label *label_;
        Size *start_off_;
        Vertex *linear_adj_list_;
        Size *label_frequency_;

        Size *core_num_;

        const std::string &filename_;
        std::ifstream fin_;
    public:
        std::vector<std::vector<Vertex>> adj_list;

        int GetEdgeIndex(Vertex u, Vertex v);

        Vertex opposite(int edge_idx, Vertex from);
    };

    inline Size Graph::GetNumLabels() const { return num_label_; }

    inline Size Graph::GetNumVertices() const { return num_vertex_; }

    inline Size Graph::GetNumEdges() const { return num_edge_; }

    inline Size Graph::GetMaxDegree() const { return max_degree_; }

    inline Label Graph::GetLabel(Vertex v) const { return label_[v]; }

    inline Size Graph::GetStartOffset(Vertex v) const { return start_off_[v]; }

    inline Size Graph::GetEndOffset(Vertex v) const { return start_off_[v + 1]; }

    inline Size Graph::GetDegree(Vertex v) const {
        return start_off_[v + 1] - start_off_[v];
    }

    inline Size Graph::GetCoreNum(Vertex v) const { return core_num_[v]; }

    inline Label Graph::GetLabelFrequency(Label l) const {
        return label_frequency_[l];
    }

    inline Vertex Graph::GetNeighbor(Size i) const { return linear_adj_list_[i]; }


    inline bool Graph::CheckEdgeExist(Vertex u, Vertex v) const {
        return edge_exists[u].find(v) != edge_exists[u].end();
    }

    inline int Graph::GetEdgeIndex(Vertex u, Vertex v) {
        if (!CheckEdgeExist(u, v)) return -1;
        return edge_index_map_[{u, v}];
    }


    inline std::vector<int>& Graph::GetIncidentEdges(Vertex v, Label l) {
        return incident_edges_[v][l];
    }

    inline Vertex Graph::opposite(int edge_idx, Vertex from) {
        return edge_info_[edge_idx].vsum - from;
    }

    inline double Graph::GetAvgDegree() const {
        return GetNumEdges() * 1.0 / GetNumVertices();
    }

}

#endif  // GRAPH_H_
