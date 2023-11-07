#ifndef DAG_H_
#define DAG_H_

#include <algorithm>

#include "global/global.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

namespace CardEst {
    class OrderedQueryGraph {
    public:
        OrderedQueryGraph(const DataGraph &data, const QueryGraph &query);

        ~OrderedQueryGraph();

        OrderedQueryGraph &operator=(const OrderedQueryGraph &) = delete;

        OrderedQueryGraph(const OrderedQueryGraph &) = delete;

        void BuildOrderedQueryGraph(Vertex root_vertex=-1);

        inline Vertex GetRoot() const;

        inline Size GetNumChildren(Vertex v) const;

        inline Size GetNumTreeChildren(Vertex v) const;

        inline Size GetNumParents(Vertex v) const;

        inline Size GetTreeChild(Vertex v, Size i) const;

        inline Size GetChild(Vertex v, Size i) const;

        inline Size GetTreeParent(Vertex v, Size i) const;

        inline Size GetParent(Vertex v, Size i) const;

        inline Size GetInitCandSize(Vertex v) const;

        inline Vertex GetVertexOrderedByBFS(Size i) const;
        inline Vertex GetVertexOrderedByTree(Size i) const;

        inline Size GetBFSLevel(Vertex v) const {
            return bfs_level_[v];
        }

        void AddTreeEdge(Vertex u, Vertex v) {
            tree_neighbors_[u].push_back(v);
            tree_neighbors_[v].push_back(u);
        }
        void BuildTree();

        inline Size GetChildIndex(Vertex v) const {
            return kth_child_[v];
        }

        std::vector<std::vector<Vertex>> tree_neighbors_;
        const DataGraph &data_;
        const QueryGraph &query_;
        std::vector<Size> bfs_level_, init_cand_size_, kth_child_;
        std::vector<Vertex> tree_sequence_, bfs_sequence_, tree_parent_;
        std::vector<std::vector<Vertex>> children_, parents_, tree_children_, dag_parent_edges_, dag_child_edges_;

        Vertex SelectRootVertex();
    };

    inline Vertex OrderedQueryGraph::GetRoot() const { return bfs_sequence_[0]; }

    inline Size OrderedQueryGraph::GetNumChildren(Vertex v) const { return children_[v].size(); }
    inline Size OrderedQueryGraph::GetNumTreeChildren(Vertex v) const { return tree_children_[v].size(); }

    inline Size OrderedQueryGraph::GetNumParents(Vertex v) const { return parents_[v].size(); }

    inline Size OrderedQueryGraph::GetTreeChild(Vertex v, Size i) const { return tree_children_[v][i]; }
    inline Size OrderedQueryGraph::GetChild(Vertex v, Size i) const { return children_[v][i]; }

    inline Size OrderedQueryGraph::GetTreeParent(Vertex v, Size i) const { return tree_parent_[v];}
    inline Size OrderedQueryGraph::GetParent(Vertex v, Size i) const { return parents_[v][i]; }

    inline Size OrderedQueryGraph::GetInitCandSize(Vertex v) const { return init_cand_size_[v]; }

    inline Vertex OrderedQueryGraph::GetVertexOrderedByBFS(Size i) const {
        return bfs_sequence_[i];
    }

    inline Vertex OrderedQueryGraph::GetVertexOrderedByTree(Size i) const {
        return tree_sequence_[i];
    }

}

#endif