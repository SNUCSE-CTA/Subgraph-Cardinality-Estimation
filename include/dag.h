#ifndef DAG_H_
#define DAG_H_

#include <algorithm>

#include "global/global.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

namespace daf {
    class DAG {
    public:
        DAG(const DataGraph &data, const QueryGraph &query);

        ~DAG();

        DAG &operator=(const DAG &) = delete;

        DAG(const DAG &) = delete;

        void BuildDAG(Vertex root_vertex=-1);

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

    private:
        const DataGraph &data_;
        const QueryGraph &query_;
        std::vector<Size> bfs_level_, init_cand_size_, kth_child_;
        std::vector<Vertex> tree_sequence_, bfs_sequence_, tree_parent_;
        std::vector<std::vector<Vertex>> children_, parents_, tree_children_, tree_neighbors_;
        Vertex SelectRootVertex();
    };

    inline Vertex DAG::GetRoot() const { return bfs_sequence_[0]; }

    inline Size DAG::GetNumChildren(Vertex v) const { return children_[v].size(); }
    inline Size DAG::GetNumTreeChildren(Vertex v) const { return tree_children_[v].size(); }

    inline Size DAG::GetNumParents(Vertex v) const { return parents_[v].size(); }

    inline Size DAG::GetTreeChild(Vertex v, Size i) const { return tree_children_[v][i]; }
    inline Size DAG::GetChild(Vertex v, Size i) const { return children_[v][i]; }

    inline Size DAG::GetTreeParent(Vertex v, Size i) const { return tree_parent_[v];}
    inline Size DAG::GetParent(Vertex v, Size i) const { return parents_[v][i]; }

    inline Size DAG::GetInitCandSize(Vertex v) const { return init_cand_size_[v]; }

    inline Vertex DAG::GetVertexOrderedByBFS(Size i) const {
        return bfs_sequence_[i];
    }

    inline Vertex DAG::GetVertexOrderedByTree(Size i) const {
        return tree_sequence_[i];
    }

}  // namespace daf

#endif  // DAG_H_
