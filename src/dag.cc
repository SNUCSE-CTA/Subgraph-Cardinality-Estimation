#include "include/dag.h"
#include <queue>
namespace CardEst {
    OrderedQueryGraph::OrderedQueryGraph(const DataGraph &data, const QueryGraph &query)
            : data_(data), query_(query) {
        
        Size qV = query_.GetNumVertices();
        bfs_sequence_.resize(qV);
        tree_sequence_.resize(qV);
        bfs_level_.resize(qV);
        children_.resize(qV);
        tree_children_.resize(qV);
        parents_.resize(qV);
        dag_parent_edges_.resize(qV);
        dag_child_edges_.resize(qV);
        tree_parent_.resize(qV);
        tree_neighbors_.resize(qV);
        init_cand_size_.resize(qV);
        kth_child_.resize(qV);
        std::fill(tree_parent_.begin(), tree_parent_.end(), -1);
    }

    OrderedQueryGraph::~OrderedQueryGraph() = default;
    void OrderedQueryGraph::BuildTree() {

        bool *visit = new bool[query_.GetNumVertices()];
        std::fill(visit, visit + query_.GetNumVertices(), false);
        std::vector<Size> v_id(query_.GetNumVertices());
        tree_children_.clear();
        tree_children_.resize(query_.GetNumVertices());

        std::fill(tree_sequence_.begin(), tree_sequence_.end(), 0);
        std::fill(tree_parent_.begin(), tree_parent_.end(), -1);
        kth_child_.resize(query_.GetNumVertices(), 0);

        tree_sequence_[0] = bfs_sequence_[0];
        std::queue <Vertex> q;
        q.push(tree_sequence_[0]);
        visit[tree_sequence_[0]] = true;
        v_id[tree_sequence_[0]] = 0;
        Size id = 1;
        while (!q.empty()) {
            Vertex v = q.front();
            q.pop();

            for (Vertex c : tree_neighbors_[v]) {
                if (!visit[c]) {
                    q.push(c);
                    visit[c] = true;
                    kth_child_[c] = tree_children_[v].size();
                    tree_children_[v].push_back(c);
                    tree_parent_[c] = v;
                    v_id[c] = id;
                    tree_sequence_[id++] = c;
                }
            }
        }

        delete[] visit;
    }
    void OrderedQueryGraph::BuildOrderedQueryGraph(Vertex root_vertex) {
        bool *visit = new bool[query_.GetNumVertices()];
        bool *popped = new bool[query_.GetNumVertices()];

        std::fill(visit, visit + query_.GetNumVertices(), false);
        std::fill(popped, popped + query_.GetNumVertices(), false);

        if (root_vertex == -1) root_vertex = SelectRootVertex();
        bfs_sequence_[0] = root_vertex;
        visit[bfs_sequence_[0]] = true;
        bfs_level_[0] = 1;

        Size begin = 0;
        Size end = 1;

        while (begin < end) {
            std::sort(bfs_sequence_.begin() + begin, bfs_sequence_.begin() + end,
                      [this](Vertex v1, Vertex v2) -> bool {
                          Size d1 = query_.GetDegree(v1);
                          Size d2 = query_.GetDegree(v2);
                          return d1 > d2;
                      });
            std::stable_sort(bfs_sequence_.begin() + begin, bfs_sequence_.begin() + end,
                             [this](Vertex v1, Vertex v2) -> bool {
                                 Label l1 = query_.GetLabel(v1);
                                 Label l2 = query_.GetLabel(v2);
                                 Size lf1 = data_.GetLabelFrequency(l1);
                                 Size lf2 = data_.GetLabelFrequency(l2);
                                 if (lf1 != lf2) return lf1 < lf2;
                                 return l1 < l2;
                             });
            Size cur_level_end = end;
            while (begin < cur_level_end) {
                Vertex parent = bfs_sequence_[begin];
                begin += 1;
                popped[parent] = true;
                for (Vertex child : query_.adj_list[parent]) {
                    if (!popped[child]) {
                        bfs_level_[child] = bfs_level_[parent] + 1;
                        children_[parent].push_back(child);
                        dag_child_edges_[parent].push_back(query_.GetEdgeIndex(parent, child));
                        parents_[child].push_back(parent);
                        dag_parent_edges_[child].push_back(query_.GetEdgeIndex(child, parent));
                        if (!visit[child]) {
                            visit[child] = true;
                            bfs_sequence_[end] = child;
                            end += 1;
                        }
                    }
                }
            }
        }
        delete[] visit;
        delete[] popped;
    }

    Vertex OrderedQueryGraph::SelectRootVertex() {
        Vertex root = 0;
        double min_rank = std::numeric_limits<double>::max();
        for (Vertex v = 0; v < query_.GetNumVertices(); ++v) {
            Label l = query_.GetLabel(v);
            Size d = query_.GetDegree(v);
            init_cand_size_[v] = data_.GetInitCandSize(l, d);

            double rank =
                    static_cast<double>(init_cand_size_[v]) / static_cast<double>(d);
            if (rank < min_rank) {
                root = v;
                min_rank = rank;
            }
        }
        return root;
    }
}  
