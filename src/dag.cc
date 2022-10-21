#include "include/dag.h"
#include <queue>
namespace daf {
    DAG::DAG(const DataGraph &data, const QueryGraph &query)
            : data_(data), query_(query) {
        // initialize
        Size qV = query_.GetNumVertices();
        bfs_sequence_.resize(qV);
        tree_sequence_.resize(qV);
        bfs_level_.resize(qV);
        children_.resize(qV);
        tree_children_.resize(qV);
        parents_.resize(qV);
        tree_parent_.resize(qV);
        tree_neighbors_.resize(qV);
        init_cand_size_.resize(qV);
        std::fill(tree_parent_.begin(), tree_parent_.end(), -1);
    }

    DAG::~DAG() = default;
    void DAG::BuildTree() {
        bool *visit = new bool[query_.GetNumVertices()];
        std::fill(visit, visit + query_.GetNumVertices(), false);
        std::vector<Size> v_id(query_.GetNumVertices());
        tree_sequence_[0] = bfs_sequence_[0];
        std::queue <Vertex> q;
        q.push(tree_sequence_[0]);
        visit[tree_sequence_[0]] = true;
        v_id[tree_sequence_[0]] = 0;
        Size id = 1;
        while (!q.empty()) {
            Vertex v = q.front();
            q.pop();
            for (Size c : tree_neighbors_[v]) {
                if (!visit[c]) {
                    q.push(c);
                    visit[c] = true;
                    v_id[c] = id;
                    tree_sequence_[id++] = c;
                }
            }
        }
        for (Size i = 0; i < query_.GetNumVertices(); i++) {
            for (Vertex c : tree_neighbors_[i]) {
                if (c > i) continue;
                if (v_id[i] < v_id[c]) {
                    tree_children_[i].push_back(c);
                    tree_parent_[c] = i;
                }
                else {
                    tree_children_[c].push_back(i);
                    tree_parent_[i] = c;
                }
            }
        }
//        for (Size i = 0; i < query_.GetNumVertices(); i++) {
//            fprintf(stderr, "%d : ",i);
//            for (Vertex cid = 0; cid < GetNumTreeChildren(i); cid++) {
//                fprintf(stderr, "%d ", tree_children_[i][cid]);
//            }
//            fprintf(stderr, "\n");
//        }
        delete[] visit;
    }
    void DAG::BuildDAG(Vertex root_vertex) {
        bool *visit = new bool[query_.GetNumVertices()];
        bool *popped = new bool[query_.GetNumVertices()];

        std::fill(visit, visit + query_.GetNumVertices(), false);
        std::fill(popped, popped + query_.GetNumVertices(), false);

        // used as queue (topologically ordered)
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

                for (Size i = query_.GetStartOffset(parent);
                     i < query_.GetEndOffset(parent); ++i) {
                    Vertex child = query_.GetNeighbor(i);

                    if (!popped[child]) {
                        // build edge from parent to child
                        bfs_level_[child] = bfs_level_[parent] + 1;
                        children_[parent].push_back(child);
                        parents_[child].push_back(parent);

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

    Vertex DAG::SelectRootVertex() {
        Vertex root = 0;
        double min_rank = std::numeric_limits<double>::max();
        for (Vertex v = 0; v < query_.GetNumVertices(); ++v) {
            Label l = query_.GetLabel(v);
            Size d = query_.GetDegree(v);
            init_cand_size_[v] = data_.GetInitCandSize(l, d);

            if (query_.GetCoreNum(v) < 2 && !query_.IsTree()) continue;

            double rank =
                    static_cast<double>(init_cand_size_[v]) / static_cast<double>(d);

            if (rank < min_rank) {
                root = v;
                min_rank = rank;
            }
        }

        return root;
    }
}  // namespace daf
