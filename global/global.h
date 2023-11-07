#ifndef GLOBAL_GLOBAL_H_
#define GLOBAL_GLOBAL_H_

#include <stack>
#include <iostream>
#include <queue>
#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include <cstring>
#include <random>
#include <istream>
#include <fstream>
#include <boost/math/distributions.hpp>
#include <boost/functional/hash.hpp>

using Size = int32_t;
using Vertex = int32_t;
using Label = int32_t;
using QueryDegree = int8_t;
using VertexPair = std::pair<Vertex, Vertex>;

static const int MAX_QUERY_VERTEX = 50, MAX_QUERY_EDGE = 250;
static std::random_device rd;
static std::mt19937 gen(rd());
constexpr Label INVALID_LB = std::numeric_limits<Label>::max();

namespace std {
    template<>
    struct hash<VertexPair> {
        auto operator()(const VertexPair &x) const -> size_t {
            std::size_t seed = 17;
            boost::hash_combine(seed, x.first);
            boost::hash_combine(seed, x.second);
            return seed;
        }
    };
}  

struct UnionFind {
    std::vector<uint32_t> par, sz;

    UnionFind(uint32_t n = 0) {
        par.resize(n);
        sz.resize(n);
        for (int i = 0; i < n; i++) {
            par[i] = i, sz[i] = 1;
        }
    }

    uint32_t find(uint32_t x) {
        return x == par[x] ? x : (par[x] = find(par[x]));
    }

    bool unite(uint32_t x, uint32_t y) {
        uint32_t u = find(x), v = find(y);
        if (u == v) return false;
        if (sz[u] > sz[v]) std::swap(u, v);
        sz[v] += sz[u];
        sz[u] = 0;
        par[u] = par[v];
        return true;
    }
};

struct BipartiteMaximumMatching {
    int *left, *right;
    int left_len, right_len;
    bool *used;
    int **adj, *adj_size;
    int **adj_index;
    bool **matchable;


    bool *bfs_visited;
    int *right_order, *inverse_right_order;
    int **lower_graph, *lower_graph_size;
    int **upper_graph, *upper_graph_size;

    int *Q, *S;
    int qright = 0, qleft = 0;
    int stkright = 0;
    int dfsn[MAX_QUERY_VERTEX], scch[MAX_QUERY_VERTEX], scc_idx[MAX_QUERY_VERTEX], ord, found_scc;
    int scc_dfs(int v) {
        S[stkright++] = v;
        int ret=dfsn[v]=++ord;
        for(int i = 0; i < upper_graph_size[v]; i++){
            int n = upper_graph[v][i];
            if(!dfsn[n])	ret=std::min(ret, scc_dfs(n));
            else if(!scch[n])	ret=std::min(ret, dfsn[n]);
        }
        if(ret==dfsn[v]){
            int u;
            do {
                u=S[--stkright];
                scch[u]=1;
                scc_idx[u] = found_scc;
            } while(u!=v);
            found_scc++;
        }
        return ret;
    }
    bool FindUnmatchableEdges(int required) {
        ord = 0;
        found_scc = 0;
        memset(dfsn, 0, sizeof (dfsn));
        memset(scch, 0, sizeof (scch));
        memset(scc_idx, -1, sizeof (scc_idx));
        int num_matched_ans = solve();
        if (num_matched_ans != required) return false;
        for (int i = 0; i < num_matched_ans; i++) {
            for (int j = 0; j < adj_size[i]; j++) {
                matchable[i][adj[i][j]] = false;
            }
        }
        for (int i = 0; i < num_matched_ans; i++) {
            matchable[i][left[i]] = true;
        }
        int num_rights = 0;
        for (int i = 0; i < num_matched_ans; i++) {
            inverse_right_order[num_rights] = left[i];
            right_order[left[i]] = num_rights++;
        }
        for (int i = 0; i < num_matched_ans; i++) {
            for (int j = 0; j < adj_size[i]; j++) {
                int r = right_order[adj[i][j]];
                if (r == -1) {
                    inverse_right_order[num_rights] = adj[i][j];
                    right_order[adj[i][j]] = num_rights++;
                    r = right_order[adj[i][j]];
                }
                if (r != i) {
                    lower_graph[r][lower_graph_size[r]++] = i;
                }
                if (r != i and r < num_matched_ans) {
                    upper_graph[i][upper_graph_size[i]++] = r;
                }
            }
        }
        for(int i=0;i<num_matched_ans;i++) {
            if (!dfsn[i])
                scc_dfs(i);
        }
        for (int i = 0; i < num_matched_ans; i++) {
            for (int j = 0; j < upper_graph_size[i]; j++) {
                int r = upper_graph[i][j];
                if (scc_idx[i] == scc_idx[r]) {
                    matchable[i][inverse_right_order[r]] = true;
                }
            }
        }

        std::memset(bfs_visited, 0, sizeof(bool) * right_len);
        for (int i = num_matched_ans; i < num_rights; i++) {
            Q[qright++] = i;
            bfs_visited[i] = true;
        }
        while (qright != qleft) {
            int u = Q[qleft++];
            for (int j = 0; j < lower_graph_size[u]; j++) {
                int v = lower_graph[u][j];
                matchable[v][inverse_right_order[u]] = true;
                if (!bfs_visited[v]) {
                    Q[qright++] = v;
                    bfs_visited[v] = true;
                }
            }
        }
        return true;
    }

    void global_initialize(int max_left, int max_right) {
        Q = new int[max_right];
        S = new int[max_right];
        qleft = qright = stkright = 0;
        left = new int[max_left];
        right = new int[max_right];
        left_len = max_left;
        right_len = max_right;
        used = new bool[max_left];
        adj = new int *[max_left];
        adj_index = new int *[max_left];
        matchable = new bool *[max_left];

        lower_graph = new int *[max_right];
        lower_graph_size = new int[max_right];

        upper_graph = new int *[max_left];
        upper_graph_size = new int[max_left];
        right_order = new int[max_right];
        inverse_right_order = new int[max_right];
        for (int i = 0; i < max_left; i++) {
            adj[i] = new int[max_right];
            adj_index[i] = new int[max_right];
            matchable[i] = new bool[max_right];
            upper_graph[i] = new int[max_left];
        }
        for (int i = 0; i < max_right; i++) {
            lower_graph[i] = new int[max_left];
        }
        adj_size = new int[max_left];
        bfs_visited = new bool[max_right];
    }

    ~BipartiteMaximumMatching() {
        delete[] left;
        delete[] right;
        delete[] used;
        for (int i = 0; i < left_len; i++) {
            delete[] adj[i];
            delete[] matchable[i];
            delete[] upper_graph[i];
        }
        delete[] matchable;
        delete[] adj;
        delete[] adj_size;
        delete[] lower_graph;
        delete[] upper_graph;
        delete[] lower_graph_size;
        delete[] upper_graph_size;
        delete[] right_order;
        delete[] inverse_right_order;
        delete[] bfs_visited;
    }

    void reset(bool reset_edges = true) {
        std::memset(left, -1, sizeof(int) * left_len);
        std::memset(right, -1, sizeof(int) * right_len);
        std::memset(used, false, sizeof(bool) * left_len);
        if (reset_edges) {
            std::memset(adj_size, 0, sizeof(int) * left_len);
        }
        std::memset(lower_graph_size, 0, sizeof(int) * right_len);
        std::memset(upper_graph_size, 0, sizeof(int) * left_len);
        std::memset(right_order, -1, sizeof(int) * right_len);
        std::memset(inverse_right_order, -1, sizeof(int) * right_len);
        qleft = qright = stkright = 0;
    }

    void add_edge(int u, int v) {
        adj_index[u][v] = adj_size[u];
        adj[u][adj_size[u]++] = v;
    }

    int remove_edge(int u, int v) {
        if (adj_size[u] > 1) {
            adj_index[u][adj[u][adj_size[u] - 1]] = adj_index[u][v];
            std::swap(adj[u][adj_size[u] - 1], adj[u][adj_index[u][v]]);
        }
        return --adj_size[u];
    }

    void revert(int *tmp_left) {
        for (int i = 0; i < left_len; i++) {
            if (left[i] == -1) continue;
            right[left[i]] = -1;
        }
        std::memcpy(left, tmp_left, sizeof(int) * left_len);
        for (int i = 0; i < left_len; i++) {
            if (left[i] == -1) continue;
            right[left[i]] = i;
        }
    }

    int solve(int ignore = -1) {
        std::memset(left, -1, sizeof(int) * left_len);
        int ans = 0;
        for (Vertex u = 0; u < left_len; u++) {
            if (u == ignore) continue;
            if (left[u] == -1) {
                std::memset(used, false, sizeof(bool) * left_len);
                if (dfs(u)) {
                    ans++;
                }
            }
        }
        return ans;
    }

    bool dfs(Vertex r) {
        if (used[r]) return false;
        used[r] = true;
        for (int i = 0; i < adj_size[r]; i++) {
            Vertex c = adj[r][i];
            Vertex k = right[c];
            if (k == -1 or dfs(k)) {
                left[r] = c;
                right[c] = r;
                return true;
            }
        }
        return false;
    }

    bool single_dfs(Vertex r) {
        if (used[r]) return false;
        used[r] = true;
        for (int i = 0; i < adj_size[r]; i++) {
            Vertex c = adj[r][i];
            Vertex k = right[c];
            if (k == -1 or single_dfs(k)) {
                return true;
            }
        }
        return false;
    }
};

namespace CardEst {
    enum STRUCTURE_FILTER {
        NO_STRUCTURE_FILTER,
        TRIANGLE_SAFETY,
        FOURCYCLE_SAFETY
    };
    enum NEIGHBORHOOD_FILTER {
        NEIGHBOR_SAFETY,
        NEIGHBOR_BIPARTITE_SAFETY,
        EDGE_BIPARTITE_SAFETY
    };
    enum REFINEMENT_ORDER {
        DAG_DP,
        PRIORITY_FIRST
    };
    struct Option {
        STRUCTURE_FILTER structure_filter = FOURCYCLE_SAFETY;
        NEIGHBORHOOD_FILTER neighborhood_filter = EDGE_BIPARTITE_SAFETY;
        REFINEMENT_ORDER refinement_order = PRIORITY_FIRST;
        int sample_size_K = 50000;
        double cutoff = 0.05;
    };
}
#endif 