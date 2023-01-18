#ifndef GLOBAL_GLOBAL_H_
#define GLOBAL_GLOBAL_H_

#include <queue>
#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include <cstring>
#include <random>
#include <boost/functional/hash.hpp>
#include "include/tsl/robin_set.h"
#include "include/tsl/robin_map.h"
#include "include/tsl/hopscotch_set.h"
#include "include/tsl/hopscotch_map.h"
#include "global/log.h"

static int functionCallCounter;
using Size = int32_t;
using Vertex = int32_t;
using Label = int32_t;
using QueryDegree = int8_t;
using VertexPair = std::pair<Vertex, Vertex>;
static uint64_t counters[1000];

static int MAX_QUERY_VERTEX = 50, MAX_QUERY_EDGE = 250;
namespace daf {
    static std::random_device rd;
    static std::mt19937 gen(rd());
}


constexpr Size INVALID_SZ = std::numeric_limits<Size>::max();
constexpr Vertex INVALID_VTX = std::numeric_limits<Vertex>::max();
constexpr Label INVALID_LB = std::numeric_limits<Label>::max();

namespace std {
    template <>
    struct hash<VertexPair> {
        auto operator()(const VertexPair &x) const -> size_t {
            std::size_t seed = 17;
            boost::hash_combine(seed, x.first);
            boost::hash_combine(seed, x.second);
            return seed;
        }
    };
}  // namespace std

struct UnionFind {
    std::vector<uint32_t> par, sz;
    UnionFind(uint32_t n = 0){
        par.resize(n);
        sz.resize(n);
        for (int i = 0; i < n; i++)
            par[i] = i, sz[i] = 1;
    }
    uint32_t find(uint32_t x) {
        return x == par[x] ? x : (par[x] = find(par[x]));
    }
    bool unite(uint32_t x, uint32_t y) {
        uint32_t u = find(x), v = find(y);
        if(u == v) return false;
        if (sz[u]>sz[v]) std::swap(u, v);
        sz[v]+=sz[u];
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

    void global_initialize(int max_left, int max_right) {
//        fprintf(stderr, "BPSolver init to %d by %d\n",max_left,max_right);
        left = new int[max_left];
        right = new int[max_right];
        left_len = max_left;
        right_len = max_right;
        used = new bool[max_left];
        adj = new int*[max_left];
        for (int i = 0; i < max_left; i++) {
            adj[i] = new int[max_right];
        }
        adj_size = new int[max_left];
    }
    ~BipartiteMaximumMatching() {
        delete[] left;
        delete[] right;
        delete[] used;
        for (int i = 0; i < left_len; i++) {
            delete[] adj[i];
        }
        delete[] adj;
        delete[] adj_size;
    }
    void reset(bool reset_edges = true) const {
        std::memset(left, -1, sizeof(int) * left_len);
        std::memset(right, -1, sizeof(int) * right_len);
        std::memset(used, false, sizeof(bool) * left_len);
        if (reset_edges) {
            std::memset(adj_size, 0, sizeof(int) * left_len);
        }
    }
    void add_edge(int u, int v) {
        adj[u][adj_size[u]++] = v;
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
    // Time: O( V(V+E) )
    int solve(int ignore = -1) {
        std::memset(left, -1, sizeof(int) * left_len);
        int ans = 0;
        for (Vertex u = 0; u < left_len; u++) {
            if (u == ignore) continue;
            if (left[u] == -1) {
                std::memset(used, false, sizeof(bool) * left_len);
                if (dfs(u))
                    ans++;
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
//        printf("Single_DFS %d\n",r);
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

    void print() {
        for (int i = 0; i < left_len; i++) {
            printf("%d[%d]\t",left[i],used[i]);
        }
        printf("\n");
    }
};
//
//struct BipartiteMaximumMatching {
//    int n_left, n_right, flow = 0;
//    std::vector<std::vector<int>> g;
//    std::vector<int> match_from_left, match_from_right;
//
//    BipartiteMaximumMatching(int _n_left, int _n_right)
//            : n_left(_n_left),
//              n_right(_n_right),
//              g(_n_left),
//              match_from_left(_n_left, -1),
//              match_from_right(_n_right, -1),
//              dist(_n_left) {}
//
//    void add_edge(int u, int v) { g[u].push_back(v); }
//
//    std::vector<int> dist;
//
//    void bfs() {
//        std::queue<int> q;
//        for (int u = 0; u < n_left; ++u) {
//            if (!~match_from_left[u])
//                q.push(u), dist[u] = 0;
//            else
//                dist[u] = -1;
//        }
//        while (!q.empty()) {
//            int u = q.front();
//            q.pop();
//            for (auto v : g[u])
//                if (~match_from_right[v] && !~dist[match_from_right[v]]) {
//                    dist[match_from_right[v]] = dist[u] + 1;
//                    q.push(match_from_right[v]);
//                }
//        }
//    }
//
//    bool dfs(int u) {
//        for (auto v : g[u])
//            if (!~match_from_right[v]) {
//                match_from_left[u] = v, match_from_right[v] = u;
//                return true;
//            }
//        for (auto v : g[u])
//            if (dist[match_from_right[v]] == dist[u] + 1 &&
//                dfs(match_from_right[v])) {
//                match_from_left[u] = v, match_from_right[v] = u;
//                return true;
//            }
//        return false;
//    }
//
//    int solve() {
//        while (true) {
//            bfs();
//            int augment = 0;
//            for (int u = 0; u < n_left; ++u)
//                if (!~match_from_left[u]) augment += dfs(u);
//            if (!augment) break;
//            flow += augment;
//        }
//        return flow;
//    }
//
//    std::vector<std::pair<int, int>> get_edges() {
//        std::vector<std::pair<int, int>> ans;
//        for (int u = 0; u < n_left; ++u)
//            if (match_from_left[u] != -1)
//                ans.emplace_back(u, match_from_left[u]);
//        return ans;
//    }
//};
#endif  // GLOBAL_GLOBAL_H_
