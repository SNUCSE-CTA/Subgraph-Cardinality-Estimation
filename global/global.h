#ifndef GLOBAL_GLOBAL_H_
#define GLOBAL_GLOBAL_H_

#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include "global/log.h"
using Size = uint32_t;
using Vertex = uint32_t;
using Label = uint32_t;
using QueryDegree = uint8_t;

constexpr Size INVALID_SZ = std::numeric_limits<Size>::max();
constexpr Vertex INVALID_VTX = std::numeric_limits<Vertex>::max();
constexpr Label INVALID_LB = std::numeric_limits<Label>::max();
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
    std::vector<Vertex> left, right;
    std::vector<bool> used;
    std::vector<std::vector<Vertex>> adj;
    BipartiteMaximumMatching(Size l = 0, Size r = 0){
        left.resize(l, -1);
        right.resize(r, -1);
        adj.resize(l);
    }
    void add_edge(size_t u, size_t v) {
        adj[u].push_back(v);
    }
    // Time: O( V(V+E) )
    size_t solve() {
        std::fill(left.begin(), left.end(), -1);
        size_t ans = 0;
        for (Vertex u = 0; u < left.size(); u++) {
            if (left[u] == -1) {
                std::fill(used.begin(), used.end(), false);
                if (dfs(u)) ans++;
            }
        }
        return ans;
    }

    bool dfs(Vertex r) {
        used[r] = true;
        for (Vertex c : adj[r]) {
            Vertex k = right[c];
            if (k == -1 or !used[k] and dfs(k)) {
                left[r] = c;
                right[c] = r;
                return true;
            }
        }
        return false;
    }
};
#endif  // GLOBAL_GLOBAL_H_
