#ifndef GLOBAL_GLOBAL_H_
#define GLOBAL_GLOBAL_H_

#include <cassert>
#include <cstdint>
#include <limits>
#include <vector>
#include "global/log.h"

namespace daf {
using Size = uint32_t;
using Vertex = uint32_t;
using Label = uint32_t;
using QueryDegree = uint8_t;

constexpr Size INVALID_SZ = std::numeric_limits<Size>::max();
constexpr Vertex INVALID_VTX = std::numeric_limits<Vertex>::max();
constexpr Label INVALID_LB = std::numeric_limits<Label>::max();
}  // namespace daf

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

#endif  // GLOBAL_GLOBAL_H_
