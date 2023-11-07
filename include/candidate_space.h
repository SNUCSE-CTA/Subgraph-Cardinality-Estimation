#ifndef CANDIDATE_SPACE_H_
#define CANDIDATE_SPACE_H_

#include <utility>
#include <vector>
#include <random>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "global/global.h"
#include "include/dag.h"
#include "include/data_graph.h"
#include "include/query_graph.h"

namespace CardEst {

    struct QueryVertex {
        double priority;
        int stage, which;
        bool operator>(const QueryVertex &o) const {
            return priority > o.priority;
        }
        bool operator<(const QueryVertex &o) const {
            return priority < o.priority;
        }
    };

    class CandidateSpace {
    public:
        CandidateSpace(DataGraph *data, Option filter_option);

        ~CandidateSpace();

        CandidateSpace &operator=(const CandidateSpace &) = delete;

        CandidateSpace(const CandidateSpace &) = delete;

        inline Size GetCandidateSetSize(Vertex u) const;

        inline Vertex GetCandidate(Vertex u, Size v_idx) const;

        std::vector<std::vector<std::vector<std::vector<int>>>> cs_edge_;

        Option opt;

        bool BuildCS(QueryGraph *query, OrderedQueryGraph *dag);
        std::vector<std::vector <Vertex>> candidate_set_;
        std::vector<int> neighbor_label_frequency;
        bool *in_neighbor_cs;
    private:
        DataGraph *data_;
        QueryGraph *query_;
        OrderedQueryGraph *dag_;

        bool **BitsetCS;
        bool **BitsetEdgeCS;

        QueryDegree *num_visit_cs_;

        bool InitializeCS();

        void ConstructCS();

        bool InitRootCandidates();

        bool BipartiteSafety(Vertex cur, Vertex cand);

        bool Filter();

        void PrepareNeighborSafety(Vertex cur);

        bool CheckNeighborSafety(Vertex cur, Vertex cand);

        bool EdgeCandidacy(int query_edge_id, int data_edge_id);

        bool TriangleSafety(int query_edge_id, int data_edge_id);

        bool FourCycleSafety(int query_edge_id, int data_edge_id);

        bool StructureFilter(int cur, int cand, int direction);

        bool NeighborhoodFilter(int cur, int cand);

        bool StructureSafety(int query_edge_id, int data_edge_id);

        bool EdgeBipartiteSafety(Vertex cur, Vertex cand);

        Size GetDAGNextCount(Vertex cur, bool topdown);
    };

    inline Size CandidateSpace::GetCandidateSetSize(Vertex u) const {
        return candidate_set_[u].size();
    }

    inline Vertex CandidateSpace::GetCandidate(Vertex u, Size v_idx) const {
        return candidate_set_[u][v_idx];
    }
}  

#endif  
