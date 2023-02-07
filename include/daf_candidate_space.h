#ifndef CANDIDATE_SPACE_H_
#define CANDIDATE_SPACE_H_

#include <utility>
#include <vector>
#include <random>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <boost/dynamic_bitset.hpp>

#include "global/global.h"
#include "include/daf_dag.h"
#include "include/daf_data_graph.h"
#include "include/daf_query_graph.h"

namespace daf {
    enum STRUCTURE_FILTER {
        NO_STRUCTURE_FILTER,
        TRIANGLE_SAFETY,
        TRIANGLE_BIPARTITE_SAFETY,
        FOURCYCLE_SAFETY
    };
    enum EGONET_FILTER {
        NEIGHBOR_SAFETY,
        NEIGHBOR_BIPARTITE_SAFETY,
        EDGE_BIPARTITE_SAFETY
    };

    struct FilterOption {
        STRUCTURE_FILTER structure_filter = FOURCYCLE_SAFETY;
        EGONET_FILTER egonet_filter = EDGE_BIPARTITE_SAFETY;
        void print() {
            fprintf(stderr, "Filtering Level : Struct[%d] Egonet[%d]\n", structure_filter, egonet_filter);
            fprintf(stdout, "Filtering Level : Struct[%d] Egonet[%d]\n", structure_filter, egonet_filter);
        }
    };


    struct variable {
        double priority;
        int stage, which;
        bool operator>(const variable &o) const {
            return priority > o.priority;
        }
        bool operator<(const variable &o) const {
            return priority < o.priority;
        }
    };

    class CandidateSpace {
    public:
        CandidateSpace(DataGraph *data, FilterOption filter_option);

        ~CandidateSpace();

        CandidateSpace &operator=(const CandidateSpace &) = delete;

        CandidateSpace(const CandidateSpace &) = delete;

        inline Size GetCandidateSetSize(Vertex u) const;

        inline Vertex GetCandidate(Vertex u, Size v_idx) const;

        std::vector<std::vector<std::vector<std::vector<int>>>> cs_edge_;

        FilterOption opt;
        void printCS();
        bool BuildCS(QueryGraph *query, DAG *dag);
        std::vector<std::vector <Vertex>> candidate_set_;
        std::vector<int> neighbor_label_frequency;
        bool *in_neighbor_cs;
    private:
        DataGraph *data_;
        QueryGraph *query_;
        DAG *dag_;


        struct online_cycle_information {
            online_cycle_information() {
                d_opp_edge_idx = -1;
                third_inc_idx = fourth_inc_idx = 0;
            }
            int q_opp_edge_idx;
            int third_inc_idx, fourth_inc_idx;
            int third, fourth;
            int q_third_edge_idx, q_fourth_edge_idx; // tex = 2->3, fex = 1->4
            int d_third_edge_idx, d_fourth_edge_idx; // tex = 2->3, fex = 1->4
            int d_opp_edge_idx; // 3->4
        };
        std::unordered_map<std::pair<int, int>, std::vector<online_cycle_information>> four_cycle_memo_old;

        bool **BitsetCS;
        bool **BitsetEdgeCS;

        QueryDegree *num_visit_cs_;
        Vertex *visited_candidates_;
        Size num_visited_candidates;

        bool FilterByTopDownWithInit();

        void ConstructCS();

        bool InitRootCandidates();

        bool BipartiteSafety(Vertex cur, Vertex cand);

        bool Filter(bool topdown);

        Vertex GetDAGNextVertex(Vertex cur, Size idx, bool topdown);

        Size GetDAGNextCount(Vertex cur, bool topdown);

        void PrepareNeighborSafety(Vertex cur);

        bool CheckNeighborSafety(Vertex cur, Vertex cand);

        bool EdgeCandidacy(int query_edge_id, int data_edge_id);

        bool TriangleSafety(int query_edge_id, int data_edge_id);
        bool TriangleBipartiteSafety(int query_edge_id, int data_edge_id);

        bool FourCycleSafety(int query_edge_id, int data_edge_id);
        bool FourCycleSafetyOnline(int query_edge_id, int data_edge_id);

        bool BipartiteEdgeSafety(Vertex cur, Vertex cand, Vertex nxt, Vertex nxt_cand);

        bool StructureFilter(int cur, int cand);

        bool EgonetFilter(int cur, int cand);

        bool StructureSafety(int query_edge_id, int data_edge_id);

        bool EdgeBipartiteSafety(Vertex cur, Vertex cand);

        void ComputeNbrInformation(Vertex u, Size *max_nbr_degree, uint64_t *nbr_label_bitset);
    };

    inline Size CandidateSpace::GetCandidateSetSize(Vertex u) const {
        return candidate_set_[u].size();
    }

    inline Vertex CandidateSpace::GetCandidate(Vertex u, Size v_idx) const {
        return candidate_set_[u][v_idx];
    }
}  // namespace daf

#endif  // CANDIDATE_SPACE_H_
