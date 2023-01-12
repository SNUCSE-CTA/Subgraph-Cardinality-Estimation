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
    class CandidateSpace {
    public:
        CandidateSpace(DataGraph *data);
        CandidateSpace(DataGraph &data, QueryGraph &query, DAG &dag);

        ~CandidateSpace();

        CandidateSpace &operator=(const CandidateSpace &) = delete;

        CandidateSpace(const CandidateSpace &) = delete;

        inline Size GetCandidateSetSize(Vertex u) const;

        inline Vertex GetCandidate(Vertex u, Size v_idx) const;

        tsl::robin_map<VertexPair, std::vector<VertexPair>> cs_edge_list_;
        tsl::robin_map<VertexPair, tsl::robin_map<Vertex, std::vector<Vertex>>> cs_adj_;


        void printCS();
        bool BuildCS(QueryGraph *query, DAG *dag);
        std::vector<std::vector <Vertex>> candidate_set_;
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
        tsl::hopscotch_map<VertexPair, tsl::hopscotch_map<int, std::pair<int, int>>> trigvertex;
        tsl::hopscotch_map<int, std::vector<VertexPair>> reverse_trigvertex;
        tsl::robin_map<std::pair<int, int>, int*> four_cycle_memo;
//        tsl::robin_map<std::pair<int, int>, std::vector<int>> four_cycle_memo;
        tsl::robin_map<std::pair<int, int>, std::vector<online_cycle_information>> four_cycle_memo_old;

        bool **BitsetCS;
        boost::dynamic_bitset<uint64_t> tmpBitset;
//        std::vector<boost::dynamic_bitset<uint64_t>> BitsetCS;

        bool **BitsetEdgeCS;

        QueryDegree *num_visit_cs_;
        Vertex *visited_candidates_;
        Size num_visited_candidates;


        int num_cs_edges_;
        std::vector<int> neighbor_label_frequency;
        std::vector<int> in_neighbor_cs;

        bool FilterByTopDownWithInit();

        void ConstructCS();

        bool InitRootCandidates();

        void ComputeNbrInformation(Vertex u, Size *max_nbr_degree,
                                   uint64_t *label_set);

        bool BipartiteSafety(Vertex cur, Vertex cand);

        bool Filter(bool topdown);

        Vertex GetDAGNextVertex(Vertex cur, Size idx, bool topdown);

        Size GetDAGNextCount(Vertex cur, bool topdown);

        void PrepareNeighborSafety(Vertex cur);

        bool CheckNeighborSafety(Vertex cur, Vertex cand);

        bool EdgeSafety(int query_edge_id, int data_edge_id);

        bool EdgeCandidacy(int query_edge_id, int data_edge_id);

        bool TriangleSafety(int query_edge_id, int data_edge_id);

        bool FourCycleSafety(int query_edge_id, int data_edge_id);
        bool FourCycleSafetyOnline(int query_edge_id, int data_edge_id);

    };

    inline Size CandidateSpace::GetCandidateSetSize(Vertex u) const {
        return candidate_set_[u].size();
    }

    inline Vertex CandidateSpace::GetCandidate(Vertex u, Size v_idx) const {
        return candidate_set_[u][v_idx];
    }
}  // namespace daf

#endif  // CANDIDATE_SPACE_H_
