#ifndef CANDIDATE_SPACE_H_
#define CANDIDATE_SPACE_H_

#include <utility>
#include <vector>
#include <random>
#include <map>
#include <unordered_map>
#include <unordered_set>

#include "global/global.h"
#include "include/daf_dag.h"
#include "include/daf_data_graph.h"
#include "include/daf_query_graph.h"

namespace daf {

    class CandidateSpace {
    public:
        CandidateSpace(const DataGraph &data, const QueryGraph &query, DAG &dag);

        ~CandidateSpace();

        CandidateSpace &operator=(const CandidateSpace &) = delete;

        CandidateSpace(const CandidateSpace &) = delete;

        bool BuildCS();

        inline Size GetCandidateSetSize(Vertex u) const;

        inline Vertex GetCandidate(Vertex u, Size v_idx) const;

        inline Size GetCandidateStartOffset(Vertex u, Size u_adj_idx,
                                            Size v_idx) const;

        inline Size GetCandidateEndOffset(Vertex u, Size u_adj_idx, Size v_idx) const;

        inline Size GetCandidateIndex(Size idx) const;
//        gp_hash_table<VertexPair, gp_hash_table<VertexPair, null_type>> cs_edge_list_;

        std::unordered_map<VertexPair, std::set<VertexPair>> cs_edge_list_;

    private:
        const DataGraph &data_;
        const QueryGraph &query_;
        DAG &dag_;
        Size *candidate_set_size_;
        Vertex **candidate_set_;

        Size ***candidate_offsets_;
        Vertex *linear_cs_adj_list_;

        QueryDegree *num_visit_cs_;
        Vertex *visited_candidates_;
        Size *cand_to_cs_idx_;
        Size num_visitied_candidates_;

        Size num_cs_edges_;

        std::vector<int> neighbor_label_frequency;
        std::vector<int> in_neighbor_cs;

        bool FilterByTopDownWithInit();

        bool FilterByBottomUp();

        bool FilterByTopDown();

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

        bool EdgeSafety(Vertex cur, Vertex cand, Vertex nxt, Vertex nxt_cand);
    };

    inline Size CandidateSpace::GetCandidateSetSize(Vertex u) const {
        return candidate_set_size_[u];
    }

    inline Vertex CandidateSpace::GetCandidate(Vertex u, Size v_idx) const {
        return candidate_set_[u][v_idx];
    }

    inline Size CandidateSpace::GetCandidateStartOffset(Vertex u, Size u_adj_idx,
                                                        Size v_idx) const {
        return candidate_offsets_[u][u_adj_idx][v_idx];
    }

    inline Size CandidateSpace::GetCandidateEndOffset(Vertex u, Size u_adj_idx,
                                                      Size v_idx) const {
        return candidate_offsets_[u][u_adj_idx][v_idx + 1];
    }

    inline Size CandidateSpace::GetCandidateIndex(Size idx) const {
        return linear_cs_adj_list_[idx];
    }

}  // namespace daf

#endif  // CANDIDATE_SPACE_H_
