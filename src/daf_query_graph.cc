#include "include/daf_query_graph.h"
using namespace daf;
QueryGraph::QueryGraph(const std::string &filename) : Graph(filename) {
    myname = filename;
}

QueryGraph::~QueryGraph() {
    delete[] NEC_map_;
    delete[] NEC_elems_;

    if (NEC_start_offs_) delete[] NEC_start_offs_;
    delete[] NEC_size_;
}

bool QueryGraph::ProcessLabeledGraph(const DataGraph &data) {
    std::fill(label_frequency_, label_frequency_ + data.GetNumLabels(), 0);

    Size cur_idx = 0;

    // transfer label & construct adj list and label frequency
    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        Label l = data.GetTransferredLabel(label_[v]);
        if (l == INVALID_LB) return false;
        if (label_frequency_[l] == 0) num_label_ += 1;
        label_[v] = l;
        verticesbyLabel[l].push_back(v);
        max_label_ = std::max(max_label_, l);
        label_frequency_[l] += 1;
        if (adj_list[v].size() > max_degree_) max_degree_ = adj_list[v].size();

        start_off_[v] = cur_idx;
        core_num_[v] = adj_list[v].size();

        std::copy(adj_list[v].begin(), adj_list[v].end(),
                  linear_adj_list_ + cur_idx);

        cur_idx += adj_list[v].size();
    }
    start_off_[GetNumVertices()] = num_edge_ * 2;

    // preprocess for query graph
    computeCoreNum();

    is_tree_ = true;
    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        if (GetCoreNum(v) > 1) {
            is_tree_ = false;
            break;
        }
    }

    ExtractResidualStructure();
    // sort back
    for (auto &it : adj_list) {
        std::sort(it.begin(), it.end());
    }


    all_incident_edges_.resize(num_vertex_);
    incident_edges_.resize(num_vertex_);
    for (int i = 0; i < GetNumVertices(); i++) {
        incident_edges_[i].resize(data.GetNumLabels());
    }
    int edge_id = 0;
    for (int i = 0; i < GetNumVertices(); i++) {
        for (int neighbor : adj_list[i]) {
            edge_index_map_[{i, neighbor}] = edge_id;
//                fprintf(stderr,"Query Edge {%u-%u} is %d\n", i, neighbor, edge_index_map_[{i, neighbor}]);
            EdgeInfo e;
            e.edge_label = edge_labels_[{i, neighbor}];
            e.vsum = i + neighbor;
            e.to = neighbor;
            e.index = edge_id;
            e.vp = {i, neighbor};
            edge_info_.push_back(e);
            edge_exists[i][neighbor] = edge_id;
            incident_edges_[i][GetLabel(neighbor)].push_back(edge_id);
            all_incident_edges_[i].push_back(edge_id);
            edge_id++;
        }
    }
    to_.resize(edge_id);
    for (int i = 0; i < edge_id; i++) {
        to_[i] = edge_info_[i].to;
    }

    opposite_edge.resize(edge_id);
    for (int i = 0; i < edge_id; i++) {
        auto p = edge_info_[i].vp;
        opposite_edge[i] = edge_exists[p.second][p.first];
    }

    edge_num_neighbors.resize(edge_info_.size(), std::vector<int>(data.GetNumLabels(), 0));
    four_cycles.resize(edge_info_.size());
    triangles.resize(edge_id, std::vector<Vertex>());
    trigvertex.resize(edge_id, tsl::hopscotch_map<int, std::pair<int, int>>());

    Size num_four_cycles = 0;
    Size num_three_cycles = 0;
    trig_empty.resize(edge_info_.size());
    quad_empty.resize(edge_info_.size());
    num_labeled_triangles_.resize(edge_id, std::vector<int>(data.GetNumLabels()));
    for (EdgeInfo e : edge_info_){
        VertexPair vp = e.vp;
        auto &va = adj_list[vp.first];
        auto &vb = adj_list[vp.second];
        std::set_intersection(va.begin(), va.end(), vb.begin(), vb.end(),
                              std::back_inserter(triangles[e.index]));

        for (int i = 0; i < GetNumLabels(); i++) {
            edge_num_neighbors[e.index][i] += GetIncidentEdges(vp.first, i).size();
            edge_num_neighbors[e.index][i] += GetIncidentEdges(vp.second, i).size();
        }

        for (auto &vc : triangles[e.index]) {
            trigvertex[e.index][vc] = {GetEdgeIndex(vp.first, vc), GetEdgeIndex(vp.second, vc)};
            num_labeled_triangles_[e.index][GetLabel(vc)]++;
            edge_num_neighbors[e.index][GetLabel(vc)]--;
        }
        num_three_cycles += triangles[e.index].size();

        edge_id = edge_index_map_[vp];
    }
    num_four_cycles_indexed = 0;
    four_cycles.resize(edge_info_.size());
    for (int i = 0; i < GetNumVertices(); i++) {
        for (int cand_edge: all_incident_edges_[i]) {
            int nxt_cand = opposite(cand_edge, i);
            for (int third_edge_idx: all_incident_edges_[nxt_cand]) {
                int third_cand = opposite(third_edge_idx, nxt_cand);
                if (third_cand == i) continue;
                for (int opp_edge_idx: all_incident_edges_[third_cand]) {
                    int fourth_cand = opposite(opp_edge_idx, third_cand);
                    if (fourth_cand == nxt_cand) continue;
                    int fourth_edge_idx = GetEdgeIndex(i, fourth_cand);
                    if (fourth_edge_idx != -1) {
                        CycleInfo c_info;
                        c_info.opp_edge_idx = opp_edge_idx;
                        c_info.third_edge_idx = third_edge_idx;
                        c_info.fourth_edge_idx = fourth_edge_idx;
                        c_info.third = third_cand;
                        c_info.fourth = fourth_cand;

                        c_info.one_three_idx = GetEdgeIndex(i, third_cand);
                        c_info.two_four_idx = GetEdgeIndex(nxt_cand, fourth_cand);

                        four_cycles[cand_edge].push_back(c_info);
                        num_four_cycles_indexed++;
                    }
                }
            }
            max_four_cycles_indexed = std::max(max_four_cycles_indexed, (int)four_cycles[cand_edge].size());
        }
    }
    for (int i = 0; i < edge_id; i++) {
        if (triangles[i].empty()) {
            trig_empty[i] = true;
        }
        if (four_cycles[i].empty()) {
            quad_empty[i] = true;
        }
    }
    fprintf(stdout, "QUERY (V, E) = (%u, %u)\n",GetNumVertices(), GetNumEdges());
    fprintf(stdout, "NUM_QUERY_CYCLES = 3[%u] 4[%u]\n",num_three_cycles,num_four_cycles_indexed);
    fflush(stdout);
    return true;
}

bool QueryGraph::LoadAndProcessGraph(const DataGraph &data) {
    LoadRoughGraph(&adj_list);

    max_degree_ = 0;
    num_label_ = 0;
    max_label_ = 0;

    label_frequency_ = new Size[data.GetNumLabels()];
    verticesbyLabel.resize(data.GetNumLabels());
    start_off_ = new Size[GetNumVertices() + 1];
    linear_adj_list_ = new Vertex[GetNumEdges()*2];
    core_num_ = new Size[GetNumVertices()];

    return ProcessLabeledGraph(data);
}

namespace {
    struct NECInfo {
        bool visit = false;
        Vertex representative;
        Size NEC_elems_idx;
    };
}  // namespace

void QueryGraph::ExtractResidualStructure() {
    NECInfo *NEC_infos_temp = new NECInfo[GetNumVertices() * (max_label_ + 1)];
    NEC_elems_ = new NECElement[GetNumVertices()];
    NEC_map_ = new Vertex[GetNumVertices()];
    NEC_size_ = new Size[GetNumVertices()];

    Size num_NEC_elems_ = 0;

    std::fill(NEC_map_, NEC_map_ + GetNumVertices(), INVALID_VTX);
    std::fill(NEC_size_, NEC_size_ + GetNumVertices(), 0);

    num_non_leaf_vertices_ = GetNumVertices();

    // construct NEC map
    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        if (GetDegree(v) == 1) {
            Vertex p = GetNeighbor(GetStartOffset(v));
            Label l = GetLabel(v);

            NECInfo &info = NEC_infos_temp[GetNumVertices() * l + p];
            if (!info.visit) {
                info = {true, v, num_NEC_elems_};
                NEC_map_[v] = v;

                for (Size nbr_idx = GetStartOffset(p); nbr_idx < GetEndOffset(p);
                     ++nbr_idx) {
                    Vertex nbr = GetNeighbor(nbr_idx);
                    if (nbr == v) {
                        NEC_elems_[num_NEC_elems_] = {l, p, v, 0,
                                                      nbr_idx - GetStartOffset(p)};
                        break;
                    }
                }
                num_NEC_elems_ += 1;
            }
            else {
                NEC_map_[v] = info.representative;
            }
            NEC_size_[info.representative] += 1;
            NEC_elems_[info.NEC_elems_idx].size += 1;
            num_non_leaf_vertices_ -= 1;
        }
        else {
            NEC_size_[v] += 1;
        }
    }

    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        NEC_size_[v] = NEC_size_[GetNECRepresentative(v)];
    }

    num_NEC_label_ = 0;
    NEC_start_offs_ = nullptr;
    if (num_NEC_elems_ > 0) {
        // sort NEC elems by label
        std::sort(NEC_elems_, NEC_elems_ + num_NEC_elems_,
                  [](const NECElement &a, const NECElement &b) -> bool {
                      return a.label < b.label;
                  });

        // construct start offsets of NEC elems for same label
        NEC_start_offs_ = new Size[GetNumVertices() + 1];
        NEC_start_offs_[0] = 0;
        num_NEC_label_ += 1;

        Label prev_label = NEC_elems_[0].label;
        for (Size i = 1; i < num_NEC_elems_; ++i) {
            if (NEC_elems_[i].label != prev_label) {
                prev_label = NEC_elems_[i].label;
                NEC_start_offs_[num_NEC_label_] = i;
                num_NEC_label_ += 1;
            }
        }
        NEC_start_offs_[num_NEC_label_] = num_NEC_elems_;
    }

    delete[] NEC_infos_temp;
}
