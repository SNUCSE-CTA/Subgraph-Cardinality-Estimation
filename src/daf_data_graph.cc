#include "include/daf_data_graph.h"
#include <cmath>
using namespace daf;
DataGraph::DataGraph(const std::string &filename) : Graph(filename) {}

DataGraph::~DataGraph() {
    delete[] true_label_;
    delete[] transferred_label_;
    delete[] adj_offs_by_label_;
    delete[] offs_by_label_;
    delete[] vertices_sorted_;
    delete[] linear_nbr_bitset_;
    delete[] max_nbr_degree_;
}

void DataGraph::ProceesLabeledGraph() {
    // transfer label & get sorted degrees (for constructing C_ini(u))
    Label max_label = 0;
    Label cur_transferred_label = 0;
    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        vertices_sorted_[v] = v;
        Label l = label_[v];
        if (l > max_label) max_label = l;
        if (transferred_label_map.find(l) == transferred_label_map.end()) {
            transferred_label_map[l] = cur_transferred_label;
            cur_transferred_label += 1;
        }
        label_[v] = transferred_label_map[l];
    }

    std::sort(vertices_sorted_, vertices_sorted_ + GetNumVertices(),
              [this](Vertex v1, Vertex v2) -> bool {
                  if (GetLabel(v1) != GetLabel(v2)) {
                      return GetLabel(v1) < GetLabel(v2);
                  }
                  else {
                      return adj_list[v1].size() > adj_list[v2].size();
                  }
              });

    num_label_ = transferred_label_map.size();

    label_frequency_ = new Size[num_label_];
    std::fill(label_frequency_, label_frequency_ + num_label_, 0);

    max_label_frequency_ = 0;

    transferred_label_ = new Label[max_label + 1];
    // transferred_label_[l] = INVALID_LB iff there is no label l in data graph
    std::fill(transferred_label_, transferred_label_ + max_label + 1, INVALID_LB);

    for (auto p: transferred_label_map) {
        transferred_label_[p.first] = p.second;
    }

    nbr_bitset_size_ = (GetNumLabels() - 1) / (sizeof(uint64_t) * CHAR_BIT) + 1;
    linear_nbr_bitset_ = new uint64_t[GetNumVertices() * nbr_bitset_size_];
    std::fill(linear_nbr_bitset_,
              linear_nbr_bitset_ + GetNumVertices() * nbr_bitset_size_, 0ull);

    max_nbr_degree_ = new Size[GetNumVertices()];
    std::fill(max_nbr_degree_, max_nbr_degree_ + GetNumVertices(), 0);

    // compute offsets & construct adjacent list and label frequency
    Size cur_idx = 0;
    max_degree_ = 0;
    linear_adj_list_ = new Vertex[GetNumEdges()*2];
    start_off_ = new Size[GetNumVertices() + 1];
    offs_by_label_ = new Size[GetNumLabels() + 1];
    adj_offs_by_label_ =
            new std::pair<Size, Size>[GetNumVertices() * GetNumLabels()];
    core_num_ = new Size[GetNumVertices()];

    Label cur_label = 0;  // min label of data graph is 0
    offs_by_label_[0] = 0;
    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        Size start = v * GetNumLabels();
        label_frequency_[GetLabel(v)] += 1;
        start_off_[v] = cur_idx;

        if (label_frequency_[GetLabel(v)] > max_label_frequency_) {
            max_label_frequency_ = label_frequency_[GetLabel(v)];
        }

        Label label_sorted = GetLabel(vertices_sorted_[v]);
        if (label_sorted != cur_label) {
            offs_by_label_[label_sorted] = v;
            cur_label = label_sorted;
        }

        // initialize core number
        core_num_[v] = adj_list[v].size();
        if (adj_list[v].size() > max_degree_) max_degree_ = adj_list[v].size();

        if (adj_list[v].size() == 0) {
            continue;
        }

        // sort by label first and degree second
        std::sort(adj_list[v].begin(), adj_list[v].end(),
                  [this](Vertex v1, Vertex v2) -> bool {
                      if (GetLabel(v1) != GetLabel(v2)) {
                          return GetLabel(v1) < GetLabel(v2);
                      }
                      else {
                          return adj_list[v1].size() > adj_list[v2].size();
                      }
                  });

        // compute adj(v, l).start, adj(v, l).end for all l in nbr_label(v)
        Label cur_adj_label = GetLabel(adj_list[v][0]);
        adj_offs_by_label_[start + cur_adj_label].first = cur_idx;
        max_nbr_degree_[v] = adj_list[adj_list[v][0]].size();
        for (Size i = 1; i < adj_list[v].size(); ++i) {
            if (cur_adj_label != GetLabel(adj_list[v][i])) {
                linear_nbr_bitset_[nbr_bitset_size_ * v +
                                   (cur_adj_label / (sizeof(uint64_t) * CHAR_BIT))] |=
                        1ull << (cur_adj_label % (sizeof(uint64_t) * CHAR_BIT));
                if (max_nbr_degree_[v] < adj_list[adj_list[v][i]].size()) {
                    max_nbr_degree_[v] = adj_list[adj_list[v][i]].size();
                }
                adj_offs_by_label_[start + cur_adj_label].second = cur_idx + i;
                cur_adj_label = GetLabel(adj_list[v][i]);
                adj_offs_by_label_[start + cur_adj_label].first = cur_idx + i;
            }
        }
        linear_nbr_bitset_[nbr_bitset_size_ * v +
                           (cur_adj_label / (sizeof(uint64_t) * CHAR_BIT))] |=
                1ull << (cur_adj_label % (sizeof(uint64_t) * CHAR_BIT));

        adj_offs_by_label_[start + cur_adj_label].second =
                cur_idx + adj_list[v].size();

        // copy adj_list to linear_adj_list
        std::copy(adj_list[v].begin(), adj_list[v].end(),
                  linear_adj_list_ + cur_idx);

        cur_idx += adj_list[v].size();
    }
    start_off_[GetNumVertices()] = num_edge_ * 2;
    offs_by_label_[GetNumLabels()] = GetNumVertices();

    // preprocess for data graph
    computeCoreNum();
}

void DataGraph::LoadAndProcessGraph() {
    LoadRoughGraph(&adj_list);
    std::cout << "Datagraph processing..." << std::endl;
    true_label_ = new Label[GetNumVertices()];
    vertices_sorted_ = new Vertex[GetNumVertices()];
    ProceesLabeledGraph();
    std::cout << "Process_finish_incidence_check..." << std::endl;

    // sort back
    for (auto &it : adj_list) {
        std::sort(it.begin(), it.end());
    }

    all_incident_edges_.resize(num_vertex_);
    incident_edges_.resize(num_vertex_);
    for (int i = 0; i < GetNumVertices(); i++) {
        incident_edges_[i].resize(GetNumLabels());
    }
    int edge_id = 0;
    for (int i = 0; i < GetNumVertices(); i++) {
        for (int neighbor : adj_list[i]) {
            edge_index_map_[{i, neighbor}] = edge_id;
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
    opposite_edge.resize(edge_id);
    for (int i = 0; i < GetNumLabels(); i++) {
        for (auto &eid : all_incident_edges_[i]) {
            int to = edge_info_[eid].to;
            opposite_edge[eid] = edge_exists[to][i];
        }
    }
    to_.resize(edge_id);
    for (int i = 0; i < edge_id; i++) {
        to_[i] = edge_info_[i].to;
    }
    for (int i = 0; i < GetNumVertices(); i++) {
        for (auto &vec : incident_edges_[i]) {
            std::sort(vec.begin(), vec.end(),[this, i](auto &a, auto &b) -> bool {
                int opp_a = edge_info_[a].vsum-i;
                int opp_b = edge_info_[b].vsum-i;
                return adj_list[opp_a].size() > adj_list[opp_b].size();
            });
        }
    }

    label_edge_offset.resize(GetNumVertices(), std::vector<int>(GetNumLabels()+1));
    for (int cand = 0; cand < GetNumVertices(); cand++) {
        for (int i = 0; i < GetNumLabels(); i++) {
            label_edge_offset[cand][i+1] = GetIncidentEdges(cand, i).size();
            if (i > 0) label_edge_offset[cand][i] += label_edge_offset[cand][i-1];
        }
    }




    std::cout << "Incidence_ok! : ";
    if (is_sparse()) {
        num_labeled_triangles_.resize(edge_id, std::vector<int>(GetNumLabels()));
        trigvertex.resize(edge_id, tsl::hopscotch_map<int, std::pair<int, int>>());
        edge_num_neighbors.resize(edge_info_.size(), std::vector<int>(GetNumLabels(), 0));
        std::cout << "Working on triangles...." << std::endl;
        std::vector<int> common_neighbor(GetNumVertices(), -1);
        int trig_count = 0;
        for (int i = 0; i < GetNumVertices(); i++) {
            for (int neighbor : adj_list[i]) {
                for (int l = 0; l < GetNumLabels(); l++) {
                    for (int fst_incident : GetIncidentEdges(i, l)) {
                        int fst_neighbor = opposite(fst_incident, i);
                        common_neighbor[fst_neighbor] = fst_incident;
                    }
                    for (int snd_incident : GetIncidentEdges(neighbor, l)) {
                        int snd_neighbor = opposite(snd_incident, neighbor);
                        if (common_neighbor[snd_neighbor] != -1) {
                            trigvertex[edge_exists[i][neighbor]][snd_neighbor] = {common_neighbor[snd_neighbor], snd_incident};
                            num_labeled_triangles_[edge_exists[i][neighbor]][GetLabel(snd_neighbor)]++;
                            edge_num_neighbors[edge_exists[i][neighbor]][GetLabel(snd_neighbor)]--;
                            trig_count++;
                        }
                    }
                    for (int fst_incident : GetIncidentEdges(i, l)) {
                        int fst_neighbor = opposite(fst_incident, i);
                        common_neighbor[fst_neighbor] = -1;
                    }
                }
            }
            if (i % 50000 == 0) {
                fprintf(stderr, "%d / %d found %d trigs\n",i,GetNumVertices(), trig_count);
            }
        }

        for (EdgeInfo e : edge_info_) {
            VertexPair vp = e.vp;
            for (int i = 0; i < GetNumLabels(); i++) {
                edge_num_neighbors[e.index][i] += GetIncidentEdges(vp.first, i).size();
                edge_num_neighbors[e.index][i] += GetIncidentEdges(vp.second, i).size();
            }
            max_num_trigs = std::max(max_num_trigs, (int)trigvertex[e.index].size());
        }

        std::cout << "Total # of trigs : " << trig_count << std::endl;

        four_cycles.resize(edge_id);

        double num_4c_cand = 0.0;
        for (int i = 0; i < GetNumVertices(); i++) {
            for (int cand_edge: all_incident_edges_[i]) {
                int nxt_cand = opposite(cand_edge, i);
                num_4c_cand += all_incident_edges_[i].size() * all_incident_edges_[nxt_cand].size();
            }
        }
        if (num_4c_cand < 1e10) {
            for (int i = 0; i < GetNumVertices(); i++) {
                for (int cand_edge : all_incident_edges_[i]) {
                    int nxt_cand = opposite(cand_edge, i);
                    for (int third_edge_idx : all_incident_edges_[nxt_cand]) {
                        int third_cand = opposite(third_edge_idx, nxt_cand);
                        if (third_cand == i) continue;
                        for (int opp_edge_idx : all_incident_edges_[third_cand]) {
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
                }
                if (i % 30000 == 0) {
                    fprintf(stderr, "%d / %d found %d quads\n",i,GetNumVertices(), num_four_cycles_indexed);
                    fflush(stderr);
                }
            }
            std::cout << "Total # of quads : " << num_four_cycles_indexed << std::endl;
        }
        else {
            std::cout << "Too many possible 4-Cycles...abandon 4cycle preprocessing" << std::endl;
        }
    }
    else {
        std::cout << "Graph is not sparse...abandon triangle preprocessing" << std::endl;
    }
    std::cout << "Finished...Compute Entropy..." << std::endl;

    // Compute Data Label Entropy
    std::vector<double> label_probability;
    for (Label i = 0; i < GetNumLabels(); i++) {
        label_probability.push_back(GetLabelFrequency(i) * 1.0 / GetNumVertices());
    }
    double ent = 0.0;
    for (auto x : label_probability) {
        ent -= x * log2(x);
    }
    fprintf(stderr, "Node Label Entropy = %.02lf bits\n", ent);
    std::map<Label, double> edge_label_probability;
    for (auto &p : edge_labels_) {
        edge_label_probability[p.second]+=1;
    }
    for (auto &p : edge_label_probability) {
        p.second /= num_edge_;
    }
    ent = 0.0;
    for (auto x : edge_label_probability) {
        ent -= x.second * log2(x.second);
    }
    fprintf(stderr, "Edge Label Entropy = %.02lf bits\n", ent);
}

