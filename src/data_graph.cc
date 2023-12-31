#include "include/data_graph.h"
#include <cmath>
using namespace CardEst;
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

    Size cur_idx = 0;
    max_degree_ = 0;
    linear_adj_list_ = new Vertex[GetNumEdges()*2];
    start_off_ = new Size[GetNumVertices() + 1];
    offs_by_label_ = new Size[GetNumLabels() + 1];
    adj_offs_by_label_ =
            new std::pair<Size, Size>[GetNumVertices() * GetNumLabels()];
    core_num_ = new Size[GetNumVertices()];

    Label cur_label = 0;  
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

        
        core_num_[v] = adj_list[v].size();
        if (adj_list[v].size() > max_degree_) max_degree_ = adj_list[v].size();

        if (adj_list[v].size() == 0) {
            continue;
        }

        
        std::sort(adj_list[v].begin(), adj_list[v].end(),
                  [this](Vertex v1, Vertex v2) -> bool {
                      if (GetLabel(v1) != GetLabel(v2)) {
                          return GetLabel(v1) < GetLabel(v2);
                      }
                      else {
                          return adj_list[v1].size() > adj_list[v2].size();
                      }
                  });

        
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

        
        std::copy(adj_list[v].begin(), adj_list[v].end(),
                  linear_adj_list_ + cur_idx);

        cur_idx += adj_list[v].size();
    }
    start_off_[GetNumVertices()] = num_edge_ * 2;
    offs_by_label_[GetNumLabels()] = GetNumVertices();

    
    computeCoreNum();
}

void DataGraph::LoadAndProcessGraph() {
    LoadRoughGraph(&adj_list);
    std::cout << "Datagraph processing..." << std::endl;
    true_label_ = new Label[GetNumVertices()];
    vertices_sorted_ = new Vertex[GetNumVertices()];
    ProceesLabeledGraph();
    
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
    to_.resize(edge_id);
    for (int i = 0; i < edge_id; i++) {
        to_[i] = edge_info_[i].to;
    }

    opposite_edge.resize(edge_id);
    for (int i = 0; i < edge_id; i++) {
        auto p = edge_info_[i].vp;
        opposite_edge[i] = edge_exists[p.second][p.first];
    }

    for (int i = 0; i < GetNumVertices(); i++) {
        for (auto &vec : incident_edges_[i]) {
            std::sort(vec.begin(), vec.end(),[this, i](auto &a, auto &b) -> bool {
                int opp_a = edge_info_[a].vsum-i;
                int opp_b = edge_info_[b].vsum-i;
                return adj_list[opp_a].size() > adj_list[opp_b].size();
            });
        }
        std::sort(all_incident_edges_[i].begin(), all_incident_edges_[i].end(), [this](auto &a, auto &b) -> bool {
            return adj_list[to_[a]].size() > adj_list[to_[b]].size();
        });
    }

    if (is_sparse()) {
        num_labeled_triangles_.resize(edge_id, std::vector<int>(GetNumLabels()));

        edge_num_neighbors.resize(edge_info_.size(), std::vector<int>(GetNumLabels(), 0));
        std::cout << "Enumerating Triangles... ";
        int trig_count = 0;
        IndexTriangles();
        max_num_trigs = 0;
        for (EdgeInfo e : edge_info_) {
            VertexPair vp = e.vp;
            for (int i = 0; i < GetNumLabels(); i++) {
                edge_num_neighbors[e.index][i] += GetIncidentEdges(vp.first, i).size();
                edge_num_neighbors[e.index][i] += GetIncidentEdges(vp.second, i).size();
            }
            trig_count += local_triangles[e.index].size();
            max_num_trigs = std::max(max_num_trigs, (int)local_triangles[e.index].size());
        }

        std::cout << trig_count << " found" << std::endl;

        double num_4c_cand = 0.0;
        for (int i = 0; i < GetNumVertices(); i++) {
            for (int cand_edge: all_incident_edges_[i]) {
                int nxt_cand = opposite(cand_edge, i);
                num_4c_cand += all_incident_edges_[i].size() * all_incident_edges_[nxt_cand].size();
            }
        }
        num_four_cycles_indexed = 0;
        if (num_4c_cand < 1e10) {
            std::cout << "Enumerating Four-Cycles... ";
            IndexFourCycles();
            std::cout << num_four_cycles_indexed << " found" << std::endl;
        }
        else {
            std::cout << "Too many possible 4-Cycles...abandon 4cycle preprocessing" << std::endl;
        }
    }
    else {
        std::cout << "Graph is not sparse...abandon triangle preprocessing" << std::endl;
    }
}

