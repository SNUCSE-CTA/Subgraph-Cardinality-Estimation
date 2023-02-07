#include "include/daf_query_graph.h"
using namespace daf;
QueryGraph::QueryGraph(const std::string &filename) : Graph(filename) {
    myname = filename;
}

QueryGraph::~QueryGraph() {}

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

    Size num_three_cycles = 0;
    num_labeled_triangles_.resize(edge_id, std::vector<int>(data.GetNumLabels()));
    IndexTriangles();
    num_four_cycles_indexed = 0;
    IndexFourCycles();
    for (int i = 0; i < edge_id; i++) {
        num_three_cycles += local_triangles[i].size();
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