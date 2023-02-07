#include "include/daf_graph.h"
using namespace daf;
using std::vector, std::string, std::deque;
deque<string> parse(string line, string del) {
    deque<string> ret;

    size_t pos = 0;
    string token;
    while((pos = line.find(del)) != string::npos)
    {
        token = line.substr(0, pos);
        ret.push_back(token);
        line.erase(0, pos + del.length());
    }
    ret.push_back(line);
    return ret;
}
Graph::Graph(const std::string &filename)
        : filename_(filename), fin_(filename) {}

Graph::~Graph() {
    delete[] start_off_;
    delete[] linear_adj_list_;
    delete[] label_;
    delete[] label_frequency_;
    delete[] core_num_;
}

void Graph::LoadRoughGraph(std::vector<std::vector<Vertex>> *graph) {
    if (!fin_.is_open()) {
        std::cerr << "Graph file " << filename_ << " not found!\n";
        exit(EXIT_FAILURE);
    }
    std::cerr << "Start reading graph file from " << filename_ << " " << fileSize(filename_.c_str()) << std::endl;

    Size v, e;
    std::string type;
    std::string ignore;

    fin_ >> type >> v >> e;

    num_vertex_ = v;
    num_edge_ = e;
    edge_exists.resize(num_vertex_);
    degrees.resize(num_vertex_);
    label_ = new Label[v];
    graph->resize(v);

    std::string line;
    // preprocessing
    Label max_vertex_label_ = 0;
    int num_lines = 0;
    while (getline(fin_, line)) {
        auto tok = parse(line, " ");
        type = tok[0];
        tok.pop_front();
        if (type[0] == 'v') {
            Vertex id = std::stoi(tok.front());
            tok.pop_front();
            Label l;
            if (tok.empty()) l = 0;
            else {
                l = std::stoi(tok.front());
                tok.pop_front();
            }
            label_[id] = l;
            max_vertex_label_ = std::max(max_vertex_label_, label_[id]);
        } else if (type[0] == 'e') {
            Vertex v1, v2;
            v1 = std::stoi(tok.front()); tok.pop_front();
            v2 = std::stoi(tok.front()); tok.pop_front();
            (*graph)[v1].push_back(v2);
            (*graph)[v2].push_back(v1);
            all_edges.push_back({v1, v2});
            all_edges.push_back({v2, v1});
            degrees[v1]++; degrees[v2]++;
            Label el;
            if (tok.empty()) el = 0;
            else {
                el = std::stoi(tok.front()); tok.pop_front();
            }
            edge_labels_[{v1, v2}] = edge_labels_[{v2, v1}] = el;
        }
        num_lines++;
    }

    fprintf(stderr, "Graph Read finish : Read %u vertices and %u edges\n",num_vertex_,num_edge_);
    fflush(stderr);
    fin_.close();
}


void Graph::computeCoreNum() {
    Size *bin = new Size[max_degree_ + 1];
    Size *pos = new Size[GetNumVertices()];
    Vertex *vert = new Vertex[GetNumVertices()];

    std::fill(bin, bin + (max_degree_ + 1), 0);

    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        bin[core_num_[v]] += 1;
    }

    Size start = 0;
    Size num;

    for (Size d = 0; d <= max_degree_; ++d) {
        num = bin[d];
        bin[d] = start;
        start += num;
    }

    for (Vertex v = 0; v < GetNumVertices(); ++v) {
        pos[v] = bin[core_num_[v]];
        vert[pos[v]] = v;
        bin[core_num_[v]] += 1;
    }

    for (Size d = max_degree_; d--;) bin[d + 1] = bin[d];
    bin[0] = 0;

    for (Size i = 0; i < GetNumVertices(); ++i) {
        Vertex v = vert[i];

        for (Size j = GetStartOffset(v); j < GetEndOffset(v); j++) {
            Vertex u = GetNeighbor(j);

            if (core_num_[u] > core_num_[v]) {
                Size du = core_num_[u];
                Size pu = pos[u];

                Size pw = bin[du];
                Vertex w = vert[pw];

                if (u != w) {
                    pos[u] = pw;
                    pos[w] = pu;
                    vert[pu] = w;
                    vert[pw] = u;
                }

                bin[du]++;
                core_num_[u]--;
            }
        }
    }

    delete[] bin;
    delete[] pos;
    delete[] vert;
}

void Graph::IndexTriangles() {
    vertex_local_triangles.resize(GetNumVertices(), 0);
    local_triangles.resize(edge_info_.size());
    int cnt = 0;
    std::vector<int> common_neighbor(GetNumVertices(), -1);
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
                        TrigInfo t;
                        t.point = snd_neighbor;
                        vertex_local_triangles[t.point]++;
                        t.fst_edge = common_neighbor[snd_neighbor];
                        t.snd_edge = snd_incident;
                        local_triangles[edge_exists[i][neighbor]].push_back(t);
                        cnt++;
                        num_labeled_triangles_[edge_exists[i][neighbor]][GetLabel(snd_neighbor)]++;
                        edge_num_neighbors[edge_exists[i][neighbor]][GetLabel(snd_neighbor)]--;
                    }
                }
                for (int fst_incident : GetIncidentEdges(i, l)) {
                    int fst_neighbor = opposite(fst_incident, i);
                    common_neighbor[fst_neighbor] = -1;
                }
            }
        }
    }
}

void Graph::IndexFourCycles() {
    four_cycles.resize(edge_info_.size());
    for (int i = 0; i < GetNumVertices(); i++) {
        for (int cand_edge : all_incident_edges_[i]) {
            int nxt_cand = opposite(cand_edge, i);
            int indexed_cnt = 0;
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
                        indexed_cnt++;
                    }
                }
            }
            num_four_cycles_indexed += indexed_cnt;
            max_four_cycles_indexed = std::max(max_four_cycles_indexed, (int)four_cycles[cand_edge].size());
        }
    }
}
