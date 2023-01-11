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
    }

//    fprintf(stderr, "Graph Read finish : Read %u vertices and %u edges\n",num_vertex_,num_edge_);
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
