#ifndef VCC_H_
#define VCC_H_

#include "Graph.h"
#include "BitMap.h"
#include <vector>
#include <string>
#include <queue>

class F_edge {
public:
    F_edge();
    F_edge(const int& a, const int& b, const int& capacity, F_edge* reverse);
    ~F_edge();
    int a_;
    int b_;
    int capacity_;
    F_edge* reverse_;
};

class F_vertex {
public:
    F_vertex();
    F_vertex(const int& id);
    ~F_vertex();
    int id_;
    std::vector<F_edge*> nbr_;
};

class VCC {
public:
    std::vector<Graph*> kvccs_;  // 设为公有
    Graph* g_;

    // 成员函数声明
    VCC();
    VCC(Graph* g);
    ~VCC();
    void print_kvccs();
    void write_kvccs(const std::string& dir, const int& k);
    void get_all_kvccs(const int& k);
    void get_all_kvccs(Graph* g, const int& k);
    int release_graph(Graph* g);
    std::vector<std::vector<int>> get_kcores(Graph* g, const int& k);
    Graph* create_directed_flow_graph(Graph* g);
    std::vector<F_vertex*> create_flow_graph(Graph* g);
    int get_source_id(Graph* g);
    int flow(const int& source_id, const int& target_id, std::vector<F_vertex*>& f_graph, const int& k);
    bool label_graph(const int& source_id, const int& target_id, std::vector<F_vertex*>& f_graph, int* label, std::vector<int>& label_count);
    int find_augmenting_path(const int& source_id, const int& target_id, std::vector<F_vertex*>& f_graph, int* label, std::vector<int>& label_count, BitMap& pru);
    std::vector<Graph*> cut_graph(const int& source_id, const int& target_id, std::vector<F_vertex*>& f_graph, Graph* g);
    std::vector<int> create_distance_order(const int& source, Graph* g);
    std::vector<Graph*> global_cut(Graph* g, const int& k);
    std::vector<Graph*> local_cut(const int& source_id, const int& target_id, Graph* g, Graph* df_graph, const int& k);
    Graph* get_sparse_certificate(Graph* g, std::vector<std::vector<int>>& groups, int* vertex_to_group, const int& k);
    void scan_first_search(const int& root, Graph* g, Graph* sparse_certificate, BitMap& mark, std::vector<std::vector<int>>& groups, int* vertex_to_group, const int& k);
    void sweep(const int& v, BitMap& pru, int* deposit, int* vertex_to_group, int* group_deposit, std::vector<std::vector<int>>& groups, Graph* g, const int& k);
    std::vector<std::vector<int>> get_all_connected_components(Graph* g);
    Graph* create_new_graph(const std::vector<int>& component, Graph* old_graph);
    Vertex* create_new_vertex(const int& old_id, Graph* old_graph, Graph* new_graph);
};

#endif // VCC_H_
