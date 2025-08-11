#ifndef _UTILITY_H_
#define _UTILITY_H_

#include "Graph.h"
#include <mpi.h>  // 添加 MPI 头文件

#include <queue>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <string>
#include <cstring>
#include <algorithm>

#define _LINUX_

#ifdef _LINUX_
    #include <sys/time.h>
    #include <unistd.h>
    #include <sys/types.h>  
    #include <sys/stat.h>
    #include <sys/resource.h>
#else
    #include <io.h>
    #include <direct.h>
#endif

using namespace std;

FILE *open_file(const char *file_name, const char *mode);

const int MAX_ID = 1000000000;

// struct used in process .txt file in class Utility
struct Arc {
    int a;
    int b;
};

class Utility {
    int* vertex_map_;
    string dir_;
    static bool edge_compare(const Arc &e1, const Arc &e2);
    bool is_merge_finished(Arc* es, int size);
    void save_tmp_edges(Arc* edges, int size, int tmp_file);
    void merge(int size);
    int get_min_edge(Arc* es, int size);
    int get_vertex_id(int u, int &num);

public:
    Utility();
    ~Utility();

    Graph* create_graph(string txt_file_name);
    Graph* load_from_binary_graph(string binary_file_dir);
    void format_graph(string binary_file_directory, string txt_file_name);
    void format_graph(string binary_file_directory, string txt_file_name, int percent, bool is_sample_edge);
    int release_graph(Graph* g);

// 添加以下函数声明
    std::vector<std::vector<int>> split_graph(Graph* g, int num_parts);
    Graph* create_local_graph(const std::vector<int>& vertices, Graph* old_graph);  // 添加参数
};

#endif
