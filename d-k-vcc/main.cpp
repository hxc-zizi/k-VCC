#include "Graph.h"
#include "Utility.h"
#include "VCC.h"
#include <cstdlib>
#include <iostream>
#include <sys/time.h>
#include <sys/resource.h>
#include <mpi.h>  // ��� MPI ͷ�ļ�

// ɾ�� #include "Test.h"
// ��������
void send_results_to_master(const std::vector<Graph*>& kvccs);
std::vector<Graph*> receive_results_from(int source);
void usage() {
    std::cout << "Usage: program graph_file k" << std::endl;
    exit(1);
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);  // ��ʼ�� MPI

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 3) {
        if (rank == 0) usage();
        MPI_Finalize();
        return 1;
    }

    Utility ut;
    Graph* g = nullptr;

    // �����̼���ͼ����
    if (rank == 0) {
        g = ut.create_graph(string(argv[1]));
    }

    // �㲥ͼ�Ķ������ͱ���
    int vertices_num = 0, edges_num = 0;
	    if (rank == 0) {
		    for (int i = 0; i < g->vertices_num_; ++i) {
		        // ���Ͷ���i���ھ�����
		        int degree = g->vertices_[i]->degree;
		        MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
		        // ���Ͷ���i���ھ��б�
		        MPI_Bcast(g->vertices_[i]->nbr.data(), degree, MPI_INT, 0, MPI_COMM_WORLD);
		    }
		} else {
		    for (int i = 0; i < vertices_num; ++i) {
		        int degree;
		        MPI_Bcast(&degree, 1, MPI_INT, 0, MPI_COMM_WORLD);
		        g->vertices_[i]->nbr.resize(degree);
		        MPI_Bcast(g->vertices_[i]->nbr.data(), degree, MPI_INT, 0, MPI_COMM_WORLD);
		        g->vertices_[i]->degree = degree;  // ���ö������
		        g->vertices_[i]->id = i;           // ���ö���ID
		    }
	}

//    // �ӽ��̳�ʼ��ͼ����
//    if (rank != 0) {
//	    g = new Graph();
//	    g->vertices_num_ = vertices_num;
//	    g->edges_num_ = edges_num;
//	    g->vertices_.resize(vertices_num);
//	    for (int i = 0; i < vertices_num; ++i) {
//	        g->vertices_[i] = new Vertex();
//	    }
//    }

    // �����̷ָ�ͼ��������ͼ���ӽ���
    vector<int> local_vertices;
	// �����̷�����ͼ����
	// �޸� main.cpp �е� MPI ͨ�Ų���
	if (rank == 0) {
	    vector<vector<int>> subgraphs = ut.split_graph(g, size);
	    for (int i = 1; i < size; i++) {
	        int num_vertices = subgraphs[i].size();
	        MPI_Send(&num_vertices, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
	        
	        // ��Ӽ�飺ȷ���������ݲ�Ϊ��
	        if (num_vertices > 0) {
	            MPI_Send(subgraphs[i].data(), num_vertices, MPI_INT, i, 0, MPI_COMM_WORLD);
	        } else {
	            fprintf(stderr, "WARNING: sending empty subgraph to process %d\n", i);
	        }
	    }
	} else {
	    int num_vertices;
	    MPI_Recv(&num_vertices, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    
	    // ��Ӽ�飺��������ݽ���
	    if (num_vertices < 0) {
	        fprintf(stderr, "ERROR: received invalid vertex count %d\n", num_vertices);
	        MPI_Abort(MPI_COMM_WORLD, 1);
	    }
	    
	    local_vertices.resize(num_vertices);
	    if (num_vertices > 0) {
	        MPI_Recv(local_vertices.data(), num_vertices, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	    }
	}	

// �����ֲ�ͼ������ԭʼͼָ�� g��
	Graph* local_graph = ut.create_local_graph(local_vertices, g);

    // �����̶������� k-VCC
    VCC vc(local_graph);
    vc.get_all_kvccs(atoi(argv[2]));
    
    

    // �ռ������������
    if (rank != 0) {
        send_results_to_master(vc.kvccs_);
    } else {
        std::vector<Graph*> global_result;
       // �ռ����ʱֱ�Ӵ洢 Graph* ָ��
		for (auto& kvcc : vc.kvccs_) {
		    global_result.push_back(kvcc);
		}
        for (int i = 1; i < size; i++) {
            auto remote_kvccs = receive_results_from(i);
            global_result.insert(global_result.end(), remote_kvccs.begin(), remote_kvccs.end());
        }
        printf("Total %lu k-VCCs\n", global_result.size());
    }

    ut.release_graph(g);
    MPI_Finalize();
    return 0;
}

// ����ʵ��
void send_results_to_master(const std::vector<Graph*>& kvccs) {
    int num_kvccs = kvccs.size();
    MPI_Send(&num_kvccs, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);  // ���� k-VCC ����

    for (auto& graph : kvccs) {
        // ����ÿ�� Graph ����Ķ�������
        int num_vertices = graph->vertices_num_;
        MPI_Send(&num_vertices, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

        // ����ÿ�� Graph ����Ķ��� ID
        for (auto& vertex : graph->vertices_) {
            int vertex_id = vertex->id;
            MPI_Send(&vertex_id, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }
}

std::vector<Graph*> receive_results_from(int source) {
    std::vector<Graph*> remote_kvccs;

    // ���� k-VCC ����
    int num_kvccs;
    MPI_Recv(&num_kvccs, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < num_kvccs; i++) {
        // ����ÿ�� Graph ����Ķ�������
        int num_vertices;
        MPI_Recv(&num_vertices, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // �����µ� Graph ����
        Graph* graph = new Graph();
        graph->vertices_num_ = num_vertices;

        // ����ÿ�� Graph ����Ķ��� ID
        for (int j = 0; j < num_vertices; j++) {
            int vertex_id;
            MPI_Recv(&vertex_id, 1, MPI_INT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            Vertex* vertex = new Vertex();
            vertex->id = vertex_id;
            graph->vertices_.push_back(vertex);
        }

        remote_kvccs.push_back(graph);
    }

    return remote_kvccs;
}
