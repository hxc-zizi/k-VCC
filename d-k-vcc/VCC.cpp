#include "VCC.h"
#include <algorithm>

// 确保所有函数的定义与声明一致
F_edge::F_edge() {
    capacity_ = 0;
    reverse_ = NULL;
}

F_edge::F_edge(const int& a, const int& b, const int& capacity, F_edge* reverse) {
    a_ = a;
    b_ = b;
    capacity_ = capacity;
    reverse_ = reverse;
}

F_edge::~F_edge() {}

F_vertex::F_vertex() {}

F_vertex::F_vertex(const int& id) {
    id_ = id;
}

F_vertex::~F_vertex() {
    for (F_edge* fe : nbr_) delete fe;
}

VCC::VCC() {}

VCC::VCC(Graph* g) : g_(g) {}

VCC::~VCC() {
    if (!kvccs_.empty()) {
        for (Graph* g : kvccs_) {
            if (g != nullptr) {
                delete g;
                g = nullptr;
            }
        }
    }
}

void VCC::print_kvccs(){
	int i = 0;
	for(Graph* g : kvccs_){
		printf("k_VCC %d:\n", ++i);
		for(Vertex* v : g->vertices_){
			printf("%d  ",g->id_new_to_old_[v->id] );
		}
		printf("\n");
	}
}

void VCC::write_kvccs(const std::string& dir, const int& k) {
	printf("Write result in the file\n");
	char path[100];
	sprintf(path,"%s%d_vcc_result.txt",dir.c_str(),k);

	FILE* f_vccs = fopen(path,"w");


	int i = 0;
	for(Graph* g : kvccs_){
		fprintf(f_vccs, "%d_vcc %d:\n",k,++i);
		
		for(Vertex* v : g->vertices_){
			fprintf(f_vccs, "%d, ",g->id_new_to_old_[v->id] );
		}
		fprintf(f_vccs, "\n\n" );
	}

	
	fclose(f_vccs);

	printf("%s created!\n",path);
}

void VCC::get_all_kvccs(const int& k) {
	get_all_kvccs(g_,k);

	printf("Obtain %lu %d_VCCs\n",kvccs_.size(),k);
	// print_kvccs();
}

void VCC::get_all_kvccs(Graph* g, const int& k) {
	
	// Core c(g_);
	vector<vector<int> > ccs = get_kcores(g,k);


	// ECC ec(g_);
	// vector<vector<int> > keccs = c.get_all_kcores(k);


	// vector<vector<int> > ccs = get_all_connected_components(g);


	if (ccs.empty()) {
        printf("No k-core found for k = %d\n", k);  // 打印未找到 k-core 的提示
        return;
    }

    for (auto cc : ccs) {
//        printf("Processing k-core with size: %lu\n", cc.size());  // 打印当前处理的 k-core 大小
		if (cc.empty()) {  // 跳过空的k-core
            printf("Skipping empty k-core\n");
            continue;
        }
        Graph* subgraph = create_new_graph(cc, g);
        if (!subgraph || subgraph->vertices_num_ == 0) { // 检查子图有效性
            printf("Failed to create valid subgraph\n");
            continue;
        }
        vector<Graph*> vs = global_cut(subgraph, k);

        if (vs.empty()) {
            printf("Found k-VCC with size: %lu\n", subgraph->vertices_num_);  // 打印找到的 k-VCC 大小
            kvccs_.push_back(subgraph);
        } else {
            get_all_kvccs(vs[0], k);
            release_graph(vs[0]);
            get_all_kvccs(vs[1], k);
            release_graph(vs[1]);
        }

        release_graph(subgraph);
    }

	
}

vector<vector<int>> VCC::get_kcores(Graph* g, const int& k) {
    int n = g->vertices_num_;
    printf("Number of vertices: %d\n", n);  // 打印顶点数量--检验 
    int* deg = new int[n];
    for (int i = 0; i < n; ++i) {
        deg[i] = g->vertices_[i]->degree;
        printf("Vertex %d (global ID %d) has degree %d\n", 
           i, g->id_new_to_old_[i], deg[i]);  // 检验 打印全局ID和度数 
    }

    BitMap mark(n);
    queue<int> q;

    // 初始化队列：将所有度数 <k 的顶点加入队列
    for (int i = 0; i < n; ++i) {
        if (deg[i] < k) {
            q.push(i);
            mark.set(i);
            printf("Vertex %d added to queue (degree < %d)\n", i, k);  // 打印加入队列的顶点
        }
    }

    // 迭代删除度数不足的顶点
    while (!q.empty()) {
        int u = q.front();
        q.pop();
        printf("Processing vertex %d\n", u);  // 打印正在处理的顶点

        // 遍历邻居，减少其度数
        for (int v : g->vertices_[u]->nbr) {
            if (mark.get(v)) continue; // 已标记的顶点不再处理
            deg[v]--;
            printf("Vertex %d degree reduced to %d\n", v, deg[v]);  // 打印邻居顶点度数更新
            if (deg[v] < k) {
                q.push(v);
                mark.set(v);
                printf("Vertex %d added to queue (degree < %d)\n", v, k);  // 打印加入队列的顶点
            }
        }
    }

    // 收集所有未标记的顶点（即 k-core）
    vector<vector<int>> kcores;
    for (int i = 0; i < n; ++i) {
        if (!mark.get(i)) {
            vector<int> core;
            queue<int> bfs_q;
            bfs_q.push(i);
            mark.set(i); // 临时标记为已访问

            while (!bfs_q.empty()) {
                int u = bfs_q.front();
                bfs_q.pop();
                core.push_back(u);

                for (int v : g->vertices_[u]->nbr) {
                    if (!mark.get(v)) {
                        mark.set(v);
                        bfs_q.push(v);
                    }
                }
            }

            printf("Found k-core with size: %lu\n", core.size());  // 打印 k-core 大小
            kcores.push_back(core);
        }
    }

    delete[] deg;
    return kcores;
}

vector<Graph*> VCC::global_cut(Graph* g, const int& k){

	int n = g->vertices_num_;

	int* deposit = new int[n];
	memset(deposit,0,sizeof(int)*n);

	BitMap pru(n);

	// Arry vertex_to_group maintains a match from a vertex to a group containing that vertex
	// The groups are generated as byproducts when computing the sparse certificate
	int* vertex_to_group = new int[n];
	// If the group id of a vertex is -1, that means there is no group containing that vertex
	memset(vertex_to_group,-1,sizeof(int)*n);

	// Each item in groups demonstrates a group and contains the ids of all vertices inside.
	vector<vector<int> > groups;

	// Compute the sparse certificate of the original graph
	Graph* sc = get_sparse_certificate(g,groups,vertex_to_group,k);

	int* group_deposit = new int[groups.size()];
	memset(group_deposit,0,sizeof(int)*groups.size());

	Graph* df_graph = create_directed_flow_graph(sc);


	int source_id = get_source_id(sc);

	sweep(source_id,pru,deposit,vertex_to_group,group_deposit,groups,g,k);

	// Sort the vertices by the distance to the source vertex
	vector<int> distance_order = create_distance_order(source_id,sc);

	// Test the connectivity between the source vertex and furthest vertex first
	for(int i = distance_order.size()-1; i > 0; --i){
		int target_id = distance_order[i];

		vector<Graph*> sides = local_cut(source_id,target_id,g,df_graph,k);

		// If sides is empty, that means the connectivity between source_id and target_id is not less than k
		if(sides.empty()){
			sweep(target_id,pru,deposit,vertex_to_group,group_deposit,groups,g,k);
			continue;
		}

		release_graph(sc);
		release_graph(df_graph);
		delete[] deposit;
		delete[] group_deposit;
		delete[] vertex_to_group;
		return sides;
	}

	delete[] deposit;
	delete[] group_deposit;
	

	// Test whether the connectivity between any two neighbors of source vertex is not less than k
	for(auto nit = sc->vertices_[source_id]->nbr.begin(); nit != sc->vertices_[source_id]->nbr.end(); ++nit){
		if(*nit == -1) continue;
		for(auto nnit = nit+1; nnit != sc->vertices_[source_id]->nbr.end(); ++nnit){
			if(*nnit == -1) continue;
			vector<int>& nit_nbr = g->vertices_[*nit]->nbr;
			if(find(nit_nbr.begin(), nit_nbr.end(), *nnit) != nit_nbr.end()) continue;

			// If two neighbors in the same side group, skip the connectivity test.
			if(vertex_to_group[*nit] != -1 && vertex_to_group[*nit] == vertex_to_group[*nnit]) continue;

			// Test the local connecitivty between *nit and *nnit.
			vector<Graph*> sides = local_cut(*nit,*nnit,g,df_graph,k);

			if(sides.empty()) continue;

			release_graph(sc);
			release_graph(df_graph);
			delete[] vertex_to_group;
			return sides;
		}

	}
	

	release_graph(sc);
	release_graph(df_graph);
	delete[] vertex_to_group;
	
	return vector<Graph*> ();
}

// Test the vertex connectivity between source_id and target_id
vector<Graph*> VCC::local_cut(const int& source_id, const int& target_id, Graph* g, Graph* df_graph, const int& k){
	
	int n = g->vertices_num_;

	int source_end = source_id + n;


	vector<F_vertex*> f_graph = create_flow_graph(df_graph);


	int f = flow(source_end,target_id,f_graph,k);


	vector<Graph*> sides;

	if(f < k){
		sides = cut_graph(source_id,target_id,f_graph,g);
	}

	for(F_vertex* fv : f_graph) delete fv;
	return sides;
}

// Partition the graph via vertex cut
vector<Graph*> VCC::cut_graph(const int& source_id, const int& target_id, vector<F_vertex*>& f_graph, Graph* g){

	int n = g->vertices_num_;

	// mark nodes during BFS, initialized by -1
	// a node in graph g will be represented as 2 nodes in directed graph
	// if BFS find 1 of 2 nodes of node u (here u is cut node), mark[u] -> 0;
	// if BFS find both 2 nodes of node u, mark[u] -> 1;
	int* visit = new int[n];
	memset(visit,-1,sizeof(int)*n);

	BitMap mark(n*2);

	queue<int> list;
	list.push(source_id);
	mark.set(source_id);

	while(!list.empty()){
		int id = list.front();
		list.pop();

		++ visit[id < n ? id : id-n];

		for(F_edge* fe : f_graph[id]->nbr_){
			if(fe->capacity_ == 0) continue;
			if(mark.get(fe->b_)) continue;
			list.push(fe->b_);
			mark.set(fe->b_);
		}
	}

	vector<int> side_0;
	vector<int> side_1;

	for(int i = 0; i < n; ++i){
		if(visit[i] != -1) side_0.push_back(i);
		if(visit[i] != 1) side_1.push_back(i);
	}


	vector<Graph*> sides;
	sides.push_back(create_new_graph(side_0,g));
	sides.push_back(create_new_graph(side_1,g));

	delete[] visit;

	return sides;
}

// Create a graph used for computing maximum flow
vector<F_vertex*> VCC::create_flow_graph(Graph* g){

	int n = g->vertices_num_;

	vector<F_vertex*> f_graph(n);

	for(int i = 0; i < n; ++i){
		f_graph[i] = new F_vertex(i);
	}

	for(int i = 0; i < n; ++i){
		Vertex* u_ptr = g->vertices_[i];
		F_vertex* fu_ptr = f_graph[i];
		for(int neighbor_id : u_ptr->nbr){
			// Create an edge from i to neighbor_id
			F_edge* fe = new F_edge(i,neighbor_id,1,NULL);
			fu_ptr->nbr_.push_back(fe);

			// Create a reverse edge from neighbor_id to i
			F_edge* fe_r = new F_edge(neighbor_id,i,0,fe);
			f_graph[neighbor_id]->nbr_.push_back(fe_r);

			// Link two edges
			fe->reverse_ = fe_r;
		}
	}

	return f_graph;
}

// Test whether the flow from source_id to target_id is less than k or not
// Return k if the exact flow is not less than k
// Return exact flow if less than k

int VCC::flow(const int& source_id, const int& target_id, vector<F_vertex*>& f_graph, const int& k){


	int* label = new int[f_graph.size()];
	memset(label, -1, sizeof(int) * f_graph.size());  // 显式初始化
	vector<int> label_count;

	int flow_val = 0;
	bool gap = false;

	while(label_graph(source_id,target_id,f_graph,label,label_count)){
		BitMap pru(f_graph.size());
		while(true){
			int res = find_augmenting_path(source_id,target_id,f_graph,label,label_count,pru);
			if(res == 1){
				++ flow_val;
				if(flow_val < k) continue;
			}else if(res == -1){
				gap = true;
			}
			break;
		}
		if(flow_val >= k || gap) break;
		pru.reset();
	}




	delete[] label;


	return flow_val;

}

int VCC::find_augmenting_path(const int& source_id, const int& target_id, vector<F_vertex*>& f_graph, int* label, vector<int>& label_count, BitMap& pru){

	if(pru.get(source_id)) return 0;

	// Find the target vertex
	if(source_id == target_id) return 1;

	int end_label = label[target_id];
	// Reach the same label of target vertex
	if(label[source_id] == end_label) return 0;

	int current_label = label[source_id];

	bool gap = true;
	for(F_edge* fe : f_graph[source_id]->nbr_){
		if(fe->capacity_ == 0) continue;

		int neighbor_id = fe->b_;

		if(label[neighbor_id] != current_label + 1) continue;

		if(gap) gap = false;

		int ap = find_augmenting_path(neighbor_id,target_id,f_graph,label,label_count,pru);

		if(ap == -1) return -1;
		else if(ap == 1){
			// Finding an angmenting path
			fe->capacity_ = 0;
			fe->reverse_->capacity_ = 1;

			return 1;
		}else{
			pru.set(neighbor_id);
		}
	}
	if(gap){
		label_count[current_label] --;
		if(label_count[current_label] == 0) return -1;
	}

	return 0;

}

// Use BFS to assign distance label for vertices
bool VCC::label_graph(const int& source_id, const int& target_id, vector<F_vertex*>& f_graph, int* label, vector<int>& label_count){
	
	int n = f_graph.size();

	for(int i = 0; i < n; ++i){
		if(label[i] != -1) label[i] = -1;
	}

	label_count.clear();

	int current_label = 0;
	int end_label = -1;

	BitMap mark(f_graph.size());

	queue<int> list;
	list.push(source_id);
	mark.set(source_id);

	while(! list.empty()){

		// if(current_label == end_label){
		// 	label[target_id] = end_label;
		// 	return true;
		// }

		int l_size = list.size();
		// Count the number of vertices in each label
		label_count.push_back(l_size);

		for( ;l_size > 0; -- l_size){
			int id = list.front();
			list.pop();

			label[id] = current_label;

			if(current_label == end_label) continue;

			for(F_edge* fe : f_graph[id]->nbr_){
				if(fe->capacity_ == 0) continue;
				int neighbor_id = fe->b_;
				if(mark.get(neighbor_id)) continue;
				mark.set(neighbor_id);

				if(neighbor_id == target_id) end_label = current_label+1;

				list.push(neighbor_id);
			}

		}

		 ++ current_label;
	}

	if(end_label == -1) return false;
	return true;

}


vector<int> VCC::create_distance_order(const int& source, Graph* g){

	vector<int> distance_order;

	queue<int> list;
	BitMap mark(g->vertices_num_);
	list.push(source);
	mark.set(source);

	while(!list.empty()){
		int id = list.front();
		list.pop();
		distance_order.push_back(id);

		for(int neighbor_id : g->vertices_[id]->nbr){
			if(neighbor_id == -1) continue;
			if(mark.get(neighbor_id)) continue;
			list.push(neighbor_id);
			mark.set(neighbor_id);
		}
	}

	return distance_order;
}


void VCC::sweep(const int& v, BitMap& pru, int* deposit, int* vertex_to_group, int* group_deposit, vector<vector<int> >& groups, Graph* g, const int& k){

	pru.set(v);

	// If v is k-connected with source vertex,
	// increase the dposit of each neighbor of v
	for(int neighbor_id : g->vertices_[v]->nbr){
		if(neighbor_id == -1 || pru.get(neighbor_id)) continue;
		deposit[neighbor_id] ++;
		// If the deposit of a vertex is not less than k, prune such vertex.
		if(deposit[neighbor_id] >= k){
			if(pru.get(neighbor_id)) continue;
			sweep(neighbor_id,pru,deposit,vertex_to_group,group_deposit,groups,g,k);
		}
	}

	if(vertex_to_group[v] != -1){
		group_deposit[vertex_to_group[v]] ++;
		if(group_deposit[vertex_to_group[v]] >= k){
			for(int vertex_id : groups[vertex_to_group[v]]){
				if(pru.get(vertex_id)) continue;

				sweep(vertex_id,pru,deposit,vertex_to_group,group_deposit,groups,g,k);
			}
		}
	}

}

int VCC::get_source_id(Graph* g){
	int min_deg = g->vertices_num_;
	int source_id = 0;
	for(Vertex* u : g->vertices_){
		if(u->degree < min_deg){
			source_id = u->id;
			min_deg = u->degree;
		}
	}

	return source_id;
}

Graph* VCC::create_directed_flow_graph(Graph* g){
	Graph* df_graph = new Graph;
	
	int n = g->vertices_num_;

	df_graph->vertices_num_ = n*2;
	
	df_graph->vertices_.resize(n*2);

	for(int i = 0; i < n; ++i){
		

		// replace a vertex by a directed edge:  a -> b
		Vertex* a = new Vertex();
		
		// for a directed graph, degree here is out degree
		a->degree = 1;
		a->id = i;
		
		df_graph->vertices_[i] = a;

		a->nbr.push_back(n+i);

		Vertex* b = new Vertex();
		
		b->degree = g->vertices_[i]->degree;
		b->id = n+i;
		df_graph->vertices_[n+i] = b;
		
		// connected b to other a
		for(int neighbor_id : g->vertices_[i]->nbr){
			if(neighbor_id == -1) continue;
			b->nbr.push_back(neighbor_id);
		}
	}

	return df_graph;
}

// The sparse certificate is computed by performing k times scan first search
// The BFS (breadth first search) is a special case of scan first search
// The details can be found in the paper: 
// "Scan-First Search and Sparse Certificates: An Improved Parallel Algorithm for k-Vertex Connectivity"
Graph* VCC::get_sparse_certificate(Graph* g, vector<vector<int> >& groups, int* vertex_to_group, const int& k){
	
	int n = g->vertices_num_;

	Graph* sc = new Graph();

	sc->vertices_num_ = n;
	sc->vertices_.resize(n);

	for(int i = 0; i < n; ++i){
		sc->vertices_[i] = new Vertex();
		sc->vertices_[i]->id = i;
		sc->vertices_[i]->degree = 0;
	}

	BitMap mark(n);
	
	for(int i = 0; i<k; ++i){

		// Select a vertex as the root to start scan first search
		for(int root = i; root < n; ++root){
			if(mark.get(root)) continue;

			if(i == k-1) scan_first_search(root,g,sc,mark,groups,vertex_to_group,k);
			else scan_first_search(root,g,sc,mark,groups,NULL,k);

		}

		if(i == k-1) continue;

		mark.reset();

	}

	return sc;
}


// BFS is a special case of scan first search and is used here.
void VCC::scan_first_search(const int& root, Graph* g, Graph* sparse_certificate, BitMap& mark, vector<vector<int> >& groups, int* vertex_to_group, const int& k){
	// vector<Node*>& gNodes = g->nodes;
	// vector<Node*>& scNodes = sc->nodes;
	
	queue<int> list;
	list.push(root);
	
	mark.set(root);

	vector<int> component;
	int new_group_id = groups.size();

	// bool hasSide = false;
	while(!list.empty()){
		int vertex_id = list.front();
		list.pop();

		// When the input parameter vertex_to_group is not NULL,
		// that means this is the last time of scan first search and 
		// we collect the side-groups here
		if(vertex_to_group){
			// Push the vertex into the group
			component.push_back(vertex_id);
			// Match the vertex id to the group id
			vertex_to_group[vertex_id] = new_group_id;
		}

		Vertex* g_vertex = g->vertices_[vertex_id];
		Vertex* sc_vertex = sparse_certificate->vertices_[vertex_id];

		
		for(int neighbor_id : g_vertex->nbr){
			if(neighbor_id == -1) continue;
			if(mark.get(neighbor_id)) continue;

			auto nit = find(sc_vertex->nbr.begin(),sc_vertex->nbr.end(),neighbor_id);
			if(nit != sc_vertex->nbr.end()) continue;

			// Add neighbor to the nighbor list
			sc_vertex->nbr.push_back(neighbor_id);
			sc_vertex->degree += 1;

			// Add vertex id to the neighbor list of neighbor
			sparse_certificate->vertices_[neighbor_id]->nbr.push_back(vertex_id);
			sparse_certificate->vertices_[neighbor_id]->degree += 1;

			mark.set(neighbor_id);
			list.push(neighbor_id);

		}

	}

	// We do not need to collect the information of group whose size is less than k
	if(vertex_to_group){
		if(component.size()>k){
			groups.push_back(component);
		}else{
			for(int id : component){
				vertex_to_group[id] = -1;
			}
		}
	}
}



// Given the vertices of a subgraph, create a new graph
// Allocate new id for each vertex in the new graph
Graph* VCC::create_new_graph(const vector<int>& component, Graph* old_graph){
	if (component.empty()) {
         fprintf(stderr, "create_new_graph ERROR: empty component\n");
        return nullptr;
    }
    for (int id : component) {
        if (id < 0 || id >= old_graph->vertices_num_) {
            fprintf(stderr, "create_new_graph ERROR: invalid vertex id=%d (max=%d)\n", 
                    id, old_graph->vertices_num_);
            return nullptr;
        }
    }
	// Create a map from old vertex id to new vertex id

	int* old_to_new = new int[old_graph->vertices_num_];
	memset(old_to_new,-1,sizeof(int)*old_graph->vertices_num_);

	Graph* new_graph = new Graph();

	// Assign the vertices number of new graph
	new_graph->vertices_num_ = component.size();

	for(int id : component){
		// Build the match from old id to new id
		old_to_new[id] = new_graph->vertices_.size();
		// Create a new vertex instance in the new graph 
		create_new_vertex(id,old_graph,new_graph);
	}

	for(int id : component){

		// Get the pointer of new vertex instance
		Vertex* new_vertex = new_graph->vertices_[old_to_new[id]];
		
		// The FOR loop below adds neighbors and updates the degree for each vertex in the new graph.
		for(int neighbor_id : old_graph->vertices_[id]->nbr){
			
			// Since the match from old id to new id has been created,
			// we can obtain the corresponding new id 
			int new_neighbor_id = old_to_new[neighbor_id];
			if(new_neighbor_id == -1) continue;

			new_vertex->degree += 1;
			new_vertex->nbr.push_back(new_neighbor_id);

		}
	}

	delete[] old_to_new;
	return new_graph;
}

vector<vector<int> > VCC::get_all_connected_components(Graph* g){
	vector<vector<int> > ccs;
	

	BitMap* mark = new BitMap(g->vertices_num_);
	

	for(int i = 0; i < g->vertices_num_; ++i){
		
		if(mark->get(i)) continue;
		
		vector<int> cc;

		// detect connected component via BFS
		queue<int> list;
		list.push(i);

		// any processed vertex will be assigned as true
		mark->set(i);

		int id;
		int nid = 0;

		while(!list.empty()){
			id = list.front();
			list.pop();

			cc.push_back(id);

			Vertex* u = g->vertices_[id];
			Vertex* n;
		
			for(auto it = u->nbr.begin(); it != u->nbr.end(); ++it){
				// This neighbor is delete for some reasons before
				// Normally, it is not used here.
				if(*it == -1) continue;

				// vertex has been processed
				if(mark->get(*it)) continue;
				mark->set(*it);

				list.push(*it);
			}
			
		}

		ccs.push_back(cc);

	}

	delete mark;
	
	return ccs;
}

// Graph* VCC::get_connected_component(int vertex_id, Graph* g, BitMap* mark){
	
// 	int n = g->vertices_num_;
// 	int* old_to_new = new int[n];
// 	memset(old_to_new,-1,sizeof(int)*n);

// 	// create a new subgraph
// 	Graph* cc = new Graph;


// 	// detect connected component via BFS
// 	queue<int> list;
// 	list.push(vertex_id);
// 	// any processed vertex will be assigned as true
// 	mark->set(vertex_id,1);

// 	int id;
// 	int nid = 0;

// 	while(!list.empty()){
// 		id = list.front();
// 		list.pop();

// 		Vertex* u = g->vertices_[id];

// 		Vertex* n;
// 		if(old_to_new[id] == -1) n = create_new_vertex(id,g,cc);
// 		else n = cc->vertices_[old_to_new[id]];
	
// 		for(auto it = u->nbr.begin(); it != u->nbr.end(); ++it){
// 			// this neighbor is delete for some reasons before
// 			if(*it == -1) continue;

// 			if(old_to_new[*it] == -1) n = create_new_vertex(id,g,cc);

// 			++ n->degree;

// 			// add neighbors
// 			n->nbr.push_back(old_to_new[*it]);

// 			// vertex has been processed
// 			if(mark->get(*it)) continue;
			
// 			mark->set(*it,1);

// 			list.push(*it);
// 		}
		
// 	}

// 	delete[] old_to_new;
// 	return cc;
// }


Vertex* VCC::create_new_vertex(const int& old_id, Graph* old_graph, Graph* new_graph){

	Vertex* n = new Vertex;

	n->id = new_graph->vertices_.size();
	n->degree = 0;

	// map the id in the new graph to the original vertex id
	new_graph->id_new_to_old_.push_back(old_graph->id_new_to_old_[old_id]);
	
	new_graph->vertices_.push_back(n);
	
	return n;
}


int VCC::release_graph(Graph* g){
	// for(Vertex* u : g->vertices_) delete u;
	delete g;
	return 1;
}


