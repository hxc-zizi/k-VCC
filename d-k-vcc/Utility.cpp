#include "Utility.h"
#include <unordered_map> 

FILE *open_file(const char *file_name, const char *mode) {
	FILE *f = fopen(file_name, mode);
	if(f == NULL) {
		printf("Can not open file: %s\n", file_name);
		exit(1);
	}

	return f;
}

Utility::Utility(){

}

Utility::~Utility(){

}

void Utility::format_graph(string dest, string src){
	format_graph(dest, src, 100, false);
}

void Utility::format_graph(string dest, string src, int percent, bool is_sample_edge){
	printf("Create binary graph from .txt file: %s\n",src.c_str());
	percent = percent/10-1;
	
	int mem_edges = 100000000;
	dir_ = dest;

	FILE* fp = open_file(src.c_str(),"r");


	Arc* edges = new Arc[mem_edges];


	vertex_map_ = new int[MAX_ID];
	memset(vertex_map_,-1,sizeof(int)*MAX_ID);

	// save every line read from txtFile
	char line[100];

	int u,v;

	unsigned long es = 0;

	int size = 0,num = 0,tmp_file = 0;
	printf("Separating: \n");
	while(fgets(line,100,fp)){
		
		if( line[0] < '0' || line[0] > '9' )
			continue;

		sscanf(line,"%d\t%d",&u,&v);

		if(u == v) continue;
		++es;
		if(is_sample_edge){
			if(es%10>percent) continue;
		}else{
			if(u%10>percent || v%10>percent) continue;
		}

		u = get_vertex_id(u,num);
		v = get_vertex_id(v,num);


		edges[size].a = u;
		edges[size].b = v;
		++size;

		edges[size].a = v;
		edges[size].b = u;
		++size;

		//检验
		printf("Loaded edge: %d <-> %d\n", u, v);  // 打印加载的边
		 
		if(size>=mem_edges){
			printf("load %d edges\n",size);
			save_tmp_edges(edges,size,tmp_file);
			size = 0;
			++tmp_file;
		}

	}
	printf("load %d edges\n",size);
	save_tmp_edges(edges,size,tmp_file);

	int* new2old = new int[num];
	int cur = 0;
	int i = 0;
	while(cur < num){
		if(vertex_map_[i] == -1){
			++i;
			continue;
		}
		new2old[vertex_map_[i]] = i;
		++cur;
		++i;
	}

	// write match file
	// printf("write match file\n");
	string mts = dest+"match.st";
	FILE* mch = fopen(mts.c_str(),"wb");
	fwrite(new2old,sizeof(int),num,mch);
	fclose(mch);
	
	delete[] new2old;
	delete[] vertex_map_;
	delete[] edges;

	fclose(fp);

	merge(tmp_file+1);
}

bool Utility::edge_compare(const Arc &e1, const Arc &e2){
	if(e1.a < e2.a){
		return true;
	}
	if( e1.a > e2.a ){
		return false;
	}
	return e1.b < e2.b;
}

// sort edges and save
void Utility::save_tmp_edges(Arc* edges, int size, int tmp_file){
	printf("sort edges...\n");
	sort(edges,edges+size,edge_compare);

	string destDir = dir_+"sort_edge_tmp";

	if(access(destDir.c_str(), 0) == -1){
		// use for LINUX
		#ifdef _LINUX_
		int flag = mkdir(destDir.c_str(),0777);
		#else
		int flag=mkdir(destDir.c_str());
		#endif
		if(flag == 0) printf("Directory \"%s\" is created\n",destDir.c_str() );
		else{
			printf("Create directory failed\n");
			exit(1);
		}
	}


	char file_name[200];
	sprintf(file_name,"%ssort_edge_tmp/edges_tmp_%d",dir_.c_str(),tmp_file);
	printf("Creating tmp_file_%d: %s\n",tmp_file,file_name);

	FILE* fo = fopen(file_name,"wb");
	for (int i = 0; i < size; ++i){
		fwrite( edges+i, sizeof(Arc), 1, fo );
		// printf("edge[%d,%d]\n",edges[i].a,edges[i].b );
	}
	printf("------\n\n");
	fclose(fo);
	
}

// merge_sort edges from all edges_tmp files
void Utility::merge(int size){
	FILE** frl = new FILE*[size];
	Arc* es = new Arc[size];
	string idx_path = dir_+"graph.idx";
	string dat_path = dir_+"graph.dat";
	FILE* f_idx = fopen(idx_path.c_str(),"wb");
	FILE* f_dat = fopen(dat_path.c_str(),"wb");

	long pos;

	for (int i = 0; i < size ; ++i){
		char file_name[200];
		sprintf(file_name,"%ssort_edge_tmp/edges_tmp_%d",dir_.c_str(),i);
		
		frl[i] = fopen(file_name,"rb");
		// get first edge of all edges_tmp files
		fread(&es[i],sizeof(Arc),1,frl[i]);
	}

	unsigned int edgeNum = 0;

	int minIndex,previousA = -1,previousB = -1;

	int i = 0;

	int max_degree = 0;

	int degree = -1;
	// printf("start merge\n");
	int x=0;
	int f = 0;
	while(is_merge_finished(es,size)){
		minIndex = get_min_edge(es,size);
		int u = es[minIndex].a;
		int v = es[minIndex].b;
		
		if(u != previousA){
			// u != previousA demonstrates that all previousA's neighbors have been writen
			if (u%100000 == 0){
				printf("[%d]\n",u );
			}
			
			if(degree != -1){

				// write the vertex degree in .idx file
				fwrite(&degree,sizeof(int),1,f_idx);
				edgeNum += degree;
				max_degree = degree>max_degree?degree:max_degree;
				// printf("max degree: %d\n",max_degree );
			}

			// write the vertex beginning position in .dat file to .idx file
			pos = ftell(f_dat);
			fwrite(&pos,sizeof(long),1,f_idx);

			degree = 1;

			fwrite(&v,sizeof(int),1,f_dat);


		}else if(v != previousB){

			fwrite(&v,sizeof(int),1,f_dat);
			++degree;

		}
		
		// if u==previousA & v==previousB, ignore edge(u,v) cause it is same as previous one.

		previousA = u;
		previousB = v;

		// replace es[minIndex] by picking up the first edge from file edges_tmp_minIndex
		if(!fread(&es[minIndex],sizeof(Arc),1,frl[minIndex])){
			es[minIndex].a = MAX_ID;
		}
		

	}

	fwrite(&degree,sizeof(int),1,f_idx);
	edgeNum += degree;
	edgeNum /= 2;
	max_degree = degree>max_degree?degree:max_degree;

	// write the vertex num and max degree
	int vertexNum = previousA+1;
	printf("vertex num: %d\n",vertexNum);
	printf("edge num: %d\n",edgeNum);

	string infoPath = dir_+"graph.info";
	FILE* fInfo = fopen(infoPath.c_str(),"wb");
	fwrite(&vertexNum,sizeof(int),1,fInfo);
	fwrite(&edgeNum,sizeof(unsigned int),1,fInfo);
	fwrite(&max_degree,sizeof(int),1,fInfo);
	fclose(fInfo);

	for (int i = 0; i < size ; ++i){
		fclose(frl[i]);
	}

	fclose(f_idx);
	fclose(f_dat);

	delete[] frl;
	delete[] es;
}

bool Utility::is_merge_finished(Arc* es, int size){
	
	for (int i = 0; i < size; ++i){
		if(es[i].a!=MAX_ID){
			return true;
		}
	}

	return false;
}

// get minimum edge from edge list
int Utility::get_min_edge(Arc* es,int size){
	int min = 0;
	for (int i = 1; i < size; ++i){
		if(es[i].a < es[min].a){
			min = i;
		}else if(es[i].a > es[min].a){
			continue;
		}else if(es[i].b < es[min].b){
			min = i;
		}
	}
	
	return min;
}

// get final vertexID from map array
int Utility::get_vertex_id(int u,int &num){
	if(vertex_map_[u]<0){
		vertex_map_[u] = num++;
	}
	return vertex_map_[u];
}


Graph* Utility::create_graph(string txt_file_name){
	
	// input example "/DW/workspace/dataset/stanford.txt"

	string::size_type pos = txt_file_name.rfind("/");
	
	// get file direcotry
	// e.g. "/DW/workspace/dataset/"
	string file_dir = txt_file_name.substr(0,pos+1);

	// get graph_name
	// e.g. "stanford"
	string graph_name = txt_file_name.substr(pos+1, txt_file_name.rfind(".")-pos-1);
	
	// create binary file directory
	// e.g. "/DW/workspace/dataset/stanford"
	string binary_dir = file_dir + graph_name;
	

	// if binary file directory does not exist, create it.
	if(access(binary_dir.c_str(), 0) == -1){
		// use for LINUX
		#ifdef _LINUX_
		int flag = mkdir(binary_dir.c_str(),0777);
		#else
		int flag=mkdir(binary_dir.c_str());
		#endif
		if(flag == 0) printf("Binary files directory \"%s\" is created\n",binary_dir.c_str() );
		else{
			printf("Create directory failed\n");
			exit(1);
		}

		// create binary file from original .txt file
		format_graph(binary_dir+"/", txt_file_name);
		
	}else{
		printf("Binary file directory exists!\n");
	}
	

	return load_from_binary_graph(binary_dir+"/");
}


Graph* Utility::load_from_binary_graph(string binary_file_dir){

	Graph* g = new Graph();


	// assign basic information to graph instance
	string info_file = binary_file_dir+"graph.info";
	FILE* f_info = fopen(info_file.c_str(),"rb");
	fread(&g->vertices_num_,sizeof(int),1,f_info);
	fread(&g->edges_num_,sizeof(int),1,f_info);
	fread(&g->max_deg_,sizeof(int),1,f_info);
	fclose(f_info);
	printf("number of vertices: %d, number of edges: %d, max degree: %d\n", g->vertices_num_, g->edges_num_, g->max_deg_ );
	

	string idx_path = binary_file_dir+"graph.idx";
	string dat_path = binary_file_dir+"graph.dat";
	FILE* f_idx = fopen(idx_path.c_str(),"rb");
	FILE* f_dat = fopen(dat_path.c_str(),"rb");

	// Initialize degree and neighbors for every vertex
	int degree;
	long pos;
	for(int i = 0;i < g->vertices_num_;++i){
		fread(&pos,sizeof(long),1,f_idx);
		fread(&degree,sizeof(int),1,f_idx);

		// construct node
		Vertex* n = new Vertex;
		n->id = i;
		n->degree = degree;
		int nbr_id;

		// read node's neighbors
		for(int j = 0;j < degree;++j){
			fread(&nbr_id,sizeof(int),1,f_dat);
			n->nbr.push_back(nbr_id);
		}

		// push node's pointer into graph instance
		g->vertices_.push_back(n);
		
//		//检验
//		printf("Vertex %d has degree %d and neighbors: ", i, degree);  // 打印顶点度数
//        for (int nbr : n->nbr) {
//            printf("%d ", nbr);
//        }
//        printf("\n"); 
	}


	fclose(f_dat);
	fclose(f_idx);

	// construct the map from new id to original id
	string match_path = binary_file_dir+"match.st";
	FILE* f_match = fopen(match_path.c_str(),"rb");
	int old_id;
	for(int i = 0; i < g->vertices_num_; ++i){
		fread(&old_id,sizeof(int),1,f_match);
		g->id_new_to_old_.push_back(old_id);
	}
	fclose(f_match);


	return g;
}

int Utility::release_graph(Graph* g){
	delete g;
	return 1;
}
// 分割图
std::vector<std::vector<int>> Utility::split_graph(Graph* g, int num_parts) {
    std::vector<std::vector<int>> subgraphs(num_parts);
    int vertices_per_part = g->vertices_num_ / num_parts;
    for (int i = 0; i < g->vertices_num_; ++i) {
        int part = i % num_parts;
        subgraphs[part].push_back(i);
           // 检验 
//	printf("Vertex %d assigned to subgraph %d\n", i, part);
    }

	
    return subgraphs;
}

// 创建局部图
Graph* Utility::create_local_graph(const std::vector<int>& vertices, Graph* old_graph) {
    Graph* local_graph = new Graph();
    local_graph->vertices_num_ = vertices.size();

    // 全局ID到局部ID的映射
    std::unordered_map<int, int> global_to_local;
    for (int i = 0; i < vertices.size(); ++i) {
        int global_id = vertices[i];
        global_to_local[global_id] = i;
    }

    // 复制顶点及其邻居
    for (int global_id : vertices) {
        Vertex* old_vertex = old_graph->vertices_[global_id];
        if (!old_vertex) {
            printf("Error: Vertex %d is null\n", global_id);
            exit(1);
        }

        Vertex* new_vertex = new Vertex();
        new_vertex->id = global_to_local[global_id];
        new_vertex->degree = old_vertex->degree;  // 确保度数正确

        // 转换邻居的全局ID为局部ID
        for (int neighbor_global_id : old_vertex->nbr) {
		    auto it = global_to_local.find(neighbor_global_id);
		    if (it != global_to_local.end()) {
		        new_vertex->nbr.push_back(it->second);
		    } else {
		        // 如果邻居不在当前子图中，可以选择忽略或记录警告
		        fprintf(stderr, "WARNING: neighbor %d not in local graph\n", neighbor_global_id);
		    }
		}

        local_graph->vertices_.push_back(new_vertex);
    }

    return local_graph;
}
