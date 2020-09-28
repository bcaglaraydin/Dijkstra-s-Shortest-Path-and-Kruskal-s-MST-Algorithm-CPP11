/*
---DIJKSTRA'S SHORTEST PATH ALGORITHM---
---KRUSKAL'S MINIMUM SPANNING TREE ALGORITHM---
Author : Berdan Çaðlar AYDIN
School : Istanbul Technical University
Date : 27/09/2020
*/

#include <iostream>
#include <vector>
#include <ctime>
#include <iomanip>
#include <limits> 
#include <string>
#include <fstream>
#include <queue>
#include <algorithm>

using namespace std;

//edge class
class edge {
public:
	edge(int iv1 = 0, int iv2 = 0, int icost = 0) {
		v1 = iv1;
		v2 = iv2;
		cost = icost;
	}
	int get_cost() const { return cost; };
	int get_v1() const { return v1; };
	int get_v2() const { return v2; };
		
private:
	int v1;
	int v2;
	int cost;
};

//using struct in order to return two vectors from function dijkstra
struct road {
	vector<double> distances;
	vector<int> shortest_path;
};

//Graph class
class Graph {
public:
	Graph(int node_number, double edge_density, double distrange_min, double distrange_max);
	Graph(string fname);
	~Graph();
	int num_V() { return vertex_num; };
	int num_E() { return edge_num; };
	bool is_adjacent(int x, int y);
	void list_neighbors(int x);
	void add_edge(int x, int y, double dist_min, double dist_max);
	void delete_edge(int x, int y);
	double get_edge_value(int x, int y) { return G[x][y]; };
	void set_edge_value(int x, int y, double a) { G[x][y] = a; };
	void print_graph();
	void path(int i, int j);
	double path_size(int i, int j);
	struct road dijkstra(int src);
	int min_dist(vector<double> d, vector<bool> s);
	void print_path(vector<int> p,int j);
	vector<edge> kruskal();
	int find_parent(int v, vector<int> parent);
private:
	double** G;
	int vertex_num = 0;
	int edge_num = 0;
	vector<edge>edge_vector;
};

//graph constructor with filename
Graph::Graph(string fname) {
	ifstream fptr(fname);
	fptr >> vertex_num;

	int num;
	int v1;
	int v2;
	int costs;

	while (!fptr.eof()) {
		for (int x = 0; x < 3; x++) {
			fptr >> num;
			switch (x) {
			case 0: v1 = num; break;
			case 1: v2 = num; break;
			case 2: costs = num; break;
			}
		}
		edge_vector.push_back(edge(v1, v2, costs));
	}


	G = new double* [vertex_num];

	for (int i = 0; i < vertex_num; i++) {
		G[i] = new double[vertex_num];
	}

	for (int i = 0; i < vertex_num; i++) {
		for (int j = 0; j < vertex_num; j++) {
			G[i][j] = 0;
		}
	}

	for (auto i = edge_vector.begin(); i != edge_vector.end(); i++) {
		int new_v1 = (*i).get_v1();
		int new_v2 = (*i).get_v2();
		int new_cost = (*i).get_cost();
		if (new_cost != 0) {
			G[new_v1][new_v2] = new_cost;
		}
	}

}

//Graph constructor with parameters
Graph::Graph(int node_number, double edge_density, double dist_min, double dist_max) {
	srand(time(0));
	//I used a adjacency matrix to represent the graph
	G = new double* [node_number];
	vertex_num = node_number;

	for (int i = 0; i < node_number; i++) {
		G[i] = new double[node_number];
	}

	for (int i = 0; i < node_number; i++) {
		for (int j = 0; j < node_number; j++) {
			G[i][j] = 0;
		}
	}


	for (int i = 0; i < node_number; i++) {
		for (int j = 0; j < node_number; j++) {

			//calculating a random edge weight
			int num = (rand() % 100);
			bool is_edge = num < (edge_density * 100);

			if (is_edge) {
				add_edge(i, j, dist_min, dist_max);
				edge_num++;
			}
			else {
				G[i][i] = 0;
			}

		}
	}

}

//Graph deconstructor
Graph::~Graph() {
	for (int i = 0; i < vertex_num; i++) {
		delete[] G[i];
	}
	delete[] G;
}

//in order to implement minHeap with edge costs
class compare_edge{
public:
	int operator() (const edge& e1, const edge& e2)
	{
		return e1.get_cost() > e2.get_cost();
	}
};

//finding the most important parent for given node
int Graph::find_parent(int v, vector<int> parent) {
	if (parent[v] == v) {
		return v;
	}
	return find_parent(parent[v], parent);
}

//implementing kruskal's algorithm
vector<edge> Graph::kruskal() {
	vector<edge> selected;

	priority_queue <edge, vector<edge>, compare_edge > pq;
	for (auto i = edge_vector.begin(); i != edge_vector.end(); i++) {
		pq.push((*i));
	}

	vector<int> parent(vertex_num);
	for (int i = 0; i < vertex_num; i++) {
		parent[i] = i;
	}

	int v_count = 0;
	while (!pq.empty() && v_count != (vertex_num - 1)) {
		edge e = pq.top();

		int p1 = find_parent(e.get_v1(), parent);
		int p2 = find_parent(e.get_v2(), parent);

		if (p1 != p2) {
			selected.push_back(e);
			parent[p1] = p2;
			v_count++;
			pq.pop();
		}
		else {
			pq.pop();
		}
	}
	return selected;
}

//print_path and path functions to print the shortest path
void Graph::path(int i, int j) {
	vector<int> p = dijkstra(i).shortest_path;
	print_path(p,j);
}

void Graph::print_path(vector<int> p, int j) {
	if (p[j] == -1) {
		return;
	}
	print_path(p,p[j]);
	cout << "->" << j;
}

//function to return the shortest path size
double Graph::path_size(int i, int j) {
	return dijkstra(i).distances[j];
}

//function to find the vertex with minimum distance
int Graph::min_dist(vector<double> d, vector<bool> s){
	int min = INT_MAX, min_index = -1;

	for (int i = 0; i < vertex_num; i++)
		if (s[i] == false && d[i] <= min)
			min = d[i], min_index = i;

	return min_index;
}

//implementing dijksta's algorithm
struct road Graph::dijkstra(int src) {
	vector<double> dist;
	vector<bool>set;
	vector<int>path(vertex_num);

	path[src] = -1;
	for (int i = 0; i < vertex_num; i++) {
		dist.push_back(INT_MAX);
		set.push_back(false);
	}

	dist[src] = 0;
	for (int i = 0; i < vertex_num - 1; i++) {
		int v = min_dist(dist, set);
		set[v] = true;
		for (int i = 0 ; i < vertex_num; i++) {
			if (!set[i] && G[v][i] && dist[v] != INT_MAX && dist[v] + G[v][i] < dist[i]) {
				path[i]= v;
				dist[i] = dist[v] + G[v][i];
			}
		}
	}

	struct road out;
	out.distances = dist;
	out.shortest_path = path;
	return out;
}

//bool function to check if two vertexes are adjacent
bool Graph::is_adjacent(int x, int y) {
	if (G[x][y] != 0) {
		return true;
	}
	else {
		return false;
	}
}

//listing the vertexes which have an edge with given vertex
void Graph::list_neighbors(int x) {
	cout << "The vertex " << x <<  " is adjacent with: ";
	for (int i = 0; i < vertex_num; i++) {
		if (is_adjacent(x,i)) {
			cout << i << " ";
		}
	}
	cout << endl;
}

//function to print the matrix
void Graph::print_graph() {
	cout << "  ";
	for (int i = 0; i < vertex_num; i++)
	{
		cout << i << "    ";
	}
	cout << endl;
	cout << endl;
	for (int i = 0; i < vertex_num; i++)
	{
		cout << i << " ";
		for (int j = 0; j < vertex_num; j++)
		{
			cout << fixed << setprecision(2) << G[i][j] << " ";
		}
		cout << endl;
		cout << endl;
	}
}

//function to add edge between given vertex
void Graph::add_edge(int x, int y, double dist_min, double dist_max) {
	if (G[x][y] == 0) {
		double edge_distance = (dist_max - dist_min) * (static_cast<double>(rand()) / static_cast<double>(RAND_MAX)) + dist_min;
		G[x][y] = edge_distance;
		G[y][x] = edge_distance;
	}
}

//function to delete the edge between given edges
void Graph::delete_edge(int x, int y) {
	if (G[x][y] != 0) {
		G[x][y] = 0;
	}
}

//reading graph from file
void graph_from_file() {
	string fname;
	cout << "Enter file name: ";
	cin >> fname;
	Graph myG(fname);
	myG.print_graph();
}

int main() {
	int sel,v;
	string fname;
	double density, min, max;
	double mst = 0;

	cout << "Construct a graph:" << endl;
	cout << "1. From a file" << endl;
	cout << "2. Random (with inputs)" << endl;

	cin >> sel;
	switch(sel){
	case 1:
		cout << "Enter file name: ";
		cin >> fname;
		break;

	case 2:
		cout << "Vertex number: ";
		cin >> v;
		cout << endl;

		cout << "Edge density (0 - 1): ";
		cin >> density;
		cout << endl;

		cout << "Path length min: ";
		cin >> min;
		cout << endl;

		cout << "Path length max: ";
		cin >> max;
		cout << endl;
		break;

	default:
		break;
	}

	Graph* myG(sel == 1 ? new Graph(fname) : new Graph(v, density, min, max));
	myG->print_graph();
	vector<edge> kruskal_sol;

	cout << "Algorithm: " << endl;
	cout << "1. DIJKSTRA'S SHORTEST PATH ALGORITHM" << endl;
	cout << "2. KRUSKAL'S MINIMUM SPANNING TREE ALGORITHM" << endl;
	cin >> sel;
	switch (sel) {
	case 1:
		int vertex;

		cout << "Starting vertex: ";
		cin >> vertex; cout << endl;

		for (int i = 0; i < myG->num_V() ; i++) {
			//showing the mininmum distances
			cout << "Distance between " << vertex << " - " << i << " is: " << myG->path_size(vertex, i) << endl;

			//showing the shortest path
			cout << "Shortest path: ";
			myG->path(vertex, i);

			cout << endl;
			cout << endl;
		}
		break;

	case 2:
		cout << "Minimum spanning tree:" << endl;
		cout << "(V1 " << "V2 " << "Weight)" << endl;

		kruskal_sol = myG->kruskal();

		for (auto i = kruskal_sol.begin(); i != kruskal_sol.end(); i++) {
			//printing out the edges
			cout << "(" << (*i).get_v1() << " " << (*i).get_v2() << " " << (*i).get_cost() << ")" << endl;
			//calculating the total cost of the minimun spanning tree
			mst = mst + (*i).get_cost();
		}
		cout << "Total cost of MST: " << mst << endl;
		break;

	default:
		break;
	}

	return 0;
}