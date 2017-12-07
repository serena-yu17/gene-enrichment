// gene enrichment dll.cpp : Defines the exported functions for the DLL application.
//
#include "main.h"

using namespace std;

void bfs(int nthread, unordered_set<pair<int, int>>* chains, mutex* chains_mutex, int node1, int node2, unordered_map<int, unordered_set<int>> &forward)
{
	unordered_set<int> used = { node1 };
	vector<vector<Node>> levels;
	levels.push_back(vector<Node>(1, Node(node1, NULL)));
	vector<Node> targets;
	while (levels.back().size() && targets.size() == 0)
	{
		levels.push_back(vector<Node>());
		vector<Node>* top = &(levels[levels.size() - 2]);
		for (unsigned i = 0; i < top->size(); i++)
		{
			Node* elem = &((*top)[i]);
			auto ids = &(forward[elem->val]);
			for (int num : *ids)
				if (used.find(num) == used.end())
				{
					Node newNode = Node(num, elem);
					levels.back().push_back(newNode);
					used.insert(num);
					if (num == node2)
						targets.push_back(newNode);
				}
		}
	}
	Node* pnode;
	for (auto tg : targets)
	{
		pnode = &tg;
		while (pnode && pnode->next)
		{
			lock_guard<std::mutex> lock(*chains_mutex);
			chains->insert(make_pair(pnode->next->val, pnode->val));
			pnode = pnode->next;
		}
	}
}

//Use Fruchterman Reingold force directed layout to assign coordinates
unordered_map<int, pair<double, double>> buildGraph(unordered_set<pair<int, int>> chains, int n)
{
	Graph graph;
	using Vertex = boost::graph_traits<Graph>::vertex_descriptor;
	unordered_map<int, Vertex> idVertex;
	unordered_map<Vertex, int> vertexID;
	for (auto& pair : chains)
	{
		Vertex v1, v2;
		if (idVertex.find(pair.first) != idVertex.end())
			v1 = idVertex[pair.first];
		else
		{
			v1 = boost::add_vertex(graph);
			vertexID[v1] = pair.first;
			idVertex[pair.first] = v1;
		}
		if (idVertex.find(pair.second) != idVertex.end())
			v2 = idVertex[pair.second];
		else
		{
			v2 = boost::add_vertex(graph);
			vertexID[v2] = pair.second;
			idVertex[pair.second] = v2;
		}
		boost::add_edge(v1, v2, graph);
	}
	using Topology = boost::square_topology<boost::mt19937>;
	using Position = Topology::point_type;
	vector<Position> positions(boost::num_vertices(graph));
	Topology topology;
	boost::random_graph_layout(
		graph,
		boost::make_iterator_property_map(positions.begin(), boost::identity_property_map{}),
		topology
	);
	boost::fruchterman_reingold_force_directed_layout(
		graph,
		boost::make_iterator_property_map(positions.begin(), boost::identity_property_map{}),
		topology,
		boost::attractive_force([](Graph::edge_descriptor, double k, double d, Graph const&) { return (d*d) / k; })
	);
	auto index_map = boost::get(boost::vertex_index, graph);
	auto begin = positions.begin();
	boost::iterator_property_map<vector<Position>::iterator, decltype(index_map)> position_map(positions.begin(), index_map);
	std::unordered_map<int, pair<double, double>> result;
	BGL_FORALL_VERTICES(v, graph, Graph)
	{
		Position pos = position_map[v];
		result[vertexID[v]] = std::make_pair(pos[0], pos[1]);
	}
	return result;
}

//Build a graph based on coordinates. vertex, coordinates, vcount, edge, ecount are output values passed by reference. 
extern "C" __declspec(dllexport) void enrichGenes(int32_t** vertex, double** coordinates, int* vcount, int32_t** edge, int* ecount, const int32_t data[], int nData, const int32_t genelist[], int n)
{
	unordered_map<int, unordered_set<int>> forward;
	unordered_set<pair<int, int>> chains;
	mutex chains_mutex;
	for (int i = 0; i < nData; i += 2)
	{
		if (forward.find(data[i]) == forward.end())
			forward[data[i]] = unordered_set<int>();
		forward[data[i]].insert(data[i + 1]);
	}
	ctpl::thread_pool thrPool(thread::hardware_concurrency());
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			if (i != j)
				thrPool.push(bfs, &chains, &chains_mutex, genelist[i], genelist[j], forward);
	//bfs(0, &chains, &chains_mutex, genelist[i], genelist[j], forward);
	thrPool.~thread_pool();
	unordered_map<int, pair<double, double>> formatedGraph = buildGraph(chains, n);
	/*cout << "geneid\tx\ty\n";
	for (auto& item : formatedGraph)
	cout << item.id << '\t' << item.x << '\t' << item.y << '\n';*/
	if (formatedGraph.size())
	{
		double minx = 100, miny = 100;
		double distance = 100;
		for (auto& key : formatedGraph)
		{
			if (key.second.first < minx)
				minx = key.second.first;
			if (key.second.second < miny)
				miny = key.second.second;
		}
		for (auto& item : chains)
		{
			pair<double, double>* coord1 = &(formatedGraph[item.first]);
			pair<double, double>* coord2 = &(formatedGraph[item.second]);
			double distanceLocal = sqrt((coord1->first - coord2->first) * (coord1->first - coord2->first) + (coord1->second - coord2->second) * (coord1->second - coord2->second));
			if (distanceLocal < distance)
				distance = distanceLocal;
		}
		*vcount = formatedGraph.size();
		*ecount = chains.size();
		*vertex = new int32_t[*vcount];
		*coordinates = new double[(*vcount) * 2];
		int i = 0, j = 0;
		for (auto& key : formatedGraph)
		{
			(*vertex)[i++] = key.first;
			(*coordinates)[j++] = (key.second.first - minx) / distance + 1;
			(*coordinates)[j++] = (key.second.second - miny) / distance + 1;
		}
		*edge = new int32_t[(*ecount) * 2];
		i = 0;
		for (auto& item : chains)
		{
			(*edge)[i++] = item.first;
			(*edge)[i++] = item.second;
		}
	}
	cout << (*edge)[0] << ", " << (*edge)[1] << '\n';
}

//clean up arrays allocated using C++ "new"
extern "C" __declspec(dllexport) void freeArray(void* arr)
{
	delete[] arr;
}




