// gene enrichment dll.cpp : Defines the exported functions for the DLL application.
//
#include "main.h"

using namespace std;

void dijkstraSearch(uint32_t nthread, unordered_set<pair<uint32_t, uint32_t>>* chains, mutex* chains_mutex, uint32_t node1, uint32_t node2, WGraph const& graph)
{
	vector<size_t> predecessors(boost::num_vertices(graph));
	vector<float> distances(boost::num_vertices(graph));
	auto indexmap = boost::get(boost::vertex_index, graph);
	auto weightmap = boost::get(boost::edge_weight, graph);
	boost::dijkstra_shortest_paths(graph,
		node1,
		&predecessors[0],
		&distances[0],
		weightmap,
		indexmap,
		less<float>(),
		boost::closed_plus<float>(),
		(std::numeric_limits<uint32_t>::max)(),
		0,
		boost::default_dijkstra_visitor()
	);
	if (distances[node2] != inf)
	{
		lock_guard<std::mutex> lock(*chains_mutex);
		for (size_t vtx = node2; vtx != node1; vtx = predecessors[vtx])
			chains->insert(make_pair(predecessors[vtx], vtx));
	}
}

//Use Fruchterman Reingold force directed layout to assign coordinates
unordered_map<uint32_t, pair<double, double>> buildGraph(unordered_set<pair<uint32_t, uint32_t>> chains, uint32_t const genelist[], uint32_t n)
{
	UGraph graph;
	unordered_map<uint32_t, size_t> idVertex;
	unordered_map<size_t, uint32_t> vertexID;
	for (auto& pair : chains)
	{
		size_t v1, v2;
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
	for (uint32_t i = 0; i < n; i++)
		if (idVertex.find(genelist[i]) == idVertex.end())
		{
			size_t vtx = boost::add_vertex(graph);
			vertexID[vtx] = genelist[i];
			idVertex[genelist[i]] = vtx;
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
		boost::attractive_force([](UGraph::edge_descriptor, double k, double d, UGraph const&) { return (d*d) / k; })
	);
	auto index_map = boost::get(boost::vertex_index, graph);
	auto begin = positions.begin();
	boost::iterator_property_map<vector<Position>::iterator, decltype(index_map)> position_map(positions.begin(), index_map);
	std::unordered_map<uint32_t, pair<double, double>> result;
	BGL_FORALL_VERTICES(v, graph, UGraph)
	{
		Position pos = position_map[v];
		result[vertexID[v]] = std::make_pair(pos[0], pos[1]);
	}
	return result;
}

//Build a graph based on coordinates. vertex, coordinates, vcount, edge, ecount are output values passed by reference. 
extern "C" __declspec(dllexport) void enrichGenes(uint32_t** vertex, double** coordinates, uint32_t* vcount, uint32_t** edge, uint32_t* ecount, uint32_t const data[], uint32_t nData, uint32_t const genelist[], uint32_t n)
{
	uint32_t maxGene = 0;
	uint32_t maxWeight = 0;
	for (uint32_t i = 0; i < nData; i++)
	{
		uint32_t j = i * 3;
		if (data[j] > maxGene)
			maxGene = data[j];
		if (data[j + 1] > maxGene)
			maxGene = data[j + 1];
		if (data[j + 2] > maxWeight)
			maxWeight = data[j + 2];
	}
	WGraph graph;
	for (uint32_t i = 0; i < maxGene + 1; i++)
		boost::add_vertex(graph);
	float fMaxWt = (float)log2(maxWeight);
	for (uint32_t i = 0; i < nData; i++)
	{
		uint32_t j = i * 3;
		boost::add_edge((Vertex)data[j], (Vertex)data[j + 1], fMaxWt - (float)log2(data[j + 2]) + 1, graph);
	}
	unordered_set<pair<uint32_t, uint32_t>> chains;
	mutex chains_mutex;
	//ctpl::thread_pool thrPool(thread::hardware_concurrency());
	for (uint32_t i = 0; i < n; i++)
		for (uint32_t j = 0; j < n; j++)
			if (i != j)
				//thrPool.push(dijkstraSearch, &chains, &chains_mutex, genelist[i], genelist[j], graph);
				dijkstraSearch(0, &chains, &chains_mutex, genelist[i], genelist[j], graph);
	//thrPool.~thread_pool();
	unordered_map<uint32_t, pair<double, double>> formatedGraph = buildGraph(chains, genelist, n);
	if (formatedGraph.size())
	{
		double minx = 100, miny = 100;
		for (auto& key : formatedGraph)
		{
			if (key.second.first < minx)
				minx = key.second.first;
			if (key.second.second < miny)
				miny = key.second.second;
		}
		for (auto& key : formatedGraph)
		{
			key.second.first = key.second.first - minx;
			key.second.second = key.second.second - miny;
		}
		*vcount = formatedGraph.size();
		*ecount = chains.size();
		*vertex = new uint32_t[*vcount];
		*coordinates = new double[(*vcount) * 2];
		uint32_t i = 0, j = 0;
		for (auto& key : formatedGraph)
		{
			(*vertex)[i++] = key.first;
			(*coordinates)[j++] = key.second.first;
			(*coordinates)[j++] = key.second.second;
		}
		*edge = new uint32_t[(*ecount) * 2];
		i = 0;
		for (auto& item : chains)
		{
			(*edge)[i++] = item.first;
			(*edge)[i++] = item.second;
		}
	}
}

//clean up arrays allocated using C++ "new"
extern "C" __declspec(dllexport) void freeArray(void* arr)
{
	delete[] arr;
}




