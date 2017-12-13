// gene enrichment dll.cpp : Defines the exported functions for the DLL application.
//
#include "main.h"

using namespace std;

void dijkstraSearch(int32_t tid, unordered_set<pair<int32_t, int32_t>>* chains, int32_t node1, int32_t const genelist[], int32_t nQuery, WGraph& graph)
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
		inf,
		0,
		boost::default_dijkstra_visitor()
	);	
	if (!running[tid])
		return;
	for (int i = 0; i < nQuery; i++)
		if (genelist[i] != node1 && distances[genelist[i]] != inf)
			for (size_t vtx = genelist[i]; vtx != node1 && running[tid]; vtx = predecessors[vtx])
				chains->insert(make_pair(predecessors[vtx], vtx));
}

//Use Fruchterman Reingold force directed layout to assign coordinates
unordered_map<int32_t, pair<float, float>> buildGraph(int32_t tid, unordered_set<pair<int32_t, int32_t>> chains)
{
	UGraph graph;
	unordered_map<int32_t, size_t> idVertex;
	unordered_map<size_t, int32_t> vertexID;
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
	using Topology = boost::square_topology<boost::mt19937>;
	using Position = Topology::point_type;
	vector<Position> positions(boost::num_vertices(graph));
	Topology topology;
	boost::random_graph_layout(
		graph,
		boost::make_iterator_property_map(positions.begin(), boost::identity_property_map{}),
		topology
	);
	std::unordered_map<int32_t, pair<float, float>> result;
	if (!running[tid])
		return result;
	boost::fruchterman_reingold_force_directed_layout(
		graph,
		boost::make_iterator_property_map(positions.begin(), boost::identity_property_map{}),
		topology
		,boost::attractive_force([](UGraph::edge_descriptor, float k, float d, UGraph const&) { return log10(d) / k; })
	);
	auto index_map = boost::get(boost::vertex_index, graph);
	auto begin = positions.begin();
	boost::iterator_property_map<vector<Position>::iterator, decltype(index_map)> position_map(positions.begin(), index_map);
	BGL_FORALL_VERTICES(v, graph, UGraph)
	{
		Position pos = position_map[v];
		result[vertexID[v]] = std::make_pair(pos[0], pos[1]);
	}
	return result;
}

//Build a graph based on coordinates. vertex, coordinates, vcount, edge, ecount are output values passed by reference. 
extern "C" __declspec(dllexport) void enrichGenes(int32_t tid, int32_t** vertex, float** coordinates, int32_t* vcount, int32_t** edge, int32_t* ecount, int32_t const data[], int32_t nData, int32_t const genelist[], int32_t nQuery)
{
	running[tid] = 1;
	int32_t maxGene = 0;
	int32_t maxWeight = 0;
	for (int32_t i = 0; i < nData; i++)
	{
		int32_t j = i * 3;
		if (data[j] > maxGene)
			maxGene = data[j];
		if (data[j + 1] > maxGene)
			maxGene = data[j + 1];
		if (data[j + 2] > maxWeight)
			maxWeight = data[j + 2];
	}
	WGraph graph;
	for (int32_t i = 0; i < maxGene + 1; i++)
		boost::add_vertex(graph);
	float fMaxWt = log2(maxWeight);
	for (int32_t i = 0; i < nData; i++)
	{
		int32_t j = i * 3;
		boost::add_edge((Vertex)data[j], (Vertex)data[j + 1], fMaxWt - log2(data[j + 2]) + 1, graph);
	}
	unordered_set<pair<int32_t, int32_t>> chains;
	for (int32_t i = 0; i < nQuery; i++)
	{
		if (!running[tid])
			return;
		dijkstraSearch(tid, &chains, genelist[i], genelist, nQuery, graph);
	}
	unordered_map<int32_t, pair<float, float>> formatedGraph = buildGraph(tid, chains);
	if (!running[tid])
		return;

	//Transform to C ABI
	if (formatedGraph.size())
	{
		float minx = 100, miny = 100;
		float avgDistance = 0, minDistance = numeric_limits<int32_t>::max();
		for (auto& key : chains)
		{
			float dist = distance(formatedGraph[key.first].first, formatedGraph[key.first].second, formatedGraph[key.second].first, formatedGraph[key.second].second);
			avgDistance += dist;
			if (dist < minDistance)
				minDistance = dist;
		}
		avgDistance /= chains.size();
		for (auto& key : formatedGraph)
		{
			if (key.second.first < minx)
				minx = key.second.first;
			if (key.second.second < miny)
				miny = key.second.second;
		}
		*vcount = formatedGraph.size();
		*ecount = chains.size();
		*vertex = new int32_t[*vcount];
		*coordinates = new float[(*vcount) * 2];
		int32_t i = 0, j = 0;
		for (auto& key : formatedGraph)
		{
			(*vertex)[i++] = key.first;
			(*coordinates)[j++] = (key.second.first - minx) / (avgDistance + minDistance) * sqrt(*vcount) + 1;
			(*coordinates)[j++] = (key.second.second - miny) / (avgDistance + minDistance) * sqrt(*vcount) + 1;
		}
		*edge = new int32_t[(*ecount) * 2];
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

extern "C" __declspec(dllexport) void terminnateProc(int32_t tid)
{
	running[tid] = 0;
}


