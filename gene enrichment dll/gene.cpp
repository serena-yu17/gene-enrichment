// gene enrichment dll.cpp : Defines the exported functions for the DLL application.
//

#include "main.h"

using namespace std;

unique_ptr<GraphData> currentGraph;

bool getRunningStatus(uint32_t tid)
{
	std::shared_lock<std::shared_mutex> lock(runningMutex);
	auto found = running.find(tid);
	if (found == running.end())
		return false;
	return found->second;
}

void dijkstraSearch(uint32_t tid, PairSet* chains, int32_t node1, int32_t const genelist[],
	int32_t nQuery, WGraph& graph,
	unordered_map<int, Vertex>& idvertex, unordered_map<Vertex, int>& vertexid)
{
	vector<size_t> predecessors(boost::num_vertices(graph));
	vector<float> distances(boost::num_vertices(graph));
	auto indexmap = boost::get(boost::vertex_index, graph);
	auto weightmap = boost::get(boost::edge_weight, graph);
	boost::dijkstra_shortest_paths(graph,
		idvertex[node1],
		&predecessors[0],
		&distances[0],
		weightmap,
		indexmap,
		less<float>(),
		boost::closed_plus<float>(),
		infinityFloat,
		0,
		boost::default_dijkstra_visitor()
	);
	if (!getRunningStatus(tid))
		return;
	for (int i = 0; i < nQuery; i++)
		if (genelist[i] != node1 && distances[idvertex[genelist[i]]] != infinityFloat)
			for (Vertex vtx = idvertex[genelist[i]]; vtx != idvertex[node1]; vtx = predecessors[vtx])
				chains->insert(make_pair(vertexid[predecessors[vtx]], vertexid[vtx]));
}

//Use Fruchterman Reingold force directed layout to assign coordinates
std::unordered_map<int32_t, std::pair<float, float>> buildGraph(uint32_t tid, PairSet& chains)
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
	using Topology = boost::square_topology<mt19937>;
	using Position = Topology::point_type;
	vector<Position> positions(boost::num_vertices(graph));
	Topology topology;
	vector<Topology::point_difference_type> displacements(boost::num_vertices(graph));
	boost::random_graph_layout(
		graph,
		boost::make_iterator_property_map(positions.begin(), boost::identity_property_map{}),
		topology
	);
	std::unordered_map<int32_t, pair<float, float>> result;
	if (!getRunningStatus(tid))
		return result;
	//try which force model works better
	boost::fruchterman_reingold_force_directed_layout(
		graph,
		boost::make_iterator_property_map(positions.begin(), boost::identity_property_map{}),
		topology,
		attractive,
		repulsive,
		boost::all_force_pairs(),
		boost::linear_cooling<double>(100),
		make_iterator_property_map(displacements.begin(),
			boost::get(boost::vertex_index, graph),
			Topology::point_difference_type()
		)
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

DLLEXP void buildGraph(int32_t const data[], int32_t nData)
{
	int32_t maxWeight = 0;
	unordered_map<int32_t, Vertex> idvertex;
	unordered_map<Vertex, int32_t> vertexid;
	for (int32_t i = 0; i < nData; i++)
		if (data[i * 3 + 2] > maxWeight)
			maxWeight = data[i * 3 + 2];
	float fMaxWt = (float)log(maxWeight);
	WGraph graph;
	for (int32_t i = 0; i < nData; i++)
	{
		int32_t gene1 = data[i * 3];
		int32_t gene2 = data[i * 3 + 1];
		Vertex v1, v2;
		if (idvertex.find(gene1) != idvertex.end())
			v1 = idvertex[gene1];
		else
		{
			v1 = boost::add_vertex(graph);
			idvertex[gene1] = v1;
			vertexid[v1] = gene1;
		}
		if (idvertex.find(gene2) != idvertex.end())
			v2 = idvertex[gene2];
		else
		{
			v2 = boost::add_vertex(graph);
			idvertex[gene2] = v2;
			vertexid[v2] = gene2;
		}
		boost::add_edge(v1, v2, fMaxWt - log(data[i * 3 + 2]) + 1, graph);
	}
	currentGraph.reset(new GraphData(std::move(graph), std::move(idvertex), std::move(vertexid)));
}

//Build a graph based on coordinates. vertex, coordinates, vcount, edge, ecount are output values passed by reference. 
DLLEXP uint32_t enrichGenes(int32_t** vertex, float** coordinates, int32_t* vcount,
	int32_t** edge, int32_t* ecount, uint32_t gid, int32_t const genelist[], int32_t nQuery)
{
	if (!currentGraph)
		return 0;

	uint32_t tid = ++currentTid;
	{
		std::unique_lock<std::shared_mutex> runningLock(runningMutex);
		if (tid == 0)
			tid = ++currentTid;
		running[tid] = true;
	}

	auto& graph = currentGraph->graph;
	auto& idvertex = currentGraph->idvertex;
	auto& vertexid = currentGraph->vertexid;

	PairSet chains;

	for (int32_t i = 0; i < nQuery; i++)
	{
		if (!getRunningStatus(tid))
			return 0;
		dijkstraSearch(tid, &chains, genelist[i], genelist, nQuery, graph, idvertex, vertexid);
	}
	unordered_map<int32_t, pair<float, float>> formatedGraph = buildGraph(tid, chains);
	if (!getRunningStatus(tid))
		return 0;

	auto graphSize = formatedGraph.size();
	//Convert the output to C ABI
	if (graphSize)
	{
		//normalize distance
		float minx = 100, miny = 100;
		float avgDistance = 0, minDistance = (numeric_limits<float>::max)();
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
		vector<int32_t> vertexVec;
		vector<float> coordVec;
		unordered_set<int32_t> usedVertex;
		for (auto& key : formatedGraph)
		{
			vertexVec.push_back(key.first);
			usedVertex.insert(key.first);
			coordVec.push_back((key.second.first - minx) / (avgDistance + minDistance) * sqrt(graphSize) + 1);
			coordVec.push_back((key.second.second - miny) / (avgDistance + minDistance) * sqrt(graphSize) + 1);
		}

		float maxYShifted = 0, maxXShifted = 0;
		for (size_t i = 0; i < coordVec.size(); i += 2)
		{
			if (coordVec[i] > maxXShifted)
				maxXShifted = coordVec[i];
			if (coordVec[i + 1] > maxYShifted)
				maxYShifted = coordVec[i + 1];
		}

		//append unused verteces from the query, list at the bottom-left corner
		float appendY = maxYShifted + 1;
		float appendX = 0;
		for (int i = 0; i < nQuery; i++)
			if (usedVertex.find(genelist[i]) == usedVertex.end())
			{
				vertexVec.push_back(genelist[i]);
				coordVec.push_back(appendX);
				coordVec.push_back(appendY);
				appendX += 1;
				if (appendX > maxXShifted - 2)
				{
					appendX = 0;
					appendY += 1;
				}
			}

		*vcount = (int32_t)vertexVec.size();
		*vertex = new int32_t[*vcount];
		*coordinates = new float[coordVec.size()];

		memcpy(*vertex, &vertexVec[0], *vcount);
		for (size_t i = 0; i < vertexVec.size(); i++)
			(*vertex)[i] = vertexVec[i];
		for (size_t i = 0; i < coordVec.size(); i++)
			(*coordinates)[i] = coordVec[i];
		*ecount = (int32_t)chains.size();
		*edge = new int32_t[(*ecount) * 2];
		int i = 0;
		for (auto& item : chains)
		{
			(*edge)[i++] = item.first;
			(*edge)[i++] = item.second;
		}
	}
	return tid;
}

DLLEXP void freeIntArr(int** arr) {
	delete[](*arr);
}

DLLEXP void freeFloatArr(float** arr) {
	delete[](*arr);
}

DLLEXP void terminateProc(uint32_t tid)
{
	std::unique_lock<std::shared_mutex> lock(runningMutex);
	running[tid] = false;
}