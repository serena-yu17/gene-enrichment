//clang gene.cpp -fPIC -shared -std=c++17 -O3 -o gene.dll
#pragma once

#if defined(_MSC_VER)
//  MSVC
#define DLLEXP extern "C" __declspec(dllexport)
#elif defined(__GNUC__)
//  GCC
#define DLLEXP extern "C"  __attribute__((visibility("default")))
#elif defined(__clang__)
//clang
#define DLLEXP extern "C"  __attribute__((visibility("default")))
#else
#define DLLEXP extern "C"
#pragma warning Unknown dynamic link import/export semantics
#endif

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>
#include <random>
#include <shared_mutex>
#include <atomic>
#include <cmath>
#include <memory>

#include <boost/system/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/functional/hash.hpp>

using WeightProperty = boost::property<boost::edge_weight_t, float>;
using WGraph = boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property, WeightProperty, boost::listS>;
using Vertex = boost::graph_traits<WGraph>::vertex_descriptor;
using UGraph = boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property, boost::listS>;
using PairSet = std::unordered_set<std::pair<int32_t, int32_t>, boost::hash<std::pair<int32_t, int32_t>>>;

constexpr float infinityFloat = (std::numeric_limits<float>::max)();

std::unordered_map<int32_t, bool> running;
std::shared_mutex runningMutex;

class GraphData {
public:
	WGraph graph;
	std::unordered_map<int32_t, Vertex> idvertex;
	std::unordered_map<Vertex, int32_t> vertexid;

	GraphData(WGraph&& graph,
		std::unordered_map<int32_t, Vertex>&& idvertex,
		std::unordered_map<Vertex, int32_t>&& vertexid)
	{
		this->graph = std::move(graph);
		this->idvertex = std::move(idvertex);
		this->vertexid = std::move(vertexid);
	}
};

inline float distance(float x1, float x2, float y1, float y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

inline float repulsive(UGraph::vertex_descriptor& descr1, UGraph::vertex_descriptor& descr2, float k, float d, UGraph const& graph)
{
	if (d < (k / 3.5))
		return infinityFloat;
	else
		return k * k / d;
}

inline float attractive(UGraph::edge_descriptor& descr, float k, float d, UGraph const& graph)
{
	if (d < (k / 3.5))
		return 0;
	else
		return d * d / k;
}