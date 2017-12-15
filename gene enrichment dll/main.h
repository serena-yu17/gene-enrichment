#pragma once

#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <limits>
#include <random>

#include <boost/system/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/functional/hash.hpp>

typedef boost::property<boost::edge_weight_t, float> WeightProperty;
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property, WeightProperty, boost::listS> WGraph;
typedef boost::graph_traits<WGraph>::vertex_descriptor Vertex;
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property, boost::listS> UGraph;
typedef std::unordered_set<std::pair<int32_t, int32_t>, boost::hash<std::pair<int32_t, int32_t>>> PairSet;

constexpr float inf = (std::numeric_limits<float>::max)();
bool running[10064]{ 0 };

namespace std
{
	template<>
	struct hash<pair<int32_t, int32_t>>
	{
		size_t operator()(pair<int32_t, int32_t> const& pair) const
		{
			return (pair.first << 16) | pair.second;				   //
		}
	};
}

inline float distance(float x1, float x2, float y1, float y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

/////
/////
////