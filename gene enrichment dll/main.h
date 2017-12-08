#pragma once

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <atomic>
#include <atomic>
#include <cmath>
#include <stdint.h>
#include <limits>
#include <memory>

#include <boost/system/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <boost/graph/topology.hpp>
#include <boost/random.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "ctpl_stl.h"

typedef boost::property<boost::edge_weight_t, float> WeightProperty;
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property, WeightProperty, boost::listS> WGraph;
typedef boost::graph_traits<WGraph>::vertex_descriptor Vertex;
typedef boost::adjacency_list <boost::vecS, boost::vecS, boost::directedS, boost::no_property, boost::no_property, boost::listS> UGraph;

std::unordered_map<std::string, int32_t> geneid;
std::vector<std::string> idgene;
const float inf = std::numeric_limits<float>::max();
std::atomic_flag spinLockChains = ATOMIC_FLAG_INIT;

struct Node
{
	int32_t val;
	Node* next;
	Node(int32_t val, Node* next)
	{
		this->val = val;
		this->next = next;
	}
	bool operator==(Node const& other) const
	{
		return this->val == other.val;
	}
	bool operator!=(Node const& other) const
	{
		return this->val != other.val;
	}
};

namespace std
{
	template<>
	struct hash<Node*>
	{
		size_t operator()(Node const* node) const
		{
			return hash<int32_t>{}(node->val);
		}
	};

	template<>
	struct hash<pair<int32_t, int32_t>>
	{
		size_t operator()(pair<int32_t, int32_t> const& pair) const
		{
			return (pair.first << 16) | pair.second;
		}
	};
}

struct ComparepNode
{
	bool operator()(Node const* node1, Node const* node2) const
	{
		return node1->val == node2->val;
	}
};

inline float distance(float x1, float x2, float y1, float y2)
{
	return sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1));
}

/////
/////
////