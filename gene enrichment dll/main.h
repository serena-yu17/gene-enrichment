#pragma once

#include <windows.h>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <atomic>
#include <mutex>
#include <cmath>
#include <stdint.h>

#include <memory>
#include <boost/system/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/fruchterman_reingold.hpp>
#include <boost/graph/random_layout.hpp>
#include <libs/graph/src/read_graphviz_new.cpp>
#include <boost/graph/topology.hpp>
#include <boost/random.hpp>

#include "ctpl_stl.h"

std::unordered_map<std::string, int> geneid;
std::vector<std::string> idgene;


struct Node
{
	int val;
	Node* next;
	Node() {}
	Node(int val, Node* next)
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
			return hash<int>{}(node->val);
		}
	};

	template<>
	struct hash<pair<int, int>>
	{
		size_t operator()(pair<int, int> const& pair)	const
		{
			return (pair.first << 16) + pair.second;
		}
	};
}

struct ComparepNode
{
	bool operator()(Node const* node1, Node const* node2)  const
	{
		return node1->val == node2->val;
	}
};

using Graph = boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS, boost::vecS>;

