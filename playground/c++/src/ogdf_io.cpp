#include "ogdf_instance.hpp"

void Instance::read(std::istream& is, std::vector<ogdf::node>& ID2node) {
	unsigned int n, m;
	std::string s;
	std::string line;
	getline(is, line);
	std::istringstream iss(line);
	iss >> s;
	OGDF_ASSERT(s == "p");
	iss >> s;
	OGDF_ASSERT(s == "ds");
	iss >> n >> m;

	clear();
	ID2node.clear();
	ID2node.reserve(n + 1);
	ID2node.push_back(nullptr);
	for (int i = 1; i <= n; i++) {
		auto n = G.newNode(i);
		ID2node.push_back(n);
		node2ID[n] = i;
	}
	for (int i = 0; i < m; i++) {
		getline(is, line);
		if (line.empty() || line[0] == 'c') {
			continue;
		}
		int u, v;
		std::istringstream iss(line);
		iss >> u >> v;
		// sources front, targets tail
		auto e = G.newEdge(ID2node[u], ogdf::Direction::before, ID2node[v], ogdf::Direction::after);
		auto f = G.newEdge(ID2node[v], ogdf::Direction::before, ID2node[u], ogdf::Direction::after);
		reverse_edge[e] = f;
		reverse_edge[f] = e;
	}
	OGDF_ASSERT(G.numberOfNodes() == n);
	OGDF_ASSERT(G.numberOfEdges() == m * 2);
}
