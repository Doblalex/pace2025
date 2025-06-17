#include "ogdf_instance.hpp"

void Instance::read(std::istream& is, std::vector<ogdf::node>& ID2node) {
	unsigned int n, m;
	std::string s;
	std::string line;
	while (line.empty() || line[0] == 'c') {
		if (!is.good()) {
			std::cerr << "Error reading input" << std::endl;
			std::exit(1);
		}
		getline(is, line);
	}
	std::istringstream iss(line);
	iss >> s;
	if (s != "p") {
		std::cerr
				<< "Bad header line not matching `p (ds|hs) [0-9]+ [0-9]+`: " << line.substr(0, 100)
				<< std::endl;
		std::exit(1);
	}
	iss >> type;
	iss >> n >> m;
	if (type == "ds") {
		read_DS(is, ID2node, n, m);
	} else if (type == "hs") {
		read_HS(is, ID2node, n, m);
	} else {
		std::cerr << "Unknown input type " << type << std::endl;
		std::exit(1);
	}
}

void Instance::read_DS(std::istream& is, std::vector<ogdf::node>& ID2node, unsigned int n,
		unsigned int m) {
	std::string line;
	clear();
	ID2node.clear();
	ID2node.reserve(n + 1);
	ID2node.push_back(nullptr);
	for (int i = 1; i <= n; i++) {
		auto n = G.newNode(i);
		ID2node.push_back(n);
		node2ID[n] = i;
	}
	maxid = n;
	for (int i = 0; i < m; i++) {
		line.clear();
		while (line.empty() || line[0] == 'c') {
			if (!is.good()) {
				std::cerr << "Error reading input" << std::endl;
				std::exit(1);
			}
			getline(is, line);
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

void Instance::read_HS(std::istream& is, std::vector<ogdf::node>& ID2node, unsigned int n,
		unsigned int m) {
	std::string line;
	clear();
	ID2node.clear();
	ID2node.reserve(n + m + 1);
	ID2node.push_back(nullptr);
	for (int i = 1; i <= n + m; i++) {
		auto node = G.newNode(i);
		ID2node.push_back(node);
		node2ID[node] = i;
		if (i <= n) {
			is_dominated[node] = true;
		} else {
			is_subsumed[node] = true;
		}
	}
	for (int i = 0; i < m; i++) {
		line.clear();
		while (line.empty() || line[0] == 'c') {
			if (!is.good()) {
				std::cerr << "Error reading input" << std::endl;
				std::exit(1);
			}
			getline(is, line);
		}
		std::istringstream iss(line);
		int u;
		while (iss >> u) {
			auto e = G.newEdge(ID2node[u], ogdf::Direction::before, ID2node[i + n + 1],
					ogdf::Direction::after);
			reverse_edge[e] = nullptr;
		}
	}
}
