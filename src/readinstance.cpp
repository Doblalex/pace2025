#include "readinstance.hpp"

Instance* read_instance(globalprops* props) {
	unsigned int n, m;
	string s;
	string line;
	getline(cin, line);
	istringstream iss(line);
	iss >> s;
	iss >> s;
	iss >> n >> m;
	props->n = n;
	props->m = m;
	props->id_to_vertex.resize(n + 1);
	Instance* instance = new Instance(props);
	instance->n = n;
	instance->m = m;
	instance->G = new Graph();
	for (VertexCount i = 0; i < n; i++) {
		props->id_to_vertex[i + 1] = boost::add_vertex(*instance->G);
		(*instance->G)[props->id_to_vertex[i + 1]].id = i + 1;
	}

	while (m) {
		getline(cin, line);
		if (line.find("c ") != string::npos || line.empty()) {
			continue;
		}
		Vertex u, v;
		istringstream iss(line);
		iss >> u >> v;
		boost::add_edge(props->id_to_vertex[u], props->id_to_vertex[v], (*instance->G));
		(*instance->G)[props->id_to_vertex[u]].cnt_undominated_neighbors++;
		(*instance->G)[props->id_to_vertex[v]].cnt_undominated_neighbors++;
		m--;
	}
	return instance;
}