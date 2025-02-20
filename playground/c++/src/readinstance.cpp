#include "readinstance.hpp"

pair<AdjacencyList, EdgeSet> read_instance()
{
    unsigned int n, m;    
    string s;
    string line;
    getline(cin, line);
    istringstream iss(line);
    iss>>s;
    iss>>s;
    iss >> n >> m;
    AdjacencyList adj(n);
    EdgeSet edges;

    while (m) {
        getline(cin, line);
        if (line.find("c ") != string::npos || line.empty())
            continue;
        Vertex u, v;
        istringstream iss(line);
        iss >> u >> v;
        adj[u-1].insert(v-1);
        adj[v-1].insert(u-1);
        edges.insert({min(u-1, v-1), max(u-1, v-1)});
        m--;
    }
    return {adj, edges};
}