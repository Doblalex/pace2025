import networkx as nx

G: nx.Graph = nx.read_adjlist("vcgraph.gr")

# tw, _ =  nx.approximation.treewidth.treewidth_min_degree(G)
# print(tw)
# print(nx.approximation.connectivity.node_connectivity(G))
# print(len(nx.approximation.vertex_cover.min_weighted_vertex_cover(G)))
# print(len(nx.approximation.dominating_set.min_weighted_dominating_set(G)))
# deg = max(G.degree(n) for n in G.nodes)
# print(deg)
# print(nx.core_number(G))
print(len(nx.approximation.max_clique(G)))
print(len(nx.approximation.maximum_independent_set(G)))
print(nx.approximation.diameter(G))
print(nx.diameter(G))