import networkx as nx
import sys
sys.setrecursionlimit(10**6)

_, _, n, m = input().split()
n, m = int(n), int(m)
G = nx.Graph()

for _ in range(m):
    u, v = map(int, input().split())
    G.add_edge(u, v)

tw, _ = nx.approximation.treewidth.treewidth_min_degree(G)
print("Treewidth approximation:", tw)
print("Is planar?", nx.is_planar(G))
print("Vertex cover approximation", len(nx.approximation.vertex_cover.min_weighted_vertex_cover(G)))
print("node connectivity:", nx.node_connectivity(G))
print("degeneracy", max(nx.algorithms.core.core_number(G).values()))
print("dominating set approximation", len(nx.approximation.dominating_set.min_weighted_dominating_set(G)))
print("independent set approximation", len(nx.approximation.maximum_independent_set(G)))
print("chromatic number", nx.algorithms.coloring.greedy_color(G))
print("clique number", nx.approximation.max_clique(G))
edges = nx.minimum_spanning_edges(G)
print("Feedback edges set size", len(G.edges) - len(list(edges)))
