import networkx as nx
import matplotlib.pyplot as plt

_, _, n, m = input().split()
n, m = int(n), int(m)
G = nx.Graph()

for _ in range(m):
    u, v = map(int, input().split())
    G.add_edge(u, v)

G = nx.ego_graph(G, 1, radius=6)

nx.draw(G, node_size=5)
plt.show()