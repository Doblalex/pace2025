import networkx as nx
import matplotlib.pyplot as plt

_, _, n, m = input().split()
n, m = int(n), int(m)
G = nx.Graph()

for _ in range(m):
    u, v = map(int, input().split())
    G.add_edge(u, v)

deg2neighs = set()
rest = set()

for n in G.nodes:
    if G.degree(n) > 800:
        for neigh in G.neighbors(n):
            if G.degree(neigh) != 2 and G.degree(neigh) < 800:
                rest.add(n)
                break
        else:
            deg2neighs.add(n)
print(len(deg2neighs), len(rest))
botha = 0
bothb = 0
diff = 0
Gvc = nx.Graph()
for n in G.nodes:
    if G.degree(n) == 2:
        neighs = list(G.neighbors(n))
        a,b = neighs[0] in deg2neighs, neighs[1] in deg2neighs

        if a and b: 
            botha += 1
            Gvc.add_edge(neighs[0], neighs[1])
        elif a != b:
            diff += 1
            if abs(neighs[0]-neighs[1]) != 1:
                print(neighs)
        else:
            bothb += 1
nx.write_adjlist(Gvc, path="vcgraph.gr")

            




# G = nx.ego_graph(G, 1, radius=6)

# nx.draw(G, node_size=5)
# plt.show()