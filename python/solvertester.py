import networkx as nx
import matplotlib.pyplot as plt
from ortools.linear_solver import pywraplp

_, _, n, m = input().split()
n, m = int(n), int(m)
G = nx.Graph()

for _ in range(m):
    u, v = map(int, input().split())
    G.add_edge(u, v)

solver = pywraplp.Solver.CreateSolver("CP-SAT")
if not solver:
    print("Could not create solver GLOP")
    exit(1)
solver.EnableOutput()

vars = {}
for u in G.nodes:
    vars[u] = solver.BoolVar("")

for u in G.nodes:
    s = [vars[u]]
    for v in G.neighbors(u):
        s.append(vars[v])
    solver.Add(solver.Sum(s) >= 1)

solver.Minimize(solver.Sum(vars.values()))
solver.Solve()

