import networkx as nx
import sys
sys.setrecursionlimit(10**6)
import os
import csv
import sys
import time
import signal
from timeout_decorator import timeout

@timeout(20)
def treewidth(G):
    tw, _ =  nx.approximation.treewidth.treewidth_min_degree(G)
    return tw

@timeout(20)
def is_planar(G):
    return nx.is_planar(G)

@timeout(20)
def vertex_cover(G):
    return len(nx.approximation.vertex_cover.min_weighted_vertex_cover(G))

@timeout(20)
def degeneracy(G):
    return max(nx.algorithms.core.core_number(G).values())

@timeout(20)
def chromatic_number(G):
    return max(nx.algorithms.coloring.greedy_color(G).values())

@timeout(20)
def feedback_edges(G):
    edges = nx.minimum_spanning_edges(G)
    return len(G.edges) - len(list(edges))

DIR = "../../instances/official/ds/exact"
csvwriter = csv.writer(open("graph_properties.csv", "w"))
csvwriter.writerow(["instance", "n", "m", "TW UB", "Planar", "VC UB", "Degeneracy", "Chromatic UB", "FES UB"])
for file in sorted(os.listdir(DIR)):
    if file.endswith(".gr"):
        print(file)
        with open(os.path.join(DIR, file)) as f:
            lines = list(f.readlines())
            _, _, n, m = lines[0].split()
            n, m = int(n), int(m)
            G = nx.Graph()

            for i in range(m):
                u, v = map(int, lines[i+1].split())
                G.add_edge(u, v)
            try: 
                tw = treewidth(G)
            except:
                tw = ""
            try: 
                ip = is_planar(G)
            except:
                ip = ""
            try:
                vc = vertex_cover(G)
            except:
                vc = ""
            try:
                deg = degeneracy(G)
            except:
                deg = ""
            try:
                cn = chromatic_number(G)
            except: 
                cn = ""
            try:
                fes = feedback_edges(G)
            except:
                fes = ""
            csvwriter.writerow([file, len(G.nodes), len(G.edges), tw, ip, vc, deg, cn, fes])
