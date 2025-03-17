import re
from pathlib import Path

from ogdf_python import * # pip install ogdf-python[quickstart]

cppinclude("ogdf/decomposition/BCTree.h")
cppinclude("ogdf/energybased/FMMMLayout.h")
cppinclude("ogdf/layered/SugiyamaLayout.h")

DIR = Path(__file__).parent / "../../instances/official/ds/exact"
DIR = DIR.absolute()

def readGR(G, p):
    N = None
    m = 0
    prob = None
    with open(p) as f:
        for l in f:
            l = l.strip()
            if not l or l[0] == "c":
                continue
            elif l[0] == "p":
                r = re.match("p (.*) ([0-9]+) ([0-9]+)", l)
                prob = r.group(1)
                N = [None] + [G.newNode() for _ in range(int(r.group(2)))]
                m = int(r.group(3))
            else:
                u, v = map(int, l.split())
                G.newEdge(N[u], N[v])
    assert G.numberOfNodes() == len(N) - 1
    assert G.numberOfEdges() == m
    return prob


def drawBCTree(BC):
    BCA = ogdf.GraphAttributes(BC.bcTree(), ogdf.GraphAttributes.all)
    for node in BC.bcTree().nodes:
        if BC.typeOfBNode(node) == ogdf.BCTree.BNodeType.BComp:
            BCA.shape[node] = ogdf.Shape.Ellipse
            # BCA.label[node] = f"{node.index()}: {BC.numberOfEdges(node)}"
            BCA.label[node] = str(BC.numberOfNodes(node))
        else:
            BCA.shape[node] = ogdf.Shape.Rhomb
            BCA.label[node] = ""  # f"{node.index()}"
    return BCA


G = ogdf.Graph()
for p in DIR.glob("*.gr"):
    G.clear()
    svg = DIR / (p.stem + ".fb.svg")
    if svg.exists(): continue
    readGR(G, p)
    print(p, "G", G.numberOfNodes(), G.numberOfEdges())
    if G.numberOfEdges() > 30_000: continue 
    BC = ogdf.BCTree(G)
    print(p, "BC", BC.numberOfBComps(), BC.numberOfCComps())
    BCA = drawBCTree(BC)
    ogdf.FMMMLayout().call(BCA)
    ogdf.GraphIO.write(BCA, str(svg))
    ogdf.GraphIO.write(BCA, str(DIR / (p.stem + ".bc.gml")))
    ogdf.SugiyamaLayout().call(BCA)
    ogdf.GraphIO.write(BCA, str(DIR / (p.stem + ".sl.svg")))
    del BCA
    del BC
