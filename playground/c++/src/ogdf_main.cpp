#include <chrono>
#include <unordered_set>

#include "ogdf/basic/Graph_d.h"
#include "ogdf/basic/GraphSets.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/fileformats/GraphIO.h"
#include "ogdf/layered/SugiyamaLayout.h"

ogdf::Logger logger;

#ifdef OGDF_DEBUG
#define log if (true) logger.lout()
#else
#define log if (false) logger.lout()
#endif

#define forAllOutAdj(v, adj) \
    OGDF_ASSERT((v)->outdeg() == 0 || (v)->adjEntries.head()->isSource()); \
    for (auto adj_it = (v)->adjEntries.begin(); adj_it != (v)->adjEntries.end() && (*adj_it)->isSource();) { \
        auto adj = *adj_it; \
        ++adj_it;

#define forAllInAdj(v, adj) \
    OGDF_ASSERT((v)->indeg() == 0 || !(v)->adjEntries.tail()->isSource()); \
    for (auto adj_it = (v)->adjEntries.rbegin(); adj_it != (v)->adjEntries.rend() && !(*adj_it)->isSource();) { \
        auto adj = *adj_it; \
        ++adj_it;

struct Instance {
    ogdf::Graph G;
    // std::vector<ogdf::node> ID2node;
    ogdf::NodeArray<int> node2ID;

    std::vector<int> DS;
    ogdf::NodeArray<bool> is_dominated;
    ogdf::NodeArray<bool> is_subsumed;
    ogdf::EdgeArray<ogdf::edge> reverse_edge;

    Instance() : node2ID(G, -1),
                 is_dominated(G, false),
                 is_subsumed(G, false),
                 reverse_edge(G, nullptr) {
    }

    OGDF_NO_COPY(Instance)
    OGDF_NO_MOVE(Instance)

    void clear() {
        G.clear();
        DS.clear();
        node2ID.init(G, -1);
        is_dominated.init(G, false);
        is_subsumed.init(G, false);
        reverse_edge.init(G, nullptr);
    }

    void read(std::istream &is) {
        unsigned int n, m;
        std::string s;
        std::string line;
        getline(is, line);
        std::istringstream iss(line);
        iss >> s;
        OGDF_ASSERT(s == "p");
        iss >> s;
        OGDF_ASSERT(s == "ds");
        iss >> n >> m;

        clear();
        std::vector<ogdf::node> ID2node; // ID2node.clear();
        ID2node.reserve(n + 1);
        ID2node.push_back(nullptr);
        for (int i = 1; i <= n; i++) {
            auto n = G.newNode(i);
            ID2node.push_back(n);
            node2ID[n] = i;
        }
        for (int i = 0; i < m; i++) {
            getline(is, line);
            if (line.empty() || line[0] == 'c') continue;
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

    void dumpBCTree() {
        ogdf::BCTree BC(G);
        ogdf::GraphAttributes BCA(BC.bcTree(), ogdf::GraphAttributes::all);
        for (auto node : BC.bcTree().nodes) {
            if (BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp) {
                BCA.shape(node) = ogdf::Shape::Ellipse;
                BCA.label(node) = std::to_string(BC.numberOfNodes(node));
                BCA.width(node) = 15 + 5 * BCA.label(node).size();
            } else {
                BCA.shape(node) = ogdf::Shape::Rhomb;
                BCA.label(node) = "";
            }
        }
        std::string stamp = std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + "-" +
            std::to_string(G.numberOfNodes());
        ogdf::FMMMLayout().call(BCA);
        ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".svg");
        ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".gml");
        ogdf::SugiyamaLayout().call(BCA);
        ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".svg");
        ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".gml");
    }

    void safeDelete(ogdf::node n, ogdf::Graph::node_iterator &it) {
        if (*it == n) ++it; // no need to check for end, as *end == nullptr
        safeDelete(n);
    }

    void safeDelete(ogdf::node n) {
        // ID2node[node2ID[n]] = nullptr;
        G.delNode(n);
        G.consistencyCheck();
    }

    void safeDelete(ogdf::edge e) {
        ogdf::edge r = reverse_edge[e];
        if (r != nullptr) {
            OGDF_ASSERT(reverse_edge[r] == e);
            reverse_edge[r] = nullptr;
        }
        G.delEdge(e);
    }

    void addToDominatingSet(ogdf::node v, ogdf::Graph::node_iterator &it) {
        OGDF_ASSERT(!is_subsumed[v]);
        DS.push_back(node2ID[v]);
        forAllOutAdj(v, adj)
            auto u = adj->twinNode();
            markDominated(u);
            if (u->outdeg() == 0) {
                safeDelete(u, it);
            }
        }
        safeDelete(v, it);
    }

    void markDominated(ogdf::node v) {
        is_dominated[v] = true;
        forAllInAdj(v, adj)
            safeDelete(adj->theEdge());
            // TODO shortcut some rules?
        }
        // G.consistencyCheck();
    }

    void markSubsumed(ogdf::node v) {
        is_subsumed[v] = true;
        forAllOutAdj(v, adj)
            safeDelete(adj->theEdge());
            // TODO shortcut some rules?
        }
        // G.consistencyCheck();
    }

    bool reductionExtremeDegrees() {
        bool reduced = false, has_undom = false;
        ogdf::node universal_in = nullptr;
        int i = G.numberOfNodes();
        for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
            OGDF_ASSERT(i >= 0); // protect against endless loops
            auto n = *it;
            ++it; // increment now so that we can safely delete n

            OGDF_ASSERT(!is_dominated[n] || n->indeg() == 0);
            OGDF_ASSERT(!is_subsumed[n] || n->outdeg() == 0);
            OGDF_ASSERT(!is_subsumed[n] || n->indeg() > 0 || is_dominated[n]);

            // isolated
            if (n->indeg() == 0 && !is_dominated[n]) {
                addToDominatingSet(n, it);
                reduced = true;
                continue;
            }
            if (n->outdeg() == 0 && is_dominated[n]) {
                safeDelete(n, it);
                reduced = true;
                continue;
            }

            // antenna
            if (n->indeg() == 1 && (n->outdeg() == 0 ||
                (n->outdeg() == 1 && n->adjEntries.head()->twinNode() == n->adjEntries.tail()->twinNode()))) {
                auto u = n->adjEntries.head()->twinNode();
                OGDF_ASSERT(u!=n);
                addToDominatingSet(u, it);
                reduced = true;
                continue;
            }
            if (n->outdeg() == 1 && is_dominated[n]) {
                OGDF_ASSERT(n->adjEntries.head()->isSource());
                auto u = n->adjEntries.head()->twinNode();
                if (!is_subsumed[u] || u->indeg() >= 2) {
                    safeDelete(n, it);
                } else {
                    addToDominatingSet(n, it);
                }
                reduced = true;
                continue;
            }

            // universal
            if (n->outdeg() == G.numberOfNodes() - 1) {
                addToDominatingSet(n, it);
                G.clear();
                return true;
            }
            if (n->indeg() == G.numberOfNodes() - 1) {
                universal_in = n;
            } else if (!has_undom && !is_dominated[n]) {
                has_undom = true;
            }
        }

        if (universal_in != nullptr) {
            auto it = G.nodes.end();
            if (has_undom) {
                safeDelete(universal_in, it);
            } else {
                addToDominatingSet(universal_in, it);
            }
            reduced = true;
        }

        return reduced;
    }

    std::list<Instance> decomposeConnectedComponents() {
        ogdf::Graph::CCsInfo CC(G);
        std::list<Instance> instances;
        if (CC.numberOfCCs() > 1) {
            ogdf::NodeArray<ogdf::node> nMap(G, nullptr); // can safely be reused for different CCs
            ogdf::EdgeArray<ogdf::edge> eMap(G, nullptr);
            for (int i = 0; i < CC.numberOfCCs(); ++i) {
                instances.emplace_back();
                Instance &I = instances.back();
                I.G.insert(CC, i, nMap, eMap);
                for (auto n : CC.nodes(i)) {
                    I.node2ID[nMap[n]] = node2ID[n];
                    I.is_dominated[nMap[n]] = is_dominated[n];
                    I.is_subsumed[nMap[n]] = is_subsumed[n];
                }
                for (auto e : CC.edges(i)) {
                    I.reverse_edge[eMap[e]] = reverse_edge[e] != nullptr ? eMap[reverse_edge[e]] : nullptr;
                }
            }
        }
        return instances;
    }
};

void reduceAndSolve(Instance &I, int d = 0) {
    bool changed = true;
    int m, n, i = 0;
    while (changed) {
        n = I.G.numberOfNodes();
        m = I.G.numberOfEdges();
        changed = false;
        log << "Reduce iteration " << i << " depth " << d << ": " << n << " nodes, " << m << " edges" << std::endl;

        // this reduction is so cheap, make sure we really have no isolated vertices before decomposing components
        while (I.reductionExtremeDegrees()) changed = true;
        // TODO reductionTwins(instance)

        std::list<Instance> comps = I.decomposeConnectedComponents();
        if (!comps.empty()) {
            log << comps.size() << " connected components" << std::endl;
            I.clear(); // save some memory

            int c = 0;
            for (auto &comp : comps) {
                log << "Connected component " << c << std::endl;
                ogdf::Logger::Indent _(logger);
                ++c;

                // run remaining reduction rules on all components
                // TODO reductionDomination(instance)
                // TODO reductionDominationPaper(instance, dominatingSet))

                // and now recurse
                reduceAndSolve(comp, d + 1);
                I.DS.insert(I.DS.end(), comp.DS.begin(), comp.DS.end());
            }
            return;
        }

        // TODO reductionDomination(instance)
        // TODO reductionDominationPaper(instance, dominatingSet))

        OGDF_ASSERT(!changed || I.G.numberOfNodes() < n|| I.G.numberOfEdges() < m);
        ++i;
    }

    if (I.G.numberOfNodes() < 1) {
        return;
    }

    // now to solving...
    I.dumpBCTree();
    // TODO:
    // debug(instance->CntCanBeDominatingSet(), instance->CntNeedsDomination(), instance->n);
    // #ifdef USE_ORTOOLS
    //     solveCPSat(instance, dominatingSet);
    // #else
    //     solveEvalMaxSat(instance, dominatingSet);
    // #endif
}

int main(int argc, char **argv) {
    auto start = std::chrono::system_clock::now();
    Instance I;
    I.read(std::cin);
    reduceAndSolve(I);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> runtime = end - start;
}
