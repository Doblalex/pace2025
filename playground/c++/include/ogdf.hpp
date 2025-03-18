#pragma once

#include "ogdf/basic/Graph.h"

extern ogdf::Logger logger;

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

    template<typename NL, typename EL>
    void initFrom(const Instance &other,
                  const NL &nodes,
                  const EL &edges,
                  const ogdf::NodeArray<ogdf::node> &nMap,
                  const ogdf::EdgeArray<ogdf::edge> &eMap) {
        for (auto n : nodes) {
            node2ID[nMap[n]] = other.node2ID[n];
            is_dominated[nMap[n]] = other.is_dominated[n];
            is_subsumed[nMap[n]] = other.is_subsumed[n];
        }
        for (auto e : edges) {
            reverse_edge[eMap[e]] = other.reverse_edge[e] != nullptr ? eMap[other.reverse_edge[e]] : nullptr;
        }
    }

    void read(std::istream &is);

    void dumpBCTree();

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

    void addToDominatingSet(ogdf::node v) {
        ogdf::Graph::node_iterator it;
        addToDominatingSet(v, it);
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

    bool reductionExtremeDegrees();

    std::list<Instance> decomposeConnectedComponents();

    bool reductionBCTree();

    std::pair<size_t, size_t> dominationStats() {
        size_t can = 0, needs = 0;
        for (ogdf::node n : G.nodes) {
            if (!is_dominated[n]) {
                ++needs;
            }
            if (!is_subsumed[n]) {
                ++can;
            }
        }
        return std::make_pair(can, needs);
    }
};

void reduceAndSolve(Instance &I, int d = 0);

void solveEvalMaxSat(Instance &I);
