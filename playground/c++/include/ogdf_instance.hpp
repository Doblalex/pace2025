#pragma once

#include "ogdf_util.hpp"

struct Instance {
private:
    bool subsumptionCondition1(const ogdf::node& u, const ogdf::node& v, const ogdf::NodeArray<bool>& adju);
    bool subsumptionCondition2(const ogdf::node& u, const ogdf::node& v, const ogdf::NodeArray<bool>& adju);
    bool subsumptionCondition3(const ogdf::node& u, const ogdf::node& v, ogdf::NodeArray<bool>& inadjv);
public:
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
        // DS.clear(); Do not clear DS! This value is still used after computing connected components
        node2ID.init(G, -1);
        is_dominated.init(G, false);
        is_subsumed.init(G, false);
        reverse_edge.init(G, nullptr);
    }

    bool checkNode(ogdf::node n) {
        OGDF_ASSERT(!is_dominated[n] || n->indeg() == 0);
        OGDF_ASSERT(!is_subsumed[n] || n->outdeg() == 0);
        OGDF_ASSERT(!is_subsumed[n] || n->indeg() > 0 || is_dominated[n]);
        OGDF_ASSERT(n->outdeg() == 0 || n->adjEntries.head()->isSource());
        OGDF_ASSERT(n->indeg() == 0 || !n->adjEntries.tail()->isSource());
#ifdef OGDF_HEAVY_DEBUG
        size_t c = 0;
        for (auto adj_it = n->adjEntries.begin(); adj_it != n->adjEntries.end(); ++c) {
            auto adj = *adj_it;
            ++adj_it;
            OGDF_ASSERT(adj->isSource() == (c < n->outdeg()));
        }
#endif
        return true;
    }

    template<typename NL, typename EL>
    void initFrom(const Instance &other,
                  const NL &nodes,
                  const EL &edges,
                  const ogdf::NodeArray<ogdf::node> &nMap,
                  const ogdf::EdgeArray<ogdf::edge> &eMap,
                  std::function<ogdf::node(ogdf::node)> originalN = internal::idn,
                  std::function<ogdf::edge(ogdf::edge)> originalE = internal::ide,
                  std::function<ogdf::edge(ogdf::edge)> copyE = internal::ide) {
        for (auto n : nodes) {
            auto tn = nMap[n];
            auto on = originalN(n);
            node2ID[tn] = other.node2ID[on];
            is_dominated[tn] = other.is_dominated[on];
            is_subsumed[tn] = other.is_subsumed[on];

            // embedding breaks when inserting by edge list, so fix it
            size_t o = 0, i = 0;
            const auto &adjs = tn->adjEntries;
            for (auto adj_it = adjs.begin(); o < tn->outdeg();) {
                OGDF_ASSERT(adj_it != adjs.end());
                auto adj = *adj_it;
                ++adj_it;
                if (adj->isSource()) {
                    ++o;
                    if (i > 0) {
                        G.moveAdj(adj, ogdf::Direction::before, adjs.head());
                    }
                } else {
                    ++i;
                }
            }

#if 0
            if (!(!is_dominated[tn] || tn->indeg() == 0)) {
                {
                    auto& l = logger.lout() << "copied edges:";
                    for (auto e : edges) {
                        l<<" "<<eMap[e]<<"["<<e<<"]";
                    }
                    l<<std::endl;
                }
                {
                    auto& l = logger.lout() << "nodes' edges:";
                    for (auto adj : tn->adjEntries) {
                        l<<" "<<adj->theEdge()<<"["<<adj<<"]";
                    }
                    l<<std::endl;
                }
            }
#endif

            OGDF_ASSERT(checkNode(tn));
        }
        for (auto e : edges) {
            auto ore = other.reverse_edge[originalE(e)];
            if (ore != nullptr) {
                reverse_edge[eMap[e]] = eMap[copyE(ore)];
            } else {
                reverse_edge[eMap[e]] = nullptr;
            }
        }
    }

    void read(std::istream &is) {
        std::vector<ogdf::node> ID2node;
        read(is, ID2node);
    }
    void read(std::istream &is, std::vector<ogdf::node>& ID2node);

    void dumpBCTree();

    void safeDelete(ogdf::node n, ogdf::Graph::node_iterator &it) {
        if (*it == n) ++it; // no need to check for end, as *end == nullptr
        safeDelete(n);
    }

    void safeDelete(ogdf::node n) {
        // ID2node[node2ID[n]] = nullptr;
        G.delNode(n);
#ifdef OGDF_DEBUG
        G.consistencyCheck();
#endif
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

    // XXX not using any shortcuts here, as they often hurt in other places (like storing the global universal_in vertex)
    void addToDominatingSet(ogdf::node v, ogdf::Graph::node_iterator &it) {
        OGDF_ASSERT(!is_subsumed[v]);
        DS.push_back(node2ID[v]);
        forAllOutAdj(v,
                     [&](ogdf::adjEntry adj) {
                         markDominated(adj->twinNode());
                         return true;
                     });
        safeDelete(v, it);
    }

    void markDominated(ogdf::node v) {
        is_dominated[v] = true;
        forAllInAdj(v,
                    [&](ogdf::adjEntry adj) {
                        safeDelete(adj->theEdge());
                        return true;
                    });
    }

    void markSubsumed(ogdf::node v) {
        is_subsumed[v] = true;
        forAllOutAdj(v,
                     [&](ogdf::adjEntry adj) {
                         safeDelete(adj->theEdge());
                         return true;
                     });
    }

    bool reductionExtremeDegrees();

    bool reductionStrongSubsumption();

    bool reductionSubsumption();

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
