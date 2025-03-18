#include <chrono>
#include <unordered_set>
#include <boost/graph/iteration_macros.hpp>

#include "ogdf/basic/Graph_d.h"
#include "ogdf/basic/GraphSets.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/cluster/sync_plan/utils/Logging.h"
#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/fileformats/GraphIO.h"
#include "ogdf/layered/SugiyamaLayout.h"

struct Instance {
    ogdf::Graph G;
    // std::vector<ogdf::node> nodes; // TODO needed?

    std::vector<int> DS;
    ogdf::NodeArray<bool> is_dominated;
    ogdf::NodeArray<bool> can_be_DS;
    ogdf::NodeArray<int> undom_neighs;

    Instance() : is_dominated(G, false),
                 can_be_DS(G, true),
                 undom_neighs(G, -1) {
    }

    OGDF_NO_COPY(Instance)
    OGDF_NO_MOVE(Instance)

    void clear() {
        G.clear();
        // nodes.clear();
        is_dominated.init(G, false);
        can_be_DS.init(G, true);
        undom_neighs.init(G, -1);
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

        G.clear();
        std::vector<ogdf::node> nodes; // nodes.clear();
        nodes.reserve(n + 1);
        nodes.push_back(nullptr);
        for (int i = 1; i <= n; i++) {
            nodes.push_back(G.newNode(i));
        }
        for (int i = 0; i < m; i++) {
            getline(is, line);
            if (line.empty() || line[0] == 'c') continue;
            int u, v;
            std::istringstream iss(line);
            iss >> u >> v;
            G.newEdge(nodes[u], nodes[v]);
        }
        for (ogdf::node n : G.nodes) {
            undom_neighs[n] = n->degree();
        }
        OGDF_ASSERT(G.numberOfNodes() == n);
        OGDF_ASSERT(G.numberOfEdges() == m);
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
        // nodes[n->index()] = nullptr;
        G.delNode(n);
    }

    void addToDominatingSet(ogdf::node v, ogdf::Graph::node_iterator &it) {
        OGDF_ASSERT(can_be_DS[v]);
        std::unordered_set<ogdf::node> to_delete;
        to_delete.emplace(v);
        DS.push_back(v->index());
        for (auto adj : v->adjEntries) {
            auto vd = adj->twinNode();

            if (!is_dominated[v]) {
                // v is now dominating itself
                --undom_neighs[vd];
            }

            if (!is_dominated[vd]) {
                for (auto adj2 : vd->adjEntries) {
                    auto vd2 = adj2->twinNode();
                    --undom_neighs[vd2]; // vd is now dominated by v

                    // shortcut "covered" rule:
                    if (undom_neighs[vd2] == 0 && is_dominated[vd2]) {
                        to_delete.emplace(vd2);
                    }
                }
            }
            is_dominated[vd] = true;

            // shortcut "covered" rule:
            if (undom_neighs[vd] == 0) {
                to_delete.emplace(vd);
            }
        }
        for (auto n : to_delete) {
            safeDelete(n, it);
        }
    }

    bool reductionIsolatedAntennaCoveredUniversal(bool only_linear = false) {
        bool reduced = false;
        bool has_undom = false;
        ogdf::node universal = nullptr;
        int i = G.numberOfNodes();
        for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
            OGDF_ASSERT(i>=0); // protect against endless loops
            auto n = *it;
            ++it; // increment now so that we can safely delete n

            // isolated
            if (n->degree() == 0) {
                if (!is_dominated[n]) {
                    DS.push_back(n->index());
                }
                safeDelete(n, it);
                reduced = true;
                continue;
            }

            // antenna
            if (n->degree() == 1) {
                auto u = n->adjEntries.head()->twinNode();
                OGDF_ASSERT(u!=n);

                if (is_dominated[n]) {
                    if (is_dominated[u] || can_be_DS[u]) {
                        safeDelete(n, it);
                        reduced = true;
                        continue;
                    } else {
                        // we might need n to dominate u, but there might also be better choices
                        // TODO only delete if we have another neighbor?
                    }
                } else {
                    if (!only_linear) {
                        if (can_be_DS[u]) {
                            addToDominatingSet(u, it);
                            reduced = true;
                            continue;
                        } else {
                            addToDominatingSet(n, it);
                            reduced = true;
                            continue;
                        }
                    }
                }
            }

            if (is_dominated[n]) {
                // already fully covered
                if (undom_neighs[n] == 0) {
                    safeDelete(n, it);
                    reduced = true;
                    continue;
                }
                // cannot and needs not be in DS
                if (!can_be_DS[n]) {
                    safeDelete(n, it);
                    reduced = true;
                    continue;
                }
            } else {
                // no neighbor can be in DS
                if (can_be_DS[n] && !only_linear) {
                    bool must_be_DS = true;
                    for (auto adj : n->adjEntries) {
                        if (can_be_DS[adj->twinNode()]) {
                            must_be_DS = false;
                            break;
                        }
                    }
                    if (must_be_DS) {
                        addToDominatingSet(n, it);
                        reduced = true;
                        continue;
                    }
                }
            }

            // checks for universal
            if (!has_undom && !is_dominated[n]) {
                has_undom = true;
            }
            if (n->degree() == G.numberOfNodes() - 1) {
                universal = n;
            }

            // consistency checks
#ifdef OGDF_DEBUG
            int exp_undom_neighs = 0;
            bool neigh_maybe_dominating = false;
            for (auto adj : n->adjEntries) {
                if (!is_dominated[adj->twinNode()]) {
                    ++exp_undom_neighs;
                }
                if (!neigh_maybe_dominating && can_be_DS[adj->twinNode()]) {
                    neigh_maybe_dominating = true;
                }
            }
            if (!is_dominated[n] && !can_be_DS[n]) {
                OGDF_ASSERT(neigh_maybe_dominating);
            }
            OGDF_ASSERT(exp_undom_neighs == undom_neighs[n]);
#endif
        }

        // trivially done cases (no undominated or one universal to dominate them all)
        if (has_undom) {
            if (universal != nullptr) {
                DS.push_back(universal->index());
                // TODO we're done
            }
        } else {
            // TODO we're done
        }
        return reduced;
    }

    std::list<Instance> decomposeConnectedComponents() {
        ogdf::Graph::CCsInfo CC(G);
        std::list<Instance> instances;
        if (CC.numberOfCCs() > 1) {
            // instances.reserve(CC.numberOfCCs());
            for (int i = 0; i < CC.numberOfCCs(); ++i) {
                instances.emplace_back();
                Instance &I = instances.back();
                ogdf::NodeArray<ogdf::node> nMap(G, nullptr);
                ogdf::EdgeArray<ogdf::edge> eMap(G, nullptr);
                I.G.insert(CC, i, nMap, eMap); // would be easier, but we need to keep node ID. ToDo:
                // m_regNodeArrays.reserveSpace(info.numberOfNodes(cc));
                // m_regEdgeArrays.reserveSpace(info.numberOfEdges(cc));
                // auto count = I.G.insert<ogdf::Array<ogdf::node>::const_iterator, std::function<bool(ogdf::edge)>, true>(
                //     CC.nodes(i).begin(),
                //     CC.nodes(i).end(),
                //     ogdf::filter_any_edge,
                //     nMap,
                //     eMap);
                // OGDF_ASSERT(count.first == CC.numberOfNodes(i));
                // OGDF_ASSERT(count.second == CC.numberOfEdges(i));

                // I.nodes.resize(nodes.size(), nullptr);
                for (auto n : CC.nodes(i)) {
                    // OGDF_ASSERT(n->index() == nMap[n]->index());
                    // I.nodes[n->index()] = nMap[n];
                    I.is_dominated[nMap[n]] = is_dominated[n];
                    I.can_be_DS[nMap[n]] = can_be_DS[n];
                    I.undom_neighs[nMap[n]] = undom_neighs[n];
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
        std::cout << "\n\nReduce iteration " << i << " depth " << d << ": " << n << " nodes, " << m << " edges" <<
            std::endl;

        if (I.reductionIsolatedAntennaCoveredUniversal(true)) changed = true;
        // TODO reductionTwins(instance)

        std::list<Instance> comps = I.decomposeConnectedComponents();
        if (!comps.empty()) {
            I.clear(); // save some memory

            int c = 0;
            for (auto &comp : comps) {
                std::cout << "\n\nConnected component " << c << std::endl;
                ++c;

                // run remaining reduction rules on all components
                comp.reductionIsolatedAntennaCoveredUniversal(false);
                // TODO reductionDomination(instance)
                // TODO reductionDominationPaper(instance, dominatingSet))

                // and now recurse
                reduceAndSolve(comp, d + 1);
                I.DS.insert(I.DS.end(), comp.DS.begin(), comp.DS.end());
            }
            return;
        }

        if (I.reductionIsolatedAntennaCoveredUniversal(false)) changed = true;
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
    //     debug(instance->CntCanBeDominatingSet(), instance->CntNeedsDomination(), instance->n);
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
