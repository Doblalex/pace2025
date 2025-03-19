#include "ogdf.hpp"

#include "ogdf/basic/GraphSets.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/fileformats/GraphIO.h"
#include "ogdf/layered/SugiyamaLayout.h"

#include <chrono>
#include <filesystem>

#define SMALL_BLOCK 100

ogdf::Logger logger;

ogdf::node internal::idn(ogdf::node n) { return n; }
ogdf::edge internal::ide(ogdf::edge n) { return n; }

void Instance::read(std::istream &is) {
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

void Instance::dumpBCTree() {
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
    log << stamp << ".svg" << std::endl;
    ogdf::FMMMLayout().call(BCA);
    std::filesystem::create_directory("out");
    ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".svg");
    ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".gml");
    ogdf::SugiyamaLayout().call(BCA);
    ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".svg");
    ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".gml");
}

bool Instance::reductionExtremeDegrees() {
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
        OGDF_ASSERT(n->outdeg() == 0 || n->adjEntries.head()->isSource());
        OGDF_ASSERT(n->indeg() == 0 || !n->adjEntries.tail()->isSource());

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
        if (has_undom) {
            safeDelete(universal_in);
        } else {
            addToDominatingSet(universal_in);
        }
        reduced = true;
    }

    return reduced;
}

std::list<Instance> Instance::decomposeConnectedComponents() {
    ogdf::Graph::CCsInfo CC(G);
    std::list<Instance> instances;
    if (CC.numberOfCCs() > 1) {
        ogdf::NodeArray<ogdf::node> nMap(G, nullptr); // can safely be reused for different CCs
        ogdf::EdgeArray<ogdf::edge> eMap(G, nullptr);
        for (int i = 0; i < CC.numberOfCCs(); ++i) {
            instances.emplace_back();
            Instance &I = instances.back();
            I.G.insert(CC, i, nMap, eMap);
            I.initFrom(*this, CC.nodes(i), CC.edges(i), nMap, eMap);
        }
    }
    return instances;
}

bool Instance::reductionBCTree() {
    ogdf::BCTree BC(G);
    ogdf::NodeSet<> nodes(BC.auxiliaryGraph());
    ogdf::NodeArray<ogdf::node> nMap1(BC.auxiliaryGraph(), nullptr);
    ogdf::EdgeArray<ogdf::edge> eMap1(BC.auxiliaryGraph(), nullptr);
    ogdf::NodeArray<ogdf::node> nMap2(BC.auxiliaryGraph(), nullptr);
    ogdf::EdgeArray<ogdf::edge> eMap2(BC.auxiliaryGraph(), nullptr);
    bool changed = false;
    for (auto node : BC.bcTree().nodes) {
        if (
            node->degree() == 1 && // TODO might also work with larger subtrees that still have few nodes
            BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp &&
            BC.numberOfNodes(node) < SMALL_BLOCK
        ) {
            log << "Processing leaf block with " << BC.numberOfNodes(node) << " nodes." << std::endl;

            nodes.clear();
            for (auto e : BC.hEdges(node)) {
                nodes.insert(e->source());
                nodes.insert(e->target());
            }

            ogdf::node h_cv;
            if (node->outdeg() > 0) {
                // deg-1 root
                h_cv = BC.hParNode(node->adjEntries.head()->twinNode());
            } else {
                h_cv = BC.hRefNode(node);
            }
            OGDF_ASSERT(h_cv != nullptr);
            OGDF_ASSERT(nodes.isMember(h_cv));
            ogdf::node cv = BC.original(h_cv);

            Instance I1;
            nMap1.fillWithDefault();
            eMap1.fillWithDefault();
            I1.G.insert(nodes, BC.hEdges(node), nMap1, eMap1);
            I1.initFrom(
                *this,
                nodes,
                BC.hEdges(node),
                nMap1,
                eMap1,
                [&BC](ogdf::node n) { return BC.original(n); },
                [&BC](ogdf::edge n) { return BC.original(n); },
                [&BC](ogdf::edge n) { return BC.rep(n); }
            );

            Instance I2;
            nMap2.fillWithDefault();
            eMap2.fillWithDefault();
            I2.G.insert(nodes, BC.hEdges(node), nMap2, eMap2);
            I2.initFrom(
                *this,
                nodes,
                BC.hEdges(node),
                nMap2,
                eMap2,
                [&BC](ogdf::node n) { return BC.original(n); },
                [&BC](ogdf::edge n) { return BC.original(n); },
                [&BC](ogdf::edge n) { return BC.rep(n); });

            log << "Solving I1 with dominated CV." << std::endl; {
                ogdf::Logger::Indent _(logger);
                I1.markDominated(nMap1[h_cv]);
                reduceAndSolve(I1, 100);
            }

            log << "Solving I2 with CV added to DS." << std::endl; {
                ogdf::Logger::Indent _(logger);
                I2.addToDominatingSet(nMap2[h_cv]);
                reduceAndSolve(I2, 200);
            }

            log << "I1 DS: " << I1.DS.size() << " I2 DS: " << I2.DS.size() << std::endl;
            if (I1.DS.size() == I2.DS.size()) {
                forAllOutAdj(cv,
                             [&](ogdf::adjEntry adj) {
                                 markDominated(adj->twinNode());
                                 return true;
                             });
                for (auto n : nodes) {
                    safeDelete(BC.original(n));
                }
                DS.insert(DS.end(), I2.DS.begin(), I2.DS.end());
                changed = true;
                break; // TODO we could continue as long as we don't process any block attached to the same CV
            } else {
                // TODO otherwise replace by gadget so that we don't reprocess the same block over all the time
            }
        }
    }
    return changed;
}

void reduceAndSolve(Instance &I, int d) {
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

                // TODO run remaining reduction rules on all components?
                // TODO reductionDomination(instance)
                // TODO reductionDominationPaper(instance, dominatingSet))

                // and now recurse
                reduceAndSolve(comp, d + 1);
                log << "DS before adding connected component " << c << ": " << I.DS.size() << std::endl;
                log << "DS of connected component " << c << ": " << comp.DS.size() << std::endl;
                I.DS.insert(I.DS.end(), comp.DS.begin(), comp.DS.end());
                log << "DS after adding connected component " << c << ": " << I.DS.size() << std::endl;
            }
            return;
        }

        // if (I.reductionBCTree()) changed = true;
        // TODO reductionDomination(instance)
        // TODO reductionDominationPaper(instance, dominatingSet))

        OGDF_ASSERT(!changed || I.G.numberOfNodes() < n|| I.G.numberOfEdges() < m);
        ++i;
    }

    if (I.G.numberOfNodes() < 1) {
        log << "Reduced instance is empty!" << std::endl;
        return;
    }
    auto [can,need] = I.dominationStats();
    log << "Reduced instance contains " << n << " nodes, " << m << " edges. "
        << need << " vertices need to be dominated, " << can << " are eligible for the DS." << std::endl;
    I.dumpBCTree();

    if (I.G.numberOfNodes() > SMALL_BLOCK) {
        log << "SKIPPING!" << std::endl;
        return;
    }

    // now to solving...
#ifdef USE_ORTOOLS
        solveCPSat(instance, dominatingSet);
#else
    solveEvalMaxSat(I);
#endif
}

int main(int argc, char **argv) {
    auto start = std::chrono::system_clock::now();
    Instance I;
    I.read(std::cin);
    reduceAndSolve(I);
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> runtime = end - start;
}
