#include "ogdf.hpp"

#include "ogdf/basic/GraphSets.h"
#include "ogdf/basic/simple_graph_alg.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/decomposition/StaticSPQRTree.h"
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
    std::string stamp = std::to_string(std::chrono::system_clock::now().time_since_epoch().count()) + "-" +
        std::to_string(G.numberOfNodes());
    log << stamp << ".svg" << std::endl;
    std::filesystem::create_directory("out");

    ogdf::BCTree BC(G);
    ogdf::GraphAttributes BCA(BC.bcTree(), ogdf::GraphAttributes::all);
    ogdf::GraphAttributes BA;
    ogdf::NodeSet<> nodes(BC.auxiliaryGraph());
    ogdf::NodeArray<ogdf::node> nMap(BC.auxiliaryGraph(), nullptr);
    ogdf::EdgeArray<ogdf::edge> eMap(BC.auxiliaryGraph(), nullptr);

    for (auto node : BC.bcTree().nodes) {
        if (BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp) {
            BCA.shape(node) = ogdf::Shape::Ellipse;
            BCA.label(node) = std::to_string(BC.numberOfNodes(node));
            BCA.width(node) = 15 + 5 * BCA.label(node).size();
        } else {
            BCA.shape(node) = ogdf::Shape::Rhomb;
            BCA.label(node) = "";
            continue;
        }

        if (BC.numberOfEdges(node) < SMALL_BLOCK) {
            continue;
        }

        nodes.clear();
        for (auto e : BC.hEdges(node)) {
            nodes.insert(e->source());
            nodes.insert(e->target());
        }
        nMap.fillWithDefault();
        eMap.fillWithDefault();
        ogdf::Graph B;
        BA.init(G, ogdf::GraphAttributes::nodeLabel);
        B.insert(nodes, BC.hEdges(node), nMap, eMap);
        for (auto n : nodes) {
            BA.label(nMap[n]) = node2ID[BC.original(n)];
        }
        auto file = "out/block-" + stamp + "-b" + std::to_string(node->index()) + "-" + std::to_string(
            BC.numberOfNodes(node));
        // ogdf::GraphIO::write(B, file + ".gml");
        ogdf::GraphIO::write(B, file + ".s6"); //

        // SPQR-Tree of B
        {
            ogdf::makeParallelFreeUndirected(B);
            // // Smooth out subdivided edges
            // ogdf::safeForEach(B.nodes,
            //                   [&B](ogdf::node n) {
            //                       if (n->degree() == 2) {
            //                           if (n->outdeg()==2 || n->indeg()==2) {
            //                               B.reverseEdge(n->adjEntries.head()->theEdge());
            //                           }
            //                           B.unsplit(n);
            //                       }
            //                   });
            // // And remove now-parallel edges
            // ogdf::makeParallelFreeUndirected(B);
            ogdf::StaticSPQRTree T(B);
            ogdf::GraphAttributes TA(T.tree(), ogdf::GraphAttributes::all);
            for (auto n : T.tree().nodes) {
                switch (T.typeOf(n)) {
                    case ogdf::SPQRTree::NodeType::SNode:
                        TA.label(n) = "S " + std::to_string(T.skeleton(n).getGraph().numberOfEdges());
                        TA.shape(n) = ogdf::Shape::Ellipse;
                        TA.strokeColor(n) = ogdf::Color(ogdf::Color::Name::Darkgreen);
                        break;
                    case ogdf::SPQRTree::NodeType::PNode:
                        TA.label(n) = "P " + std::to_string(T.skeleton(n).getGraph().numberOfEdges());
                        TA.shape(n) = ogdf::Shape::Rhomb;
                        TA.strokeColor(n) = ogdf::Color(ogdf::Color::Name::Darkblue);
                        break;
                    case ogdf::SPQRTree::NodeType::RNode:
                        TA.label(n) = "R " + std::to_string(T.skeleton(n).getGraph().numberOfEdges()) + "/" +
                            std::to_string(T.skeleton(n).getGraph().numberOfNodes());
                        TA.shape(n) = ogdf::Shape::Rect;
                        TA.strokeColor(n) = ogdf::Color(ogdf::Color::Name::Darkred);
                        break;
                }
                TA.strokeWidth(n) = 2;
                TA.width(n) = 15 + TA.label(n).size() * 5;
                TA.height(n) = 20 + TA.label(n).size() * 3;
            }
            ogdf::FMMMLayout().call(TA);
            ogdf::GraphIO::write(TA, file + "-spqr-fmmm.svg");
            ogdf::SugiyamaLayout().call(TA);
            ogdf::GraphIO::write(TA, file + "-spqr-sugi.svg");
        }
    }

    ogdf::FMMMLayout().call(BCA);
    ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".svg");
    // ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".gml");
    ogdf::SugiyamaLayout().call(BCA);
    ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".svg");
    // ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".gml");
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
    ogdf::NodeArray<bool> replaced(BC.bcTree(), false);
    ogdf::NodeArray<int> subtree_size(BC.bcTree(), 0); //
    {
        std::vector<ogdf::node> todo;
        ogdf::NodeArray<int> subtrees_seen(BC.bcTree(), 0);
        for (auto n : BC.bcTree().nodes) {
            OGDF_ASSERT(n->outdeg() <= 1);
            if (n->indeg() == 0) {
                todo.push_back(n);
            }
            if (BC.typeOfBNode(n) == ogdf::BCTree::BNodeType::BComp) {
                subtree_size[n] = BC.numberOfNodes(n);
            } else {
                subtree_size[n] = -n->degree() + 1; // discount copies of C-nodes
            }
        }
        while (!todo.empty()) {
            auto n = todo.back();
            todo.pop_back();
            ogdf::node p = BC.parent(n);

            subtree_size[p] += subtree_size[n];
            subtrees_seen[p] += 1;
            if (subtrees_seen[p] == p->indeg()) {
                todo.push_back(p);
            }
        }
    }

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
            BC.numberOfNodes(node) < SMALL_BLOCK &&
            !replaced[node]
        ) {
            ogdf::node h_cv;
            ogdf::node parent = BC.parent(node);
            bool is_root = false;
            if (parent == nullptr) {
                // deg-1 root
                h_cv = BC.hParNode(node->adjEntries.head()->twinNode());
                parent = node->adjEntries.head()->twinNode();
                is_root = true;
            } else {
                h_cv = BC.hRefNode(node);
            }
            OGDF_ASSERT(parent!=nullptr);
            OGDF_ASSERT(h_cv != nullptr);
            OGDF_ASSERT(nodes.isMember(h_cv));
            ogdf::node cv = BC.original(h_cv);

            if (replaced[parent]) continue;
            log << "Processing leaf block with " << BC.numberOfNodes(node) << " nodes." << std::endl;

            nodes.clear();
            for (auto e : BC.hEdges(node)) {
                nodes.insert(e->source());
                nodes.insert(e->target());
            }

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
            replaced[node] = true;
            changed = true;
            if (I1.DS.size() == I2.DS.size()) {
                replaced[parent] = true;
                int old_size = subtree_size[parent];
                for (auto p = BC.parent(parent); p != nullptr && !is_root; p = BC.parent(p)) {
                    subtree_size[p] -= old_size;
                }

                forAllOutAdj(cv,
                             [&](ogdf::adjEntry adj) {
                                 markDominated(adj->twinNode());
                                 return true;
                             });
                for (auto n : nodes) {
                    safeDelete(BC.original(n));
                }
                DS.insert(DS.end(), I2.DS.begin(), I2.DS.end());
            } else {
                // TODO otherwise replace by gadget so that we don't reprocess the same block over all the time
                int size_diff = subtree_size[node] - 5; // TODO gadget size
                for (auto p = parent; p != nullptr && !is_root; p = BC.parent(p)) {
                    subtree_size[p] -= size_diff;
                }
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

    // if (I.G.numberOfNodes() > SMALL_BLOCK) {
    log << "SKIPPING!" << std::endl;
    return;
    // }

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
