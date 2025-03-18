#include "util.hpp"
#include "cxxopts.h"
#include "readinstance.hpp"
#include "preprocessing.hpp"
#include "heuristics.hpp"
#include "exactsolver.hpp"

#ifdef USE_OGDF
#include "ogdf/basic/Graph_d.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/fileformats/GraphIO.h"
#include "ogdf/layered/SugiyamaLayout.h"
#endif

void dumpBCTree(Instance *instance, VertexList &dominatingSet) {
#ifdef USE_OGDF
    ogdf::Graph G;
    std::vector<ogdf::node> map;
    boost::property_map<Graph, boost::vertex_index_t>::type index_map =
        boost::get(boost::vertex_index, (*instance->G));
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        auto idxv = (*instance->G)[v].id;
        if (idxv >= map.size()) {
            map.resize(idxv + 1, nullptr);
        }
        if (map[idxv] == nullptr) {
            map[idxv] = G.newNode(idxv);
        }
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            auto idxw = (*instance->G)[w].id;
            if (idxw >= map.size()) {
                map.resize(idxw + 1, nullptr);
            }
            if (map[idxw] == nullptr) {
                map[idxw] = G.newNode(idxw);
            }
            if (idxv < idxw) {
                G.newEdge(map[idxv], map[idxw]);
            }
        }
    }
    OGDF_ASSERT(G.numberOfNodes() == instance->n);
    OGDF_ASSERT(G.numberOfEdges() == instance->m);

    std::vector<ogdf::node> deg1;
    for (auto n : G.nodes) {
        if (n->degree() <= 1) {
            deg1.push_back(n);
        }
    }
    for (auto n : deg1) {
        G.delNode(n);
    }

    ogdf::BCTree BC(G);
    ogdf::GraphAttributes BCA(BC.bcTree(), ogdf::GraphAttributes::all);
    for (auto node : BC.bcTree().nodes) {
        if (BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp) {
            BCA.shape(node) = ogdf::Shape::Ellipse;
            BCA.label(node) = to_string(BC.numberOfNodes(node));
            BCA.width(node) = 15 + 5 * BCA.label(node).size();
        } else {
            BCA.shape(node) = ogdf::Shape::Rhomb;
            BCA.label(node) = "";
        }
    }
    ogdf::FMMMLayout().call(BCA);
    ogdf::GraphIO::write(BCA, "fmmm.svg");
    ogdf::GraphIO::write(BCA, "fmmm.gml");
    ogdf::SugiyamaLayout().call(BCA);
    ogdf::GraphIO::write(BCA, "sugi.svg");
    ogdf::GraphIO::write(BCA, "sugi.gml");
    cout << "dumped" << endl;
#endif
}

void recursiveReduction(Instance* instance, VertexList& dominatingSet) {
    // debug(instance->n, instance->m);

    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if (!(*instance->G)[v].is_dominated) {
            bool ok = (*instance->G)[v].can_be_dominating_set;
            BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
                if ((*instance->G)[w].can_be_dominating_set) {
                    ok = true;
                }
            }
            if (!ok) {
                debug("something went wrong");
                exit(1);
            }
        }
    }

    // TODO: reduction when there is only one neighbor that can dominate me

    if (instance->n == 0) {
        return;
    }

    if (reductionIsolated(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }
    if (reduceUniversal(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }
    if (reductionCovered(instance)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    if (reductionTwins(instance)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    auto subinstances = decomposeConnectedComponents(instance);
    if (subinstances.size() > 1) {
        for (auto subinstance: subinstances) {
            recursiveReduction(subinstance, dominatingSet);
        }
        return;
    }

    instance = subinstances.front();

    // non-linear reductions
    if (reductionDegree1(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    if (reductionDomination(instance)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    if (reductionByCanBeDominatingSet(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }
    if (reductionDominationPaper(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }


    // dumpBCTree(instance, dominatingSet);
    debug(instance->CntCanBeDominatingSet(), instance->CntNeedsDomination(), instance->n);
    #ifdef USE_ORTOOLS
    solveCPSat(instance, dominatingSet);
    #else
    solveEvalMaxSat(instance, dominatingSet);
    #endif
    delete instance;
}

int main(int argc, char** argv) {
    globalprops* props = new globalprops();
    auto start = std::chrono::system_clock::now();
    auto instance = read_instance(props);
    // VertexList greedysol;
    // Greedy(instance, greedysol);
    // debug(greedysol.size());

    debug(instance->n, instance->m);
    VertexList dominatingSet;
    recursiveReduction(instance, dominatingSet);
    debug(dominatingSet.size());

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> runtime = end-start;
    #ifndef MYLOCAL
    cout<<"Solution size: "<<dominatingSet.size()<<endl;
    cout<<"Solution time: "<<runtime.count()<<endl;
    // for (auto v: dominatingSet) {
    //     cout << v << endl;
    // }
    #endif
}
