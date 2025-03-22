#include "util.hpp"
#include "cxxopts.h"
#include "readinstance.hpp"
#include "preprocessing.hpp"
#include "heuristics.hpp"
#include "exactsolver.hpp"

void recursiveReduction(Instance* instance, VertexList& dominatingSet) {
    // debug(instance->n, instance->m);

    // BGL_FORALL_VERTICES(v, *instance->G, Graph) {
    //     if (!(*instance->G)[v].is_dominated) {
    //         bool ok = (*instance->G)[v].can_be_dominating_set;
    //         BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
    //             if ((*instance->G)[w].can_be_dominating_set) {
    //                 ok = true;
    //             }
    //         }
    //         if (!ok) {
    //             debug("something went wrong");
    //             exit(1);
    //         }
    //     }
    // }

    // // TODO: reduction when there is only one neighbor that can dominate me

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
