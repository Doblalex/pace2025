#include "util.hpp"
#include "cxxopts.h"
#include "readinstance.hpp"
#include "preprocessing.hpp"
#include "heuristics.hpp"
#include "exactsolver.hpp"

void recursiveReduction(Instance* instance, VertexList& dominatingSet) {
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

    // TODO: bugs in the reduction below
    // if (reductionDegree1(instance, dominatingSet)) {
    //     recursiveReduction(instance, dominatingSet);
    //     return;
    // }
    
    if (reductionDomination(instance)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    if (reductionByCanBeDominatingSet(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    // TODO: bugs in the reduction below
    // if (reductionDominationPaper(instance, dominatingSet)) {
    //     recursiveReduction(instance, dominatingSet);
    //     return;
    // }  

    
    
    debug(instance->CntCanBeDominatingSet(), instance->CntNeedsDomination(), instance->n);
    solveEvalMaxSat(instance, dominatingSet);
    delete instance;
}

int main(int argc, char** argv) {
    globalprops* props = new globalprops();

    auto instance = read_instance(props);
    // VertexList greedysol;
    // Greedy(instance, greedysol);
    // debug(greedysol.size());
    
    debug(instance->n, instance->m);
    VertexList dominatingSet;
    recursiveReduction(instance, dominatingSet);    
    debug(dominatingSet.size());
    cout<<dominatingSet.size()<<endl;
}