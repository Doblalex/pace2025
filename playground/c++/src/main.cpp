#include "util.hpp"
#include "cxxopts.h"
#include "readinstance.hpp"
#include "preprocessing.hpp"
#include "heuristics.hpp"
#include "exactsolver.hpp"

void recursiveReduction(Instance* instance, VertexList& dominatingSet) {
    // debug(instance->n, instance->m);
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

    
    
    debug(instance->CntCanBeDominatingSet(), instance->CntNeedsDomination(), instance->n);
    solveCPSat(instance, dominatingSet);
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
    reductionTwins(instance); // It seems very unsafe to do twin reduction rule once there are dominated vertices. 
    recursiveReduction(instance, dominatingSet);    
    debug(dominatingSet.size());
    cout<<dominatingSet.size()<<endl;
}