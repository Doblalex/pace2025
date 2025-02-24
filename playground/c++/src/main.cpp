#include "util.hpp"
#include "cxxopts.h"
#include "readinstance.hpp"
#include "preprocessing.hpp"

void recursiveReduction(Instance* instance, VertexList& dominatingSet) {
    if (instance->n == 0) {
        return;
    }
    if (reductionIsolated(instance, dominatingSet)) {
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
    if (reductionDomination(instance)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }

    if (reductionByCanBeDominatingSet(instance, dominatingSet)) {
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

    if (reductionDegree1(instance, dominatingSet)) {
        recursiveReduction(instance, dominatingSet);
        return;
    }
    
    debug(instance->CntCanBeDominatingSet(), instance->CntNeedsDomination(), instance->n);
    delete instance;
}

int main(int argc, char** argv) {
    globalprops* props = new globalprops();

    auto instance = read_instance(props);
    debug(instance->n, instance->m);
    VertexList dominatingSet;
    recursiveReduction(instance, dominatingSet);    
    debug(dominatingSet.size());
}