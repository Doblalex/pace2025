#ifndef PREPROCESSING_H
#define PREPROCESSING_H


#include "util.hpp"
#include "instance.hpp"
#include "partition.hpp"

bool reductionIsolated(Instance* instance, VertexList& dominatingSet);

bool reduceUniversal(Instance* instance, VertexList& dominatingSet);

bool reductionCovered(Instance* instance);

bool reductionDegree1(Instance* instance, VertexList& dominatingSet);

bool reductionDomination(Instance* instance);

bool reductionDominationPaper(Instance* instance, VertexList& dominatingSet);

bool reductionTwins(Instance* instance);

list<Instance*> decomposeConnectedComponents(Instance* instance);

bool reductionByCanBeDominatingSet(Instance* instance, VertexList& dominatingSet);

#endif