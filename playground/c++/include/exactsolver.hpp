#include "util.hpp"
#include "instance.hpp"

// void solveScipExactILP(Instance* instance, VertexList& dominatingSet);

void solveGurobiExactIlp(Instance* instance, VertexList& dominatingSet);

void solveEvalMaxSat(Instance* instance, VertexList& dominatingSet);