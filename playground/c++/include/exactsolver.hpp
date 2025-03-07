#include "util.hpp"
#include "instance.hpp"

// void solveScipExactILP(Instance* instance, VertexList& dominatingSet);

#ifdef USE_GUROBI
void solveGurobiExactIlp(Instance* instance, VertexList& dominatingSet);
#endif

void solveEvalMaxSat(Instance* instance, VertexList& dominatingSet);

#ifdef USE_ORTOOLS
void solveCPSat(Instance* instance, VertexList& dominatingSet);
#endif