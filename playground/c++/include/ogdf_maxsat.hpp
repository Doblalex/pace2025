#pragma once

#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"

// void solveGurobiExactGurobi(Instance& instance);

void solveEvalMaxSat(Instance& I);

#ifdef USE_ORTOOLS
void solvecpsat(Instance& I);
#endif