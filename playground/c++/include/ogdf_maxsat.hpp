#pragma once

#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"

void solveEvalMaxSat(Instance& I);
#define USE_ORTOOLS
#ifdef USE_ORTOOLS
void solvecpsat(Instance& I);
#endif