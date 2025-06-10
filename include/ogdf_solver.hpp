#pragma once

#include "ogdf_instance.hpp"

void reduceAndSolve(Instance& I, int d = 0);

void solveGreedy(Instance& I);

#ifndef PACE_EMS_FACTOR
#define PACE_EMS_FACTOR 1
#endif

#ifdef USE_EVALMAXSAT
void solveEvalMaxSat(Instance& I);
#endif

#ifdef USE_ORTOOLS
void solvecpsat(Instance& I);
#endif

#ifdef USE_UWRMAXSAT
void solveIPAMIR(Instance& I);
#endif

#ifdef USE_GUROBI
void solveGurobiExactGurobi(Instance& instance);
#endif

#ifdef SAT_CACHE
uint64_t hash_clauses(const std::vector<std::vector<int>>& clauses);
void dump_sat(const std::string& file, const std::vector<std::vector<int>>& clauses);
bool is_same_sat(const std::string& file, const std::vector<std::vector<int>>& clauses);
std::string get_filename(uint64_t hash, const std::string& ext, const std::string& dir = "cache/");
bool try_load_solution(Instance& I, std::vector<std::vector<int>>& hclauses, std::string& filename);
#endif
