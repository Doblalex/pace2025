#include "EvalMaxSAT.h"
#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"

void solveEvalMaxSat(Instance& I) {
	EvalMaxSAT solver;
	int before = I.DS.size();
	log << "Solving EvalMaxSat with " << I.G.numberOfNodes() << " nodes" << std::endl;
	ogdf::NodeArray<int> varmap(I.G, -1);
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		int var = solver.newSoftVar(true, -1);
		varmap[v] = var;
	}
#ifdef SAT_CACHE
	std::vector<std::vector<int>> hclauses;
	hclauses.reserve(I.G.numberOfNodes());
#endif
	std::vector<int> clause;
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) { // might have hidden edges
			continue;
		}
		clause.clear();
		// clause.reserve(v->indeg() + 1);
#ifdef SAT_CACHE
		hclauses.emplace_back();
		hclauses.back().reserve(v->indeg() + 1);
#endif
		if ((!I.is_subsumed[v] && !I.is_dominated[v]) || I.is_hidden_loop[v]) {
			clause.push_back(varmap[v]);
#ifdef SAT_CACHE
			hclauses.back().push_back(I.node2ID[v]);
#endif
		}
		auto add_neigh = [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!adj->isSource() && !I.is_subsumed[w]) {
				clause.push_back(varmap[w]);
#ifdef SAT_CACHE
				hclauses.back().push_back(I.node2ID[w]);
#endif
			}
			return true;
		};
		forAllInAdj(v, add_neigh);
		ogdf::safeForEach(I.hidden_edges.adjEntries(v), add_neigh);
		if (clause.empty()) {
			continue;
		}
		solver.addClause(clause); //hard clause
#ifdef SAT_CACHE
		std::sort(hclauses.back().begin(), hclauses.back().end());
#endif
	}
#ifdef SAT_CACHE
	std::string filename;
	if (try_load_solution(I, hclauses, filename)) {
		return;
	}
	std::ofstream f(filename);
#endif

	solver.setTargetComputationTime(20 * 60);
	solver.setBoundRefTime(0.1, 10);
	solver.setCoef(1, 0.1);
#ifndef PACE_LOG
	std::cout.setstate(std::ios::failbit); // https://stackoverflow.com/a/8246430
#endif
	// solver.cost = 2800;
	// solver.cost = 800;
	bool solved = solver.solve();
#ifndef PACE_LOG
	std::cout.clear();
#endif
	if (!solved) {
		throw std::runtime_error("EvalMaxSAT didn't find optimal result!");
	}

	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : I.G.nodes) {
		if (!I.is_subsumed[v]) {
			if (solver.getValue(varmap[v])) {
				I.DS.insert(I.node2ID[v]);
				l << " " << I.node2ID[v];
#ifdef SAT_CACHE
				f << " " << I.node2ID[v];
#endif
			}
		}
	}
	l << "\n";
	log << "Updated DS (solved " << solved << "): " << before << "+" << (I.DS.size() - before)
		<< "=" << I.DS.size() << std::endl;
}
