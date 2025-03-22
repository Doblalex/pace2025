#include "ogdf_maxsat.hpp"

#include "EvalMaxSAT.h"

void solveEvalMaxSat(Instance& I) {
	EvalMaxSAT solver;
	const int n = I.G.numberOfNodes();
	log << "Solving MaxSat with " << n << " nodes" << std::endl;
	ogdf::NodeArray<int> varmap(I.G, -1);
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		auto var = solver.newVar();
		varmap[v] = var;
		solver.addClause({-var}, 1); // soft clause
	}
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
		std::vector<int> clause;
		if (!I.is_subsumed[v]) {
			clause.push_back(varmap[v]);
		}
		forAllInAdj(v, [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!I.is_subsumed[w]) {
				clause.push_back(varmap[w]);
			}
			return true;
		});
		solver.addClause(clause); //hard clause
	}
	solver.setTargetComputationTime(10 * 60);
	std::cout.setstate(std::ios::failbit); // https://stackoverflow.com/a/8246430
	bool solved = solver.solve();
	std::cout.clear();
	int before = I.DS.size();
	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : I.G.nodes) {
		if (!I.is_subsumed[v]) {
			if (solver.getValue(varmap[v])) {
				I.DS.push_back(I.node2ID[v]);
				l << " " << I.node2ID[n];
			}
		}
	}
	l << "\n";
	log << "Updated DS (solved " << solved << "): " << before << "+" << (I.DS.size() - before)
		<< "=" << I.DS.size() << std::endl;
}
