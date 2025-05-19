#include "EvalMaxSAT.h"
#include "matching.hpp"
#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"

void createSolver(Instance& I, EvalMaxSAT<Solver_cadical>* solver, ogdf::NodeArray<int>& varmap,
		std::vector<std::vector<int>>& hclauses,
		std::vector<std::pair<ogdf::node, ogdf::node>>& size_2_clauses) {
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		int var = solver->newSoftVar(true, -1);
		varmap[v] = var;
	}

#ifdef SAT_CACHE
	hclauses.reserve(I.G.numberOfNodes());
#endif
	std::vector<int> clause;
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
		std::vector<ogdf::node> nodes_in_clause;
		clause.clear();
		// clause.reserve(v->indeg() + 1);
#ifdef SAT_CACHE
		hclauses.emplace_back();
		hclauses.back().reserve(v->indeg() + 1);
#endif
		if ((!I.is_subsumed[v] && !I.is_dominated[v]) || I.is_hidden_loop[v]) {
			clause.push_back(varmap[v]);
			nodes_in_clause.push_back(v);
#ifdef SAT_CACHE
			hclauses.back().push_back(I.node2ID[v]);
#endif
		}
		auto add_neigh = [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!adj->isSource() && !I.is_subsumed[w]) {
				clause.push_back(varmap[w]);
				nodes_in_clause.push_back(w);
#ifdef SAT_CACHE
				hclauses.back().push_back(I.node2ID[w]);
#endif
			}
			return true;
		};
		forAllInAdj(v, add_neigh);
		if (clause.empty()) {
			continue;
		}
		solver->addClause(clause); //hard clause

		if (clause.size() == 2) {
			auto uu = nodes_in_clause[0];
			auto vv = nodes_in_clause[1];
			size_2_clauses.push_back({uu, vv});
		}

#ifdef SAT_CACHE
		std::sort(hclauses.back().begin(), hclauses.back().end());
#endif
	}
}

void solveEvalMaxSat(Instance& I) {
	EvalMaxSAT<Solver_cadical>* solver = new EvalMaxSAT();
	int before = I.DS.size();
	log << "Solving EvalMaxSat with " << I.G.numberOfNodes() << " nodes" << std::endl;
	ogdf::NodeArray<int> varmap(I.G, -1);
	ogdf::NodeArray<int> mapBlossom(I.G, -1);
	std::vector<ogdf::node> mapBlossomInv;
	std::vector<std::vector<int>> hclauses;
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		mapBlossomInv.push_back(v);
		mapBlossom[v] = mapBlossomInv.size() - 1;
	}
	std::vector<std::vector<int>> blossomadj(mapBlossomInv.size());
	std::vector<std::pair<ogdf::node, ogdf::node>> size_2_clauses;

	createSolver(I, solver, varmap, hclauses, size_2_clauses);

#ifdef SAT_CACHE
	std::string filename;
	if (try_load_solution(I, hclauses, filename)) {
		return;
	}
	std::ofstream f(filename);
#endif
	for (auto& p : size_2_clauses) {
		auto u = mapBlossom[p.first];
		auto v = mapBlossom[p.second];
		blossomadj[u].push_back(v);
		blossomadj[v].push_back(u);
	}
	Blossom matching_blossom(blossomadj);
	auto ans = matching_blossom.solve();
	auto matching_size = ans.second;

#ifndef PACE_LOG
	std::cout.setstate(std::ios::failbit); // https://stackoverflow.com/a/8246430
#endif
	solver->adapt_am1_exact();
	solver->adapt_am1_FastHeuristicV7();

	{
		// Try if matching plus heuristic gives better lower bound
		EvalMaxSAT<Solver_cadical>* solver2 = new EvalMaxSAT();
		hclauses.clear();
		size_2_clauses.clear();
		ogdf::NodeArray<int> varmap2(I.G, -1);
		createSolver(I, solver2, varmap2, hclauses, size_2_clauses);
		for (int i = 0; i < ans.first.size(); i++) {
			if (ans.first[i] != -1) {
				auto u = mapBlossomInv[i];
				auto v = mapBlossomInv[ans.first[i]];
				if (u > v) {
					continue;
				}
				auto varu = varmap2[u];
				auto varv = varmap2[v];
				solver2->processAtMostOne({-varu, -varv});
			}
		}
		solver2->adapt_am1_exact();
		solver2->adapt_am1_FastHeuristicV7();
		if (solver2->cost >= solver->cost) {
			log << "Matching found better or at least as good lower bound" << std::endl;
			log << "Matching size: " << matching_size << std::endl;
			log << "Solver2 cost: " << solver2->cost << std::endl;
			log << "Solver cost: " << solver->cost << std::endl;
			delete solver;
			varmap = varmap2;
			solver = solver2;
		} else {
			delete solver2;
		}
	}
	// set parameters as default divided by divide_by, otherwise each core finding takes too long
	double divide_by =
			2; // as evalmaxsat is made for 1 hour, this seems reasonable as we have 30 minutes.
	solver->setTargetComputationTime(30 * 60); // this seems to work over multiple solvers
	solver->setBoundRefTime(5 / divide_by, (5 * 60) / divide_by);
	solver->setCoef(10, 1.66);
	bool solved = solver->solve();
#ifndef PACE_LOG
	std::cout.clear();
#endif
	if (!solved) {
		throw std::runtime_error("EvalMaxSAT didn't find optimal result!");
	}

	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : I.G.nodes) {
		if (!I.is_subsumed[v]) {
			if (solver->getValue(varmap[v])) {
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
	delete solver;
}
