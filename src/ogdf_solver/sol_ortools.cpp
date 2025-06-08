#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"
#include "ortools/linear_solver/linear_expr.h"
#include "ortools/linear_solver/linear_solver.h"

void solvecpsat(Instance& I) {
	using namespace operations_research;
	std::unique_ptr<MPSolver> solver(MPSolver::CreateSolver("CP-SAT"));
	int before = I.DS.size();
	log << "Solving CP-SAT with " << I.G.numberOfNodes() << " nodes" << std::endl;
	if (!solver) {
		std::cerr << "solver not available" << std::endl;
		exit(1);
	}
	ogdf::NodeArray<MPVariable*> varmap(I.G, nullptr);
	LinearExpr obj;
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		auto var = solver->MakeBoolVar("");
		varmap[v] = var;
		obj += var;
	}
#ifdef SAT_CACHE
	std::vector<std::vector<int>> hclauses;
	hclauses.reserve(I.G.numberOfNodes());
#endif
	MPObjective* const objective = solver->MutableObjective();
	objective->MinimizeLinearExpr(obj);
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
		LinearExpr expr;
#ifdef SAT_CACHE
		hclauses.emplace_back();
		hclauses.back().reserve(v->indeg() + 1);
#endif
		if (!I.is_subsumed[v]) {
			expr += varmap[v];
#ifdef SAT_CACHE
			hclauses.back().push_back(I.node2ID[v]);
#endif
		}
		auto add_neigh = [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!adj->isSource() && !I.is_subsumed[w]) {
				expr += varmap[w];
#ifdef SAT_CACHE
				hclauses.back().push_back(I.node2ID[w]);
#endif
			}
			return true;
		};
		forAllInAdj(v, add_neigh);
		ogdf::safeForEach(I.hidden_edges.adjEntries(v), add_neigh);
		solver->MakeRowConstraint(expr >= 1); //hard clause
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

	// solver->EnableOutput();
	const MPSolver::ResultStatus result_status = solver->Solve();
	if (result_status != MPSolver::OPTIMAL) {
		std::cerr << "result_status " << result_status << std::endl;
		throw std::runtime_error("MPSolver didn't find optimal result!");
	}

	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : I.G.nodes) {
		if (!I.is_subsumed[v]) {
			if (varmap[v]->solution_value() > 0.5) {
				I.DS.insert(I.node2ID[v]);
				l << " " << I.node2ID[v];
#ifdef SAT_CACHE
				f << " " << I.node2ID[v];
#endif
			}
		}
	}
	l << "\n";
	log << "Updated DS (solved " << result_status << "): " << before << "+"
		<< (I.DS.size() - before) << "=" << I.DS.size() << std::endl;
	solver->Clear();
}
