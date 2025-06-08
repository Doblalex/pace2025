#include "gurobi_c++.h"
#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"

void solveGurobiExactGurobi(Instance& instance) {
	GRBEnv env;
	GRBModel model(env);

	model.set(GRB_IntParam_LogToConsole, 1); // Ensure logging is enabled

	// model.set(GRB_IntParam_PoolSearchMode, 2); // Store multiple solutions
	// model.set(GRB_DoubleParam_Heuristics, 0.5); // Increase heuristic effort (optional)
	// model.set(GRB_DoubleParam_NoRelHeurTime, 10); // Allow NoRel heuristic extra time


	log << "Solving ILP witn number of nodes" << instance.G.numberOfNodes() << std::endl;

	ogdf::NodeArray<GRBVar> varmap(instance.G);
	std::vector<GRBVar> vars;

	for (auto v : instance.G.nodes) {
		if (instance.is_subsumed[v]) {
			continue;
		}
		GRBVar var;
		var = model.addVar(0, 1, 1, GRB_BINARY);
		varmap[v] = var;
		vars.push_back(var);
	}
	for (auto v : instance.G.nodes) {
		if (instance.is_dominated[v]) {
			continue;
		}
		GRBLinExpr expr;
		instance.forAllCanBeDominatedBy(v, [&](ogdf::node w) {
			expr += varmap[w];
			return true;
		});
		model.addConstr(expr >= 1);
	}

	model.optimize();

	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : instance.G.nodes) {
		if (instance.is_subsumed[v]) {
			continue;
		}
		if (varmap[v].get(GRB_DoubleAttr_X) > 0.5) {
			instance.addToDominatingSet(v);
			l << " " << instance.node2ID[v];
		}
	}
}
