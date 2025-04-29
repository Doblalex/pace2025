#include "ogdf_maxsat.hpp"
// #include "gurobi_c++.h"
#include "EvalMaxSAT.h"
#include "ogdf_util.hpp"

#ifdef USE_ORTOOLS
#	include "ortools/linear_solver/linear_expr.h"
#	include "ortools/linear_solver/linear_solver.h"
#endif



// void solveGurobiExactGurobi(Instance& instance) {
// 	GRBEnv env;
// 	GRBModel model(env);

// 	model.set(GRB_IntParam_LogToConsole, 1); // Ensure logging is enabled

// 	// model.set(GRB_IntParam_PoolSearchMode, 2); // Store multiple solutions
// 	// model.set(GRB_DoubleParam_Heuristics, 0.5); // Increase heuristic effort (optional)
// 	// model.set(GRB_DoubleParam_NoRelHeurTime, 10); // Allow NoRel heuristic extra time


// 	log<<"Solving ILP witn number of nodes"<< instance.G.numberOfNodes()<<std::endl;

// 	ogdf::NodeArray<GRBVar> varmap(instance.G);
// 	std::vector<GRBVar> vars;

// 	for (auto v: instance.G.nodes) {
// 		if (instance.is_subsumed[v]) {
// 			continue;
// 		}
// 		GRBVar var;
// 		var = model.addVar(0, 1, 1, GRB_BINARY);
// 		varmap[v] = var;
// 		vars.push_back(var);
// 	}
// 	for (auto v: instance.G.nodes) {
// 		if (instance.is_dominated[v]) {
// 			continue;
// 		}
// 		GRBLinExpr expr;
// 		instance.forAllCanBeDominatedBy(v, [&](ogdf::node w) {
// 			expr += varmap[w];
// 			return true;
// 		});
// 		model.addConstr(expr >= 1);
// 	}

// 	model.optimize();

// 	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
// 	for (auto v: instance.G.nodes) {
// 		if (instance.is_subsumed[v]) {
// 			continue;
// 		}
// 		if (varmap[v].get(GRB_DoubleAttr_X) > 0.5) {
// 			instance.addToDominatingSet(v);
// 			l << " " << instance.node2ID[v];
// 		}
// 	}
// }

#ifdef EMS_CACHE
uint64_t hash_clauses(const std::vector<std::vector<int>>& clauses) {
	uint64_t hash = FNV1a_64_SEED;
	OGDF_ASSERT(clauses.size() > 0);
	FNV1a_64_update(hash, clauses.size());
	for (auto& clause : clauses) {
		OGDF_ASSERT(clause.size() > 0);
		FNV1a_64_update(hash, clause.size());
		for (auto v : clause) {
			FNV1a_64_update(hash, v);
		}
	}
	return hash;
}

void dump_sat(const std::string& file, const std::vector<std::vector<int>>& clauses) {
	std::ofstream f(file);
	for (auto& clause : clauses) {
		bool first = true;
		for (auto x : clause) {
			if (first) {
				first = false;
			} else {
				f << " ";
			}
			f << x;
		}
		f << "\n";
	}
}

bool is_same_sat(const std::string& file, const std::vector<std::vector<int>>& clauses) {
	std::ifstream f(file);
	int a;
	for (auto& clause : clauses) {
		for (auto x : clause) {
			if (!(f >> a) || a != x) {
				return false;
			}
		}
	}
	if (f >> a) {
		return false;
	} else {
		return true;
	}
}

std::string get_filename(uint64_t hash, const std::string& ext, const std::string& dir = "cache/") {
	std::stringstream fnstr;
	fnstr << dir;
	fnstr << std::setfill('0') << std::setw(sizeof(uint64_t) * 2) << std::hex << hash;
	fnstr << ext;
	return fnstr.str();
}

bool try_load_solution(Instance& I, std::vector<std::vector<int>>& hclauses, std::string& filename) {
	std::sort(hclauses.begin(), hclauses.end());
	uint64_t hash = hash_clauses(hclauses);
	filename = get_filename(hash, ".sol");
	std::string filename_sat = get_filename(hash, ".sat");

	int before = I.DS.size();
	if (std::filesystem::exists(filename)) {
		log << "Found cached solution " << filename << std::endl;
		bool can_load = true;
		if (!is_same_sat(filename_sat, hclauses)) {
			log << "Hash collision! Solution " << filename << " was generated for different input!"
				<< std::endl;
			can_load = false;
			for (int i = 0; i < 100; ++i) {
				filename = get_filename(hash, ".sol.col" + std::to_string(i));
				filename_sat = get_filename(hash, ".sat.col" + std::to_string(i));
				if (!std::filesystem::exists(filename)) {
					break;
				} else if (is_same_sat(filename_sat, hclauses)) {
					log << "Will load solution from " << filename
						<< " that was generated for the same input." << std::endl;
					can_load = true;
					break;
				}
			}
		}
		if (can_load) {
			auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
			std::ifstream f(filename);
			int id;
			while (f >> id) {
				I.DS.push_back(id);
				l << " " << id;
			}
			l << "\n";
			if (I.DS.size() > before) {
				log << "Updated DS (cached): " << before << "+" << (I.DS.size() - before) << "="
					<< I.DS.size() << std::endl;
				return true;
			} else {
				log << "Cached file seems empty, discarding!" << std::endl;
			}
		}
	}

	log << "Will cache solution in " << filename << std::endl;
	std::filesystem::create_directory("cache");
	dump_sat(filename_sat, hclauses);
	return false;
}
#endif


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
#ifdef EMS_CACHE
	std::vector<std::vector<int>> hclauses;
	hclauses.reserve(I.G.numberOfNodes());
#endif
	std::vector<int> clause;
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
		clause.clear();
		clause.reserve(v->indeg() + 1);
#ifdef EMS_CACHE
		hclauses.emplace_back();
		hclauses.back().reserve(v->indeg() + 1);
#endif
		if (!I.is_subsumed[v]) {
			clause.push_back(varmap[v]);
#ifdef EMS_CACHE
			hclauses.back().push_back(I.node2ID[v]);
#endif
		}
		forAllInAdj(v, [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!I.is_subsumed[w]) {
				clause.push_back(varmap[w]);
#ifdef EMS_CACHE
				hclauses.back().push_back(I.node2ID[w]);
#endif
			}
			return true;
		});
		solver.addClause(clause); //hard clause
#ifdef EMS_CACHE
		std::sort(hclauses.back().begin(), hclauses.back().end());
#endif
	}
#ifdef EMS_CACHE
	std::string filename;
	if (try_load_solution(I, hclauses, filename)) {
		return;
	}
	std::ofstream f(filename);
#endif

	solver.setTargetComputationTime(std::min(10*60, I.G.numberOfNodes() / 10 + 1));
	std::cout.setstate(std::ios::failbit); // https://stackoverflow.com/a/8246430
	bool solved = solver.solve();
	std::cout.clear();
	if (!solved) {
		throw std::runtime_error("EvalMaxSAT didn't find optimal result!");
	}

	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : I.G.nodes) {
		if (!I.is_subsumed[v]) {
			if (solver.getValue(varmap[v])) {
				I.DS.push_back(I.node2ID[v]);
				l << " " << I.node2ID[v];
#ifdef EMS_CACHE
				f << " " << I.node2ID[v];
#endif
			}
		}
	}
	l << "\n";
	log << "Updated DS (solved " << solved << "): " << before << "+" << (I.DS.size() - before)
		<< "=" << I.DS.size() << std::endl;
}

#ifdef USE_ORTOOLS
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
#	ifdef EMS_CACHE
	std::vector<std::vector<int>> hclauses;
	hclauses.reserve(I.G.numberOfNodes());
#	endif
	MPObjective* const objective = solver->MutableObjective();
	objective->MinimizeLinearExpr(obj);
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
		LinearExpr expr;
#	ifdef EMS_CACHE
		hclauses.emplace_back();
		hclauses.back().reserve(v->indeg() + 1);
#	endif
		if (!I.is_subsumed[v]) {
			expr += varmap[v];
#	ifdef EMS_CACHE
			hclauses.back().push_back(I.node2ID[v]);
#	endif
		}
		forAllInAdj(v, [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!I.is_subsumed[w]) {
				expr += varmap[w];
#	ifdef EMS_CACHE
				hclauses.back().push_back(I.node2ID[w]);
#	endif
			}
			return true;
		});
		solver->MakeRowConstraint(expr >= 1); //hard clause
#	ifdef EMS_CACHE
		std::sort(hclauses.back().begin(), hclauses.back().end());
#	endif
	}
#	ifdef EMS_CACHE
	std::string filename;
	if (try_load_solution(I, hclauses, filename)) {
		return;
	}
	std::ofstream f(filename);
#	endif

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
				I.DS.push_back(I.node2ID[v]);
				l << " " << I.node2ID[v];
#	ifdef EMS_CACHE
				f << " " << I.node2ID[v];
#	endif
			}
		}
	}
	l << "\n";
	log << "Updated DS (solved " << result_status << "): " << before << "+"
		<< (I.DS.size() - before) << "=" << I.DS.size() << std::endl;
	solver->Clear();
}
#endif
