#include "ogdf_maxsat.hpp"

#include "EvalMaxSAT.h"

#define EMS_CACHE

#ifdef EMS_CACHE
static const uint64_t FNV1a_64_SEED = 0xcbf29ce484222325UL;

// https://stackoverflow.com/a/77342581
inline void FNV1a_64_update(uint64_t& h, uint64_t v) {
	h ^= v;
	h *= 0x00000100000001B3UL;
}
#endif

void solveEvalMaxSat(Instance& I) {
	EvalMaxSAT solver;
	int before = I.DS.size();
	log << "Solving MaxSat with " << I.G.numberOfNodes() << " nodes" << std::endl;
	ogdf::NodeArray<int> varmap(I.G, -1);
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		auto var = solver.newVar();
		varmap[v] = var;
		solver.addClause({-var}, 1); // soft clause
	}
#ifdef EMS_CACHE
	std::vector<std::vector<int>> clauses;
#endif
	std::vector<int> clause;
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
		clause.clear();
		clause.reserve(v->indeg() + 1);
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
#ifndef EMS_CACHE
	}
#else
		std::sort(clause.begin(), clause.end());
		clauses.push_back({});
		std::swap(clause, clauses.back());
	}
	std::sort(clauses.begin(), clauses.end());
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
	std::string filename;
	{
		std::stringstream fnstr;
		fnstr << "cache/";
		fnstr << std::setfill('0') << std::setw(sizeof(uint64_t) * 2) << std::hex << hash;
		fnstr << ".txt";
		filename = fnstr.str();
	}
	if (std::filesystem::exists(filename)) {
		log << "Found cached solution " << filename << std::endl;
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
			return;
		} else {
			log << "Cached file seems empty, discarding!" << std::endl;
		}
	}
	log << "Will cache solution in " << filename << std::endl;
	std::filesystem::create_directory("cache");
	std::ofstream f(filename);
#endif

	solver.setTargetComputationTime(10 * 60);
	std::cout.setstate(std::ios::failbit); // https://stackoverflow.com/a/8246430
	bool solved = solver.solve();
	std::cout.clear();

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
