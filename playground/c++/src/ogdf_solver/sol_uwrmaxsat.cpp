#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"
#include <uwrmaxsat/ipamir.h>

void solveIPAMIR(Instance& I) {
	void* ipamir = ipamir_init();
	int before = I.DS.size();
	log << "Solving IPAMIR " << std::string(ipamir_signature()) << " with " << I.G.numberOfNodes()
		<< " nodes" << std::endl;
	if (!ipamir) {
		std::cerr << "solver not available" << std::endl;
		exit(1);
	}
	for (auto v : I.G.nodes) {
		if (I.is_subsumed[v]) {
			continue;
		}
		ipamir_add_soft_lit(ipamir, I.node2ID[v], 1);
	}
#ifdef SAT_CACHE
	std::vector<std::vector<int>> hclauses;
	hclauses.reserve(I.G.numberOfNodes());
#endif
	for (auto v : I.G.nodes) {
		if (I.is_dominated[v]) {
			continue;
		}
#ifdef SAT_CACHE
		hclauses.emplace_back();
		hclauses.back().reserve(v->indeg() + 1);
#endif
		if (!I.is_subsumed[v]) {
			ipamir_add_hard(ipamir, I.node2ID[v]);
#ifdef SAT_CACHE
			hclauses.back().push_back(I.node2ID[v]);
#endif
		}
		auto add_neigh = [&](ogdf::adjEntry adj) {
			auto w = adj->twinNode();
			if (!adj->isSource() && !I.is_subsumed[w]) {
				ipamir_add_hard(ipamir, I.node2ID[w]);
#ifdef SAT_CACHE
				hclauses.back().push_back(I.node2ID[w]);
#endif
			}
			return true;
		};
		forAllInAdj(v, add_neigh);
		ogdf::safeForEach(I.hidden_edges.adjEntries(v), add_neigh);
		ipamir_add_hard(ipamir, 0);
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

	int result = ipamir_solve(ipamir);
	if (result != 30) {
		ipamir_release(ipamir);
		std::cerr << "result_status " << result << std::endl;
		throw std::runtime_error(
				"IPAMIR solver " + std::string(ipamir_signature()) + " didn't find optimal result!");
	}

	auto& l = logger.lout(ogdf::Logger::Level::Minor) << "Add to DS:";
	for (auto v : I.G.nodes) {
		if (!I.is_subsumed[v]) {
			if (ipamir_val_lit(ipamir, I.node2ID[v]) > 0) {
				I.DS.insert(I.node2ID[v]);
				l << " " << I.node2ID[v];
#ifdef SAT_CACHE
				f << " " << I.node2ID[v];
#endif
			}
		}
	}
	l << "\n";
	log << "Updated DS (solved " << result << "): " << before << "+" << (I.DS.size() - before)
		<< "=" << I.DS.size() << " (weight " << ipamir_val_obj(ipamir) << ")" << std::endl;
	OGDF_ASSERT(ipamir_val_obj(ipamir) == I.DS.size() - before);
	ipamir_release(ipamir);
}
