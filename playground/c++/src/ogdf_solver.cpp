#include "ogdf_solver.hpp"

#include "ogdf_greedy.hpp"
#include "ogdf_maxsat.hpp"

void reduceAndSolve(Instance& I, int d) {
	bool changed = true;
	int m, n, i = 0;
	logger.localLogLevel(ogdf::Logger::Level::Default);
	while (changed) {
		n = I.G.numberOfNodes();
		m = I.G.numberOfEdges();
		changed = false;
		log << "Reduce iteration " << i << " depth " << d << ": " << n << " nodes, " << m
			<< " edges" << std::endl;

		// this reduction is so cheap, make sure we really have no isolated vertices before decomposing components
		while (I.reductionExtremeDegrees()) {
			changed = true;
		}

		std::list<Instance> comps = I.decomposeConnectedComponents();
		if (!comps.empty()) {
			log << comps.size() << " connected components" << std::endl;
			I.clear(); // save some memory

			int c = 0;
			for (auto& comp : comps) {
				log << "Connected component " << c << std::endl;
				ogdf::Logger::Indent _(logger);
				++c;

				// and now recurse
				reduceAndSolve(comp, d + 1);
				I.addToDominatingSet(comp.DS.begin(), comp.DS.end(),
						"connected component " + std::to_string(c));
			}
			return;
		}

		if (I.reductionBCTree(d)) {
			changed = true;
		} else {
			if (I.reductionSubsumption()) {
				changed = true; // TODO: this rule is still slow for large graphs
			}
			if (I.reductionStrongSubsumption()) {
				changed = true; // TODO: this rule is still slow for large graphs
			}
		}

		OGDF_ASSERT(!changed || I.G.numberOfNodes() < n || I.G.numberOfEdges() < m);
		++i;
	}

	if (I.G.numberOfNodes() < 1) {
		log << "Reduced instance is empty with DS " << I.DS.size() << "!" << std::endl;
		return;
	}
#ifdef OGDF_DEBUG
	auto [can, need] = I.dominationStats();
	log << "Reduced instance contains " << n << " nodes, " << m << " edges. " << need
		<< " vertices need to be dominated, " << can << " are eligible for the DS." << std::endl;
	// I.dumpBCTree();
#endif

	// if (I.G.numberOfNodes() > SMALL_BLOCK) {
	//     log << "Using greedy approximation for large block!" << std::endl;
	//     solveGreedy(I);
	//     return;
	// }

	// now to solving...
#ifdef USE_ORTOOLS
	solvecpsat(I);
#else
	solveEvalMaxSat(I);
#endif
}
