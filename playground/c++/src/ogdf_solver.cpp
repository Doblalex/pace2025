#include "ogdf_solver.hpp"

#include "ogdf_greedy.hpp"
#include "ogdf_maxsat.hpp"
#include "ogdf_treewidth.h"

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

		if (I.reductionNeighborhoodSubsets()) {
			changed = true;
		} else if (I.reductionContraction()) {
			changed = true;
		} else if (I.reductionSpecial1()) {
			changed = true;
		} else if (I.reductionBCTree(d)) {
			changed = true;
		}

		// some RRs only add subsumed / dominated nodes without actually deleting sth (or rather only their imaginary selfloop)
		//OGDF_ASSERT(!changed || I.G.numberOfNodes() < n || I.G.numberOfEdges() < m);
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
#endif

	// now to solving...

	ReductionTreeDecomposition rtd(I.G, I);
	rtd.computeDecomposition();
	int ansdp = -1;
	if (rtd.decomposition != nullptr) {
		log << "Decomposition found with treewidth " << rtd.treewidth << std::endl;
		if (rtd.treewidth <= 13) {
			log << "Solving with DP" << std::endl;
			ansdp = rtd.solveDPExact();
			return;
		}
	}

#ifdef USE_ORTOOLS
	solvecpsat(I);
// #elif USE_GUROBI
// 	solveGurobiExactGurobi(I);
#else
	int anssat = I.DS.size();
	solveEvalMaxSat(I);
	anssat = I.DS.size() - anssat;
	if (ansdp != -1) {
		log << anssat << " vs. " << ansdp << std::endl;
		OGDF_ASSERT(anssat == ansdp);
	}
#endif
}
