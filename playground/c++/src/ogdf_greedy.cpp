#include "ogdf_greedy.hpp"

void solveGreedy(Instance& I) {
	log << "Solving Greedy with " << I.G.numberOfNodes() << " nodes" << std::endl;
	log << "Old DS size: " << I.DS.size() << std::endl;
	while (!I.G.empty()) {
		ogdf::node bestVertex = nullptr;
		int bestValue = 0;
		for (auto v : I.G.nodes) {
			int wouldbecovered = v->outdeg();
			if (!I.is_dominated[v]) {
				wouldbecovered++;
			}
			if (wouldbecovered > bestValue) {
				bestValue = wouldbecovered;
				bestVertex = v;
			}
		}
		if (bestVertex == nullptr) {
			break;
		}
		I.addToDominatingSet(bestVertex);
	}
	log << "New DS size: " << I.DS.size() << std::endl;
}
