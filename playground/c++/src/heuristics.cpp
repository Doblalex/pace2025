#include "heuristics.hpp"

void Greedy(Instance* instance, VertexList& dominatingSet) {
	while (instance->n > 0) {
		VD bestVertex = nullptr;
		int bestValue = -1;
		BGL_FORALL_VERTICES(v, *instance->G, Graph) {
			int wouldbecovered = (*instance->G)[v].cnt_undominated_neighbors;
			if (!(*instance->G)[v].is_dominated) {
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
		unordered_set<VD> toDelete;
		dominatingSet.push_back((*instance->G)[bestVertex].id);
		instance->toDominatingSet(bestVertex, toDelete);
		instance->deleteVertices(toDelete);
	}
}