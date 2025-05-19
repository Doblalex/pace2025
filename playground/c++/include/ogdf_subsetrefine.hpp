#pragma once
#include <set>

#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"

enum class RefineType { Subsume, Dominate };

class SubsetRefine {
	RefineType type;
	Instance& instance;

	ogdf::Graph& G;
	ogdf::NodeArray<std::vector<ogdf::node>> needRefineBy;
	ogdf::NodeArray<size_t> needRefineBySize;
	std::set<std::pair<size_t, ogdf::node>> needTouch;
	ogdf::Graph refineG;
	ogdf::NodeArray<ogdf::node> bagof; // node(G) -> node(refineG)
	ogdf::NodeArray<std::vector<ogdf::node>> bagNodeVec; // node(refineG) -> list of contained nodes from G
	ogdf::NodeArray<size_t> vecIndex; // node(G) -> index of node in its vector for swap and pop
	ogdf::NodeArray<ogdf::node> refinedBag; // node(refineG) -> pointer to node which contains the elements that got refined
	size_t cntedgesadded = 0;
	ogdf::NodeArray<size_t> cnttouched;
	size_t cntreduced = 0;

public:
	SubsetRefine(Instance& instance, RefineType type)
		: instance(instance)
		, G(instance.G)
		, bagof(G, nullptr)
		, bagNodeVec(refineG)
		, vecIndex(G, 0)
		, refinedBag(refineG, nullptr)
		, type(type)
		, cnttouched(G, 0) { }

	void init() {
		auto initbag = refineG.newNode();
		needRefineBy.init(G, std::vector<ogdf::node>());
		needRefineBySize.init(G, 0);
		for (auto u : G.nodes) {
			if ((!instance.is_dominated(u) && type == RefineType::Dominate)
					|| (!instance.is_subsumed(u) && type == RefineType::Subsume)) {
				bagof[u] = initbag;
				bagNodeVec[initbag].push_back(u);
				vecIndex[u] = bagNodeVec[initbag].size() - 1;

				if (type == RefineType::Subsume) {
					instance.forAllCanDominate(u, [&](ogdf::node adj) {
						needRefineBy[u].push_back(adj);
						return true;
					});
				} else {
					instance.forAllCanBeDominatedBy(u, [&](ogdf::node adj) {
						needRefineBy[u].push_back(adj);
						return true;
					});
				}
				needTouch.insert({needRefineBy[u].size(), u});
				needRefineBySize[u] = needRefineBy[u].size();
			}
		}
	}

	bool doReduce(ogdf::node u) {
		// bool reduce = type == RefineType::Subsume || instance.is_subsumed[u];
		bool reduce = true;
		if (reduce) {
			if (type == RefineType::Subsume) {
				instance.markSubsumed(u);
			} else {
				instance.markDominated(u, false);
			}
			cntreduced++;
		}
		auto& vec = bagNodeVec[bagof[u]];
		auto oldindex = vecIndex[u];
		vecIndex[vec[vec.size() - 1]] = oldindex;
		std::swap(vec[oldindex], vec[vec.size() - 1]);
		vec.pop_back();
		bagof[u] = nullptr;
		if (needTouch.find({needRefineBySize[u], u}) != needTouch.end()) {
			needTouch.erase({needRefineBySize[u], u});
		}
		return reduce;
	}

	size_t doRefinementReduction();

	void refineByNode(const ogdf::node& u);
};
