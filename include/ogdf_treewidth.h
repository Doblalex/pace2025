#pragma once

#include <chrono>
#include <csignal>
#include <future>
#include <iostream>
#include <memory>
#include <queue>
#include <stack>

#include <ogdf/basic/GraphSets.h>
#include <htd/main.hpp>

#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"

typedef size_t TW_SIGNATURE_TYPE;
#define TW_WAITING 0
#define TW_DOMINATED 1
#define TW_INDS 2
#define TW_QM 1

struct hash_pair final {
	template<class TFirst, class TSecond>
	size_t operator()(const std::pair<TFirst, TSecond>& p) const noexcept {
		uintmax_t hash = std::hash<TFirst> {}(p.first);
		hash <<= sizeof(uintmax_t) * 4;
		hash ^= std::hash<TSecond> {}(p.second);
		return std::hash<uintmax_t> {}(hash);
	}
};

class ReductionTreeDecomposition {
public:
	ogdf::Graph& G;
	Instance& I;
	htd::ITreeDecomposition* decomposition = nullptr;
	htd::IMutableGraph* graph;
	std::vector<std::unordered_map<TW_SIGNATURE_TYPE, std::tuple<size_t, TW_SIGNATURE_TYPE, TW_SIGNATURE_TYPE>>>
			DP; // TODO: instead of unordered_map we could use vector, this would cost more memory
	// std::vector<std::unordered_map<std::pair<TW_SIGNATURE_TYPE, size_t>, u_int64_t>> DPCNT;
	std::vector<size_t> minkappa;
	std::vector<std::vector<ogdf::node>> bag_nodes;
	ogdf::NodeArray<size_t> nodeid;
	std::vector<ogdf::node> idnode;
	ogdf::NodeSet currBagNodes;
	ogdf::NodeArray<size_t> currSigIndex;
	std::vector<TW_SIGNATURE_TYPE> chosensig;
	int treewidth = -1;

	ReductionTreeDecomposition(ogdf::Graph& G, Instance& instance)
		: G(G)
		, I(instance)
		, nodeid(G, 0)
		, idnode(G.numberOfNodes() + 1)
		, currBagNodes(G)
		, currSigIndex(G) { }

	~ReductionTreeDecomposition() {
		if (decomposition) {
			delete decomposition;
		}
		if (graph) {
			delete graph;
		}
	}

	void computeDecomposition();

	int solveDPExact();

private:
	// size_t getCntDP(size_t bag, TW_SIGNATURE_TYPE sig, size_t cnt) {
	// 	auto it = DPCNT[bag].find({sig, cnt});
	// 	if (it != DPCNT[bag].end()) {
	// 		return it->second;
	// 	}
	// 	return 0;
	// }

	// void addCntDP(size_t bag, TW_SIGNATURE_TYPE sig, size_t cnt, size_t val) {
	// 	auto it = DPCNT[bag].find({sig, cnt});
	// 	if (it != DPCNT[bag].end()) {
	// 		it->second += val;
	// 	} else {
	// 		DPCNT[bag][{sig, cnt}] = val;
	// 	}
	// }

	void handleLeaf(htd::vertex_t curbag);

	void handleForgetNode(htd::vertex_t curbag);

	void handleIntroduceNode(htd::vertex_t curbag);

	void handleJoinNode(htd::vertex_t curbag);

	void handleCopyNode(htd::vertex_t curbag);

	void backTrackJoinNode(htd::vertex_t curbag);
};

class FitnessFunction : public htd::ITreeDecompositionFitnessFunction {
public:
	FitnessFunction(void) { }

	~FitnessFunction() { }

	htd::FitnessEvaluation* fitness(const htd::IMultiHypergraph& graph,
			const htd::ITreeDecomposition& decomposition) const {
		HTD_UNUSED(graph)

		/**
              * Here we specify the fitness evaluation for a given decomposition.
              * In this case, we select the maximum bag size and the height.
              */
		return new htd::FitnessEvaluation(2, -(double)(decomposition.maximumBagSize()),
				-(double)(decomposition.height()));
	}

	FitnessFunction* clone(void) const { return new FitnessFunction(); }
};
