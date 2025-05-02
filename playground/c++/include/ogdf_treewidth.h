#pragma once

#include <chrono>
#include <csignal>
#include <future>
#include <iostream>
#include <memory>
#include <queue>
#include <stack>
#include <thread>

#include "ogdf/basic/GraphSets.h"
#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"
#include <htd/main.hpp>

typedef size_t TW_SIGNATURE_TYPE;
#define TW_WAITING 0
#define TW_DOMINATED 1
#define TW_INDS 2

class ReductionTreeDecomposition {
public:
	ogdf::Graph& G;
	Instance& I;
	htd::ITreeDecomposition* decomposition = nullptr;
	htd::IMutableMultiGraph* graph;
	ogdf::NodeArray<size_t> nodeid;
	std::vector<ogdf::node> idnode;
	int treewidth = -1;

	ReductionTreeDecomposition(ogdf::Graph& G, Instance& instance)
		: G(G), I(instance), nodeid(G, 0), idnode(G.numberOfNodes() + 1) { }

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