#include "ogdf_treewidth.h"

static int terminate_counter = 0;
static std::mutex mtx;
std::unique_ptr<htd::LibraryInstance> manager(htd::createManagementInstance(htd::Id::FIRST));

void ReductionTreeDecomposition::computeDecomposition() {
	std::srand(0);

	// Create a new graph instance which can handle (multi-)hyperedges.
	graph = manager->graphFactory().createInstance();
	graph->addVertices(G.numberOfNodes());
	size_t cnt = 0;
	for (ogdf::node v : G.nodes) {
		nodeid[v] = ++cnt;
		idnode[cnt] = v;
	}

	for (auto u : G.nodes) {
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			auto v = adj->twinNode();
			// if (u < v || I.reverse_edge[adj->theEdge()] == nullptr) {
			graph->addEdge(nodeid[u], nodeid[v]);
			// }
			return true;
		});
	}


	FitnessFunction fitnessFunction;
	htd::TreeDecompositionOptimizationOperation* operation =
			new htd::TreeDecompositionOptimizationOperation(manager.get(), fitnessFunction.clone());
	operation->setManagementInstance(manager.get());
	operation->setVertexSelectionStrategy(new htd::RandomVertexSelectionStrategy(10));
	operation->addManipulationOperation(
			new htd::NormalizationOperation(manager.get(), true, false, true, true));
	manager->orderingAlgorithmFactory().setConstructionTemplate(
			new htd::MinFillOrderingAlgorithm(manager.get()));
	htd::ITreeDecompositionAlgorithm* baseAlgorithm =
			manager->treeDecompositionAlgorithmFactory().createInstance();
	baseAlgorithm->addManipulationOperation(operation);
	htd::IterativeImprovementTreeDecompositionAlgorithm algorithm(manager.get(), baseAlgorithm,
			fitnessFunction.clone());
	algorithm.setIterationCount(10);
	algorithm.setNonImprovementLimit(3);

	std::packaged_task<int(int)> task([&](int cnt) {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		mtx.lock();
		int x = terminate_counter;
		if (x == cnt) {
			manager->terminate();
		}
		mtx.unlock();
		return 0;
	});
	mtx.lock();
	std::future<int> f1 = task.get_future(); // get a future
	std::thread t(std::move(task), terminate_counter); // launch on a thread
	t.detach();
	mtx.unlock();

	std::size_t optimalBagSize = (std::size_t)-1;
	decomposition = algorithm.computeDecomposition(*graph,
			[&](const htd::IMultiHypergraph& graph, const htd::ITreeDecomposition& decomposition,
					const htd::FitnessEvaluation& fitness) {
				// Disable warnings concerning unused variables.
				HTD_UNUSED(graph)
				HTD_UNUSED(decomposition)

				std::size_t bagSize = -fitness.at(0);

				/**
				 *  After each improvement we print the current optimal
				 *  width + 1 and the time when the decomposition was found.
				 */
				if (bagSize < optimalBagSize) {
					// optimalBagSize = bagSize;

					// std::chrono::milliseconds::rep msSinceEpoch =
					// 		std::chrono::duration_cast<std::chrono::milliseconds>(
					// 				std::chrono::system_clock::now().time_since_epoch())
					// 				.count();

					// std::cout << "c status " << optimalBagSize << " " << msSinceEpoch << std::endl;
				}
			});
	mtx.lock();
	terminate_counter++;
	mtx.unlock();

	if (decomposition != nullptr) {
		if (!manager->isTerminated() || algorithm.isSafelyInterruptible()) {
			// Print the height of the decomposition to stdout.
			// std::cout << decomposition->height() << std::endl;

			// Print the size of the largest bag of the decomposition to stdout.
			// std::cout << decomposition->maximumBagSize() << std::endl;
			treewidth = decomposition->maximumBagSize();
		} else {
			delete decomposition;
			decomposition = nullptr;
		}
	} else {
		log << "Interrupted, no decomposition found!" << std::endl;
	}
}

void ReductionTreeDecomposition::solveDPExact() {
	OGDF_ASSERT(decomposition != nullptr);
	htd::PostOrderTreeTraversal traversal;
	std::vector<std::unordered_map<size_t, size_t>> DP;

	traversal.traverse(*decomposition, [&](htd::vertex_t v, htd::vertex_t parent, std::size_t depth) {
		OGDF_ASSERT(decomposition->isVertex(v));

		auto bagvertices = decomposition->bagContent(v);
		if (DP.size() <= std::max(v, parent)) {
			DP.resize(std::max(v, parent) + 1);
		}

		// Initialize DP for leaf
		if (decomposition->isLeaf(v)) {
			ogdf::node voriginal = idnode[bagvertices[0]];
			OGDF_ASSERT(I.is_subsumed.graphOf() == I.is_dominated.graphOf());
			OGDF_ASSERT(I.is_subsumed.graphOf() == voriginal->graphOf());
			log << I.G.numberOfNodes() << std::endl;
			if (!I.is_subsumed[voriginal]) {
				DP[v][2] = 1;
			}
			if (I.is_dominated[voriginal]) {
				DP[v][1] = 0;
			}

			DP[v][0] = 0;
		} else if (decomposition->isForgetNode(v)) {
			auto child = decomposition->childAtPosition(v, 0);
			auto childvertices = decomposition->bagContent(child);
			OGDF_ASSERT(childvertices.size() == bagvertices.size() + 1);
			OGDF_ASSERT(decomposition->forgottenVertices(v).size() == 1);
			auto forgottenvertex = *decomposition->forgottenVertices(v).begin();
		}
	});
}