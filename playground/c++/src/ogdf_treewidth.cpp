#include "ogdf_treewidth.h"

static int terminate_counter = 0;
static std::mutex mtx;
std::unique_ptr<htd::LibraryInstance> manager(htd::createManagementInstance(htd::Id::FIRST));

void ReductionTreeDecomposition::computeDecomposition() {
	std::srand(0);

	// Create a new graph instance which can handle (multi-)hyperedges.
	graph = manager->graphFactory().createInstance(); // Use Multigraph! Graph checks if parallel edges, which is super slow!
	graph->addVertices(G.numberOfNodes());
	size_t cnt = 0;
	for (ogdf::node v : G.nodes) {
		nodeid[v] = ++cnt;
		idnode[cnt] = v;
	}
	ogdf::NodeSet<true> added(G);
	size_t i = 0;
	for (auto u : G.nodes) {
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			auto v = adj->twinNode();
			if (u < v) {
				graph->addEdgeWithoutCheck(nodeid[u], nodeid[v]);
				added.insert(v);
			}
			return true;
		});
		forAllInAdj(u, [&](ogdf::adjEntry adj) {
			auto v = adj->twinNode();
			if (u < v && !added.isMember(v)) {
				graph->addEdgeWithoutCheck(nodeid[u], nodeid[v]);
			}
			return true;
		});
		added.clear();
	}


	FitnessFunction fitnessFunction;
	htd::TreeDecompositionOptimizationOperation* operation =
			new htd::TreeDecompositionOptimizationOperation(manager.get(), fitnessFunction.clone());
	operation->setManagementInstance(manager.get());
	operation->setVertexSelectionStrategy(new htd::RandomVertexSelectionStrategy(10));
	operation->addManipulationOperation(
			new htd::NormalizationOperation(manager.get(), false, false, true, true));
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
		log << "Waiting for termination signal..." << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(1));
		mtx.lock();
		int x = terminate_counter;
		if (x == cnt) {
			log << "Terminating tree decomposition computation!" << std::endl;
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

TW_SIGNATURE_TYPE fromSigWithInsert(const std::vector<TW_SIGNATURE_TYPE>& sig, const size_t index,
		TW_SIGNATURE_TYPE status) {
	TW_SIGNATURE_TYPE id = 0;
	if (index == sig.size()) {
		id = status;
	}
	for (int i = sig.size() - 1; i >= 0; i--) {
		id *= 3;
		id += sig[i];
		if (i == index) {
			id *= 3;
			id += status;
		}
	}
	return id;
}

TW_SIGNATURE_TYPE fromSigVec(const std::vector<TW_SIGNATURE_TYPE>& sig) {
	TW_SIGNATURE_TYPE id = 0;
	for (int i = sig.size() - 1; i >= 0; i--) {
		id *= 3;
		id += sig[i];
	}
	return id;
}

std::vector<TW_SIGNATURE_TYPE> toSigVec(TW_SIGNATURE_TYPE id, size_t size) {
	std::vector<TW_SIGNATURE_TYPE> sig(size);
	for (size_t i = 0; i < size; i++) {
		sig[i] = id % 3;
		id /= 3;
	}
	return sig;
}

TW_SIGNATURE_TYPE sigAt(TW_SIGNATURE_TYPE id, size_t i) {
	return (id / (size_t)pow(3, i)) % 3; // TODO: maybe precompute powers
}

TW_SIGNATURE_TYPE stripSigAt(TW_SIGNATURE_TYPE id, size_t i) {
	return id % (size_t)pow(3, i) + // bits before index
			+(id / (size_t)pow(3, i + 1)) * (size_t)pow(3, i); // bits after index
}

size_t cntUndominated(TW_SIGNATURE_TYPE id, size_t vecsize) {
	size_t cnt = 0;
	for (size_t i = 0; i < vecsize; i++) {
		if (id % 3 == TW_WAITING) {
			cnt++;
		}
		id /= 3;
	}
	return cnt;
}

size_t cntInDs(TW_SIGNATURE_TYPE id) {
	size_t cnt = 0;
	while (id > 0) {
		if (id % 3 == TW_INDS) {
			cnt++;
		}
		id /= 3;
	}
	return cnt;
}

inline void updateTable(std::unordered_map<TW_SIGNATURE_TYPE,
								std::tuple<size_t, TW_SIGNATURE_TYPE, TW_SIGNATURE_TYPE>>& DP,
		TW_SIGNATURE_TYPE sig, size_t val, TW_SIGNATURE_TYPE sig2, TW_SIGNATURE_TYPE sig3) {
	if (DP.find(sig) == DP.end()) {
		DP[sig] = {val, sig2, sig3};
	} else if (std::get<0>(DP[sig]) > val) {
		DP[sig] = {val, sig2, sig3};
	}
}

int ReductionTreeDecomposition::solveDPExact() {
	OGDF_ASSERT(decomposition != nullptr);
	htd::PostOrderTreeTraversal traversal;
	std::vector<std::unordered_map<TW_SIGNATURE_TYPE, std::tuple<size_t, TW_SIGNATURE_TYPE, TW_SIGNATURE_TYPE>>>
			DP; // TODO: instead of unordered_map we could use vector, this would cost more memory
	std::vector<std::vector<ogdf::node>> bag_nodes;
	ogdf::NodeSet<true> currBagNodes(G);
	ogdf::NodeArray<size_t> currSigIndex(G);
	std::set<htd::vertex_t> visited;
	traversal.traverse(*decomposition,
			[&](htd::vertex_t curbag, htd::vertex_t bagparent, std::size_t depth) {
				if (visited.find(curbag) != visited.end()) {
					return;
				}
				visited.insert(curbag);
				OGDF_ASSERT(decomposition->isVertex(curbag));
				if (DP.size() <= curbag) {
					DP.resize(curbag + 1);
					bag_nodes.resize(curbag + 1);
				}
				for (auto u : decomposition->bagContent(curbag)) {
					ogdf::node node = idnode[u];
					bag_nodes[curbag].push_back(node);
				}
				std::sort(bag_nodes[curbag].begin(), bag_nodes[curbag].end());

				if (decomposition->isLeaf(curbag)) {
					// BASE CASE: As this is a nice tree decomposition, bag should only contain one node
					OGDF_ASSERT(bag_nodes[curbag].size() == 1);
					ogdf::node voriginal = bag_nodes[curbag][0];
					OGDF_ASSERT(I.is_subsumed.graphOf() == I.is_dominated.graphOf());
					OGDF_ASSERT(I.is_subsumed.graphOf() == voriginal->graphOf());
					if (!I.is_dominated[voriginal]) {
						DP[curbag][TW_WAITING] = {0, 0, 0};
					} else {
						DP[curbag][TW_DOMINATED] = {0, 0, 0};
					}
					if (!I.is_subsumed[voriginal]) {
						DP[curbag][TW_INDS] = {1, 0, 0};
					}
				} else if (decomposition->isForgetNode(curbag)) {
					// FORGET NODE
					OGDF_ASSERT(decomposition->children(curbag).size() == 1);
					auto bagchild = decomposition->childAtPosition(curbag, 0);
					OGDF_ASSERT(bag_nodes[curbag].size() + 1 == bag_nodes[bagchild].size());
					OGDF_ASSERT(decomposition->forgottenVertices(curbag).size() == 1);
					auto forgotten = *decomposition->forgottenVertices(curbag).begin();
					ogdf::node forgottenvertex = idnode[forgotten];
					size_t forgottenindex = lower_bound(bag_nodes[bagchild].begin(),
													bag_nodes[bagchild].end(), forgottenvertex)
							- bag_nodes[bagchild].begin();
					OGDF_ASSERT(forgottenindex < bag_nodes[bagchild].size());
					OGDF_ASSERT(bag_nodes[bagchild][forgottenindex] == forgottenvertex);

					for (auto& [sig, val] : DP[bagchild]) {
						if (sigAt(sig, forgottenindex) != TW_WAITING) {
							// forgetting vertex is fine because it is already dominated
							TW_SIGNATURE_TYPE newsig = stripSigAt(sig, forgottenindex);
							updateTable(DP[curbag], newsig, std::get<0>(val), sig, 0);
						}
					}
				} else if (decomposition->isIntroduceNode(curbag)) {
					// INTRODUCE NODE
					OGDF_ASSERT(decomposition->children(curbag).size() == 1);
					auto bagchild = decomposition->childAtPosition(curbag, 0);
					OGDF_ASSERT(bag_nodes[curbag].size() == bag_nodes[bagchild].size() + 1);
					OGDF_ASSERT(decomposition->introducedVertices(curbag).size() == 1);
					auto introduced = *decomposition->introducedVertices(curbag).begin();
					ogdf::node introducedvertex = idnode[introduced];
					for (auto u : bag_nodes[bagchild]) {
						currBagNodes.insert(u);
					}
					size_t introducedindex = lower_bound(bag_nodes[curbag].begin(),
													 bag_nodes[curbag].end(), introducedvertex)
							- bag_nodes[curbag].begin();
					OGDF_ASSERT(bag_nodes[curbag][introducedindex] == introducedvertex);
					for (size_t i = 0; i < bag_nodes[bagchild].size(); i++) {
						currSigIndex[bag_nodes[bagchild][i]] = i;
					}
					for (auto& [sig, val] : DP[bagchild]) {
						auto sigvec = toSigVec(sig, bag_nodes[bagchild].size());

						if (I.is_dominated[introducedvertex]) {
							auto newsig = fromSigWithInsert(sigvec, introducedindex, TW_DOMINATED);
							updateTable(DP[curbag], newsig, std::get<0>(val), sig, 0);
						} else {
							bool anyininds = false;
							forAllInAdj(introducedvertex, [&](ogdf::adjEntry adj) {
								auto v = adj->twinNode();
								if (currBagNodes.isMember(v) && sigvec[currSigIndex[v]] == TW_INDS) {
									anyininds = true;
									return false;
								}
								return true;
							});
							if (anyininds) {
								auto newsig =
										fromSigWithInsert(sigvec, introducedindex, TW_DOMINATED);
								updateTable(DP[curbag], newsig, std::get<0>(val), sig, 0);
							} else {
								auto newsig = fromSigWithInsert(sigvec, introducedindex, TW_WAITING);
								updateTable(DP[curbag], newsig, std::get<0>(val), sig, 0);
							}
						}

						if (!I.is_subsumed[introducedvertex]) {
							forAllOutAdj(introducedvertex, [&](ogdf::adjEntry adj) {
								auto v = adj->twinNode();
								if (currBagNodes.isMember(v)
										&& sigvec[currSigIndex[v]] == TW_WAITING) {
									sigvec[currSigIndex[v]] = TW_DOMINATED;
								}
								return true;
							});
							auto newsig = fromSigWithInsert(sigvec, introducedindex, TW_INDS);
							updateTable(DP[curbag], newsig, std::get<0>(val) + 1, sig, 0);
						}
					}

					for (auto u : bag_nodes[bagchild]) {
						currBagNodes.remove(u);
					}
				} else if (decomposition->isJoinNode(curbag)) {
					// JOIN NODE
					OGDF_ASSERT(decomposition->children(curbag).size() == 2);
					auto bagchild1 = decomposition->childAtPosition(curbag, 0);
					auto bagchild2 = decomposition->childAtPosition(curbag, 1);
					for (auto& [sig, val] : DP[bagchild1]) {
						auto sigvec = toSigVec(sig, bag_nodes[bagchild1].size());
						for (auto& [sig2, val2] : DP[bagchild2]) {
							auto sigvec2 = toSigVec(sig2, bag_nodes[bagchild2].size());
							size_t cnt_both_in = 0;
							for (size_t i = 0; i < sigvec.size(); i++) {
								if (sigvec[i] == TW_INDS && sigvec2[i] == TW_INDS) {
									cnt_both_in++;
								}
								sigvec2[i] = std::max(sigvec[i], sigvec2[i]);
							}
							updateTable(DP[curbag], fromSigVec(sigvec2),
									std::get<0>(val) + std::get<0>(val2) - cnt_both_in, sig, sig2);
						}
					}
				} else {
					// There seem to be nodes with children having the exact same bag
					OGDF_ASSERT(decomposition->children(curbag).size() == 1);
					auto bagchild = decomposition->childAtPosition(curbag, 0);
					OGDF_ASSERT(bag_nodes[curbag].size() == bag_nodes[bagchild].size());
#ifdef OGDF_DEBUG
					for (int i = 0; i < bag_nodes[curbag].size(); i++) {
						OGDF_ASSERT(bag_nodes[curbag][i] == bag_nodes[bagchild][i]);
					}
#endif
					for (auto& [sig, val] : DP[bagchild]) {
						auto sigvec = toSigVec(sig, bag_nodes[bagchild].size());
						DP[curbag][sig] = {std::get<0>(val), sig, 0};
					}
				}
			});
	int ans = std::numeric_limits<int>::max();
	TW_SIGNATURE_TYPE anssig;
	for (auto& [sig, val] : DP[decomposition->root()]) {
		if (cntUndominated(sig, decomposition->bagContent(decomposition->root()).size()) == 0) {
			if ((int)std::get<0>(val) < ans) {
				anssig = sig;
				ans = (int)std::get<0>(val);
			}
		}
	}
	std::vector<TW_SIGNATURE_TYPE> chosensig;
	chosensig.resize(decomposition->root() + 1);
	chosensig[decomposition->root()] = anssig;
	htd::PreOrderTreeTraversal traversalpre;
	int sizebefore = I.DS.size();
	std::set<ogdf::node> addednodes;
	traversalpre.traverse(*decomposition,
			[&](htd::vertex_t curbag, htd::vertex_t bagparent, std::size_t depth) {
				OGDF_ASSERT(decomposition->isVertex(curbag));
				if (decomposition->isLeaf(curbag) || decomposition->isIntroduceNode(curbag)) {
					htd::vertex_t introduced;
					if (decomposition->isLeaf(curbag)) {
						introduced = decomposition->bagContent(curbag)[0];
					} else {
						introduced = *decomposition->introducedVertices(curbag).begin();
					}
					ogdf::node introducedvertex = idnode[introduced];
					size_t introducedindex = lower_bound(bag_nodes[curbag].begin(),
													 bag_nodes[curbag].end(), introducedvertex)
							- bag_nodes[curbag].begin();
					auto status = sigAt(chosensig[curbag], introducedindex);
					if (status == TW_INDS && addednodes.find(introducedvertex) == addednodes.end()) {
						addednodes.insert(introducedvertex);
						I.DS.push_back(I.node2ID(introducedvertex));
					}
				}
				size_t i = 0;
				for (auto u : decomposition->children(curbag)) {
					if (chosensig.size() <= u) {
						chosensig.resize(u + 1);
					}
					auto t = DP[curbag][chosensig[curbag]];
					if (i == 0) {
						chosensig[u] = std::get<1>(t);
					} else {
						chosensig[u] = std::get<2>(t);
					}
					i++;
				}
			});

	int cntadded = I.DS.size() - sizebefore;
	OGDF_ASSERT(cntadded == ans);
	return ans;
}