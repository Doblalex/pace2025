#include "ogdf_treewidth.h"

static int terminate_counter = 0;
static std::mutex mtx;
std::unique_ptr<htd::LibraryInstance> manager(htd::createManagementInstance(htd::Id::FIRST));

void ReductionTreeDecomposition::computeDecomposition() {
	std::srand(0);

	manager->reset();

	// Create a new graph instance which can handle (multi-)hyperedges.
	graph = manager->graphFactory().createInstance(); // Use Multigraph! Graph checks if parallel edges, which is super slow!
	graph->addVertices(G.numberOfNodes());
	size_t cnt = 0;
	for (ogdf::node v : G.nodes) {
		nodeid[v] = ++cnt;
		idnode[cnt] = v;
	}
	ogdf::NodeSet added(G);
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
	operation->setVertexSelectionStrategy(new htd::RandomVertexSelectionStrategy(3));
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
	
	std::thread([](int cnt) {
		log << "Waiting for termination signal..." << std::endl;
		std::this_thread::sleep_for(std::chrono::seconds(1));
		mtx.lock();
		if (terminate_counter == cnt) {
			log << "Terminating tree decomposition computation!" << std::endl;
			manager->terminate();
		}
		mtx.unlock();
	},
	terminate_counter).detach();
	// auto t = std::async(
	// 		std::launch::async,
	// 		[](int cnt) {
	// 			log << "Waiting for termination signal..." << std::endl;
	// 			std::this_thread::sleep_for(std::chrono::seconds(1));
	// 			mtx.lock();
	// 			if (terminate_counter == cnt) {
	// 				log << "Terminating tree decomposition computation!" << std::endl;
	// 				manager->terminate();
	// 			}
	// 			mtx.unlock();
	// 		},
	// 		terminate_counter);

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
	return (id / (size_t)pow((size_t)3, i)) % 3; // TODO: maybe precompute powers
}

TW_SIGNATURE_TYPE sigAtSet(TW_SIGNATURE_TYPE id, size_t i, TW_SIGNATURE_TYPE status) {
	return status * pow((size_t)3, i) + id % (size_t)pow((size_t)3, i)
			+ (id / (size_t)pow((size_t)3, i + 1)) * (size_t)pow((size_t)3, i + 1); // bits after index
}

TW_SIGNATURE_TYPE stripSigAt(TW_SIGNATURE_TYPE id, size_t i) {
	return id % (size_t)pow((size_t)3, i) + // bits before index
			+(id / (size_t)pow((size_t)3, i + 1)) * (size_t)pow((size_t)3, i); // bits after index
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

void ReductionTreeDecomposition::handleLeaf(htd::vertex_t curbag) {
	// BASE CASE: As this is a nice tree decomposition, bag should only contain one node
	OGDF_ASSERT(bag_nodes[curbag].size() == 1);
	ogdf::node voriginal = bag_nodes[curbag][0];
	OGDF_ASSERT(I.is_subsumed.graphOf() == I.is_dominated.graphOf());
	OGDF_ASSERT(I.is_subsumed.graphOf() == voriginal->graphOf());
	if (!I.is_dominated[voriginal]) {
		DP[curbag][TW_WAITING] = {0, 0, 0};
		// DPCNT[curbag][{TW_WAITING, 0}] = 1;
	} else {
		DP[curbag][TW_DOMINATED] = {0, 0, 0};
		// DPCNT[curbag][{TW_DOMINATED, 0}] = 1;
	}
	if (!I.is_subsumed[voriginal]) {
		// DPCNT[curbag][{TW_INDS, 1}] = 0;
		DP[curbag][TW_INDS] = {1, 0, 0};
	}
	minkappa[curbag] = 0;
}

void ReductionTreeDecomposition::handleForgetNode(htd::vertex_t curbag) {
	// FORGET NODE
	OGDF_ASSERT(decomposition->children(curbag).size() == 1);
	auto bagchild = decomposition->childAtPosition(curbag, 0);
	OGDF_ASSERT(bag_nodes[curbag].size() + 1 == bag_nodes[bagchild].size());
	OGDF_ASSERT(decomposition->forgottenVertices(curbag).size() == 1);
	auto forgotten = *decomposition->forgottenVertices(curbag).begin();
	ogdf::node forgottenvertex = idnode[forgotten];
	size_t forgottenindex =
			lower_bound(bag_nodes[bagchild].begin(), bag_nodes[bagchild].end(), forgottenvertex)
			- bag_nodes[bagchild].begin();
	OGDF_ASSERT(forgottenindex < bag_nodes[bagchild].size());
	OGDF_ASSERT(bag_nodes[bagchild][forgottenindex] == forgottenvertex);

	for (auto& [sig, val] : DP[bagchild]) {
		if (sigAt(sig, forgottenindex) != TW_WAITING) {
			// forgetting vertex is fine because it is already dominated
			TW_SIGNATURE_TYPE newsig = stripSigAt(sig, forgottenindex);
			updateTable(DP[curbag], newsig, std::get<0>(val), sig, 0);
			minkappa[curbag] = std::min(minkappa[curbag], std::get<0>(val));
		}
	}

	// for (auto& [sig_cnt, val] : DPCNT[bagchild]) {
	// 	if (sigAt(sig_cnt.first, forgottenindex) != TW_WAITING) {
	// 		// forgetting vertex is fine because it is already dominated
	// 		TW_SIGNATURE_TYPE newsig = stripSigAt(sig_cnt.first, forgottenindex);
	// 		addCntDP(curbag, newsig, sig_cnt.second, val);
	// 	}
	// }
}

void ReductionTreeDecomposition::handleIntroduceNode(htd::vertex_t curbag) {
	OGDF_ASSERT(decomposition->children(curbag).size() == 1);
	auto bagchild = decomposition->childAtPosition(curbag, 0);
	OGDF_ASSERT(bag_nodes[curbag].size() == bag_nodes[bagchild].size() + 1);
	OGDF_ASSERT(decomposition->introducedVertices(curbag).size() == 1);
	auto introduced = *decomposition->introducedVertices(curbag).begin();
	minkappa[curbag] = minkappa[bagchild];
	ogdf::node introducedvertex = idnode[introduced];
	for (auto u : bag_nodes[bagchild]) {
		currBagNodes.insert(u);
	}
	size_t introducedindex =
			lower_bound(bag_nodes[curbag].begin(), bag_nodes[curbag].end(), introducedvertex)
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
			auto newsig = fromSigWithInsert(sigvec, introducedindex, TW_WAITING);
			updateTable(DP[curbag], newsig, std::get<0>(val), sig,
					0); // need this for fast join nodes!
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
				auto newsig = fromSigWithInsert(sigvec, introducedindex, TW_DOMINATED);
				updateTable(DP[curbag], newsig, std::get<0>(val), sig, 0);
			}
		}

		if (!I.is_subsumed[introducedvertex]) {
			forAllOutAdj(introducedvertex, [&](ogdf::adjEntry adj) {
				auto v = adj->twinNode();
				if (currBagNodes.isMember(v) && sigvec[currSigIndex[v]] == TW_WAITING) {
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
}

std::unordered_map<std::pair<TW_SIGNATURE_TYPE, size_t>, u_int64_t, hash_pair> getDPCNTatBag(
		std::unordered_map<TW_SIGNATURE_TYPE,
				std::tuple<size_t, TW_SIGNATURE_TYPE, TW_SIGNATURE_TYPE>>& DPatbag) {
	std::unordered_map<std::pair<TW_SIGNATURE_TYPE, size_t>, u_int64_t, hash_pair> DPCNTatBag;
	for (auto& [sig, val] : DPatbag) {
		DPCNTatBag[{sig, std::get<0>(val)}] = 1;
	}
	return DPCNTatBag;
}

void transformDPCNTatBag(
		std::unordered_map<std::pair<TW_SIGNATURE_TYPE, size_t>, u_int64_t, hash_pair>& DPCNTatBag,
		size_t sigsize) {
	for (int i = 0; i < sigsize; i++) {
		for (auto [sig_cnt, val] : std::vector(DPCNTatBag.begin(),
					 DPCNTatBag.end())) { // copy to avoid concurrent modification
			auto sig = sig_cnt.first;
			auto cnt = sig_cnt.second;
			if (sigAt(sig, i) == TW_WAITING) {
				auto sigset = sigAtSet(sig, i, TW_QM);
				if (DPCNTatBag.find({sigset, cnt}) == DPCNTatBag.end()) {
					DPCNTatBag[{sigset, cnt}] = 0;
				}
				DPCNTatBag[{sigset, cnt}] += val;
			}
		}
	}
}

void transformDPCNTatBagBack(
		std::unordered_map<std::pair<TW_SIGNATURE_TYPE, size_t>, u_int64_t, hash_pair>& DPCNTatBag,
		size_t sigsize) {
	for (int i = 0; i < sigsize; i++) {
		for (auto [sig_cnt, val] : std::vector(DPCNTatBag.begin(),
					 DPCNTatBag.end())) { // copy to avoid concurrent modification
			auto sig = sig_cnt.first;
			auto cnt = sig_cnt.second;
			if (sigAt(sig, i) == TW_QM) {
				auto sigset = sigAtSet(sig, i, TW_WAITING);
				if (DPCNTatBag.find({sigset, cnt}) != DPCNTatBag.end()) {
					auto sub = DPCNTatBag[{sigset, cnt}];
					OGDF_ASSERT(DPCNTatBag[sig_cnt] >= sub);
					DPCNTatBag[sig_cnt] -= sub;
				}
			}
		}
	}
}

void ReductionTreeDecomposition::handleJoinNode(htd::vertex_t curbag) {
	// JOIN NODE
	// log << "Join node " << curbag << std::endl;
	OGDF_ASSERT(decomposition->children(curbag).size() == 2);
	auto bagchild1 = decomposition->childAtPosition(curbag, 0);
	auto bagchild2 = decomposition->childAtPosition(curbag, 1);
	auto DPCNTleft = getDPCNTatBag(DP[bagchild1]);
	auto DPCNTright = getDPCNTatBag(DP[bagchild2]);
	transformDPCNTatBag(DPCNTleft, bag_nodes[bagchild1].size());
	transformDPCNTatBag(DPCNTright, bag_nodes[bagchild2].size());
	std::unordered_map<std::pair<TW_SIGNATURE_TYPE, size_t>, u_int64_t, hash_pair> DPCNTans;
	for (auto [sig_cnt, val] : DPCNTleft) {
		auto sig = sig_cnt.first;
		auto cnt = sig_cnt.second;
		for (int kappa = std::max((int)0, (int)minkappa[bagchild2] - (int)treewidth - 1);
				kappa <= minkappa[bagchild2] + treewidth + 1; kappa++) {
			// TODO: check if right range
			if (DPCNTright.find({sig, kappa}) != DPCNTright.end()) {
				if (DPCNTans.find({sig, cnt + kappa - cntInDs(sig)}) == DPCNTans.end()) {
					DPCNTans[{sig, cnt + kappa - cntInDs(sig)}] = 0;
				}
				DPCNTans[{sig, cnt + kappa - cntInDs(sig)}] += val * DPCNTright[{sig, kappa}];
			}
		}
	}
	transformDPCNTatBagBack(DPCNTans, bag_nodes[bagchild1].size());
	for (auto& [sig_cnt, val] : DPCNTans) {
		if (val == 0) {
			continue;
		}
		minkappa[curbag] = std::min(minkappa[curbag], sig_cnt.second);
		auto sig = sig_cnt.first;
		auto cnt = sig_cnt.second;
		// the solution of children will be computed in the backtracking step
		updateTable(DP[curbag], sig, cnt, 0, 0);
	}
	// log << "join node " << curbag << " done" << std::endl;
	// for (auto& [sig, val] : DP[bagchild1]) {
	// 	auto sigvec = toSigVec(sig, bag_nodes[bagchild1].size());
	// 	for (auto& [sig2, val2] : DP[bagchild2]) {
	// 		auto sigvec2 = toSigVec(sig2, bag_nodes[bagchild2].size());

	// 		size_t cnt_both_in = 0;
	// 		bool same = true;
	// 		for (size_t i = 0; i < sigvec.size(); i++) {
	// 			if (sigvec[i] == TW_INDS && sigvec2[i] == TW_INDS) {
	// 				cnt_both_in++;
	// 			}
	// 			if (sigvec[i] == TW_INDS && sigvec2[i] != TW_INDS) {
	// 				same = false;
	// 			}
	// 			if (sigvec[i] != TW_INDS && sigvec2[i] == TW_INDS) {
	// 				same = false;
	// 			}
	// 			sigvec2[i] = std::max(sigvec[i], sigvec2[i]);
	// 		}
	// 		if (!same) {
	// 			continue;
	// 		}
	// 		OGDF_ASSERT(DP[curbag].find(fromSigVec(sigvec2)) != DP[curbag].end());
	// 		auto curval = std::get<0>(DP[curbag][fromSigVec(sigvec2)]);
	// 		OGDF_ASSERT(curval <= std::get<0>(val) + std::get<0>(val2) - cnt_both_in);
	// 		// updateTable(DP[curbag], fromSigVec(sigvec2),
	// 		// 		std::get<0>(val) + std::get<0>(val2) - cnt_both_in, sig, sig2);
	// 	}
	// }
}

void ReductionTreeDecomposition::handleCopyNode(htd::vertex_t curbag) {
	OGDF_ASSERT(decomposition->children(curbag).size() == 1);
	auto bagchild = decomposition->childAtPosition(curbag, 0);
	minkappa[curbag] = minkappa[bagchild];
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

void ReductionTreeDecomposition::backTrackJoinNode(htd::vertex_t curbag) {
	auto sig = chosensig[curbag];
	auto val = std::get<0>(DP[curbag][sig]);

	auto child1 = decomposition->childAtPosition(curbag, 0);
	auto child2 = decomposition->childAtPosition(curbag, 1);

	// try to find two sigs for child1 and child2 that match sig
	size_t cntdom = cntInDs(sig);
	size_t l = bag_nodes[curbag].size();

	auto sigvec = toSigVec(sig, l);
	std::vector<TW_SIGNATURE_TYPE> sigvecempty(l);
	std::vector<size_t> inds;
	for (size_t i = 0; i < l; i++) {
		if (sigvec[i] != TW_INDS) {
			inds.push_back(i);
		} else if (sigvec[i] == TW_WAITING) {
			sigvecempty[i] = TW_WAITING;
		} else {
			sigvecempty[i] = TW_INDS;
		}
	}

	auto cntbits = l - cntdom;
	size_t to = pow((size_t)2, cntbits);
	for (size_t bitmask = 0; bitmask < to; bitmask++) {
		// log << "bitmask: " << bitmask << std::endl;
		auto sigvecl = sigvecempty;
		auto j = bitmask;
		for (auto ind : inds) {
			sigvecl[ind] = j & 1;
			j >>= 1;
		}
		if (DP[child1].find(fromSigVec(sigvecl)) != DP[child1].end()) {
			std::vector<size_t> indsr;
			auto sigvecemptyr = sigvecempty;
			for (size_t i = 0; i < l; i++) {
				if (sigvec[i] == TW_DOMINATED && sigvecl[i] == TW_WAITING) {
					sigvecemptyr[i] = TW_DOMINATED;
				} else if (sigvec[i] == TW_DOMINATED) {
					indsr.push_back(i);
				}
			}
			size_t tor = pow((size_t)2, indsr.size());
			for (size_t bitmaskr = 0; bitmaskr < tor; bitmaskr++) {
				// log << "bitmaskr: " << bitmaskr << std::endl;
				auto sigvecr = sigvecemptyr;
				auto j = bitmaskr;
				for (auto ind : indsr) {
					sigvecr[ind] = j & 1;
					j >>= 1;
				}
				if (DP[child2].find(fromSigVec(sigvecr)) == DP[child2].end()) {
					continue;
				}
				auto sigl = fromSigVec(sigvecl);
				auto sigr = fromSigVec(sigvecr);
				auto val1 = std::get<0>(DP[child1][sigl]);
				auto val2 = std::get<0>(DP[child2][sigr]);
				if (val == val1 + val2 - cntInDs(sig)) {
					if (chosensig.size() <= std::max(child1, child2)) {
						chosensig.resize(std::max(child1, child2) + 1);
					}
					chosensig[child1] = sigl;
					chosensig[child2] = sigr;
					return;
				}
			}
		}
	}
	OGDF_ASSERT(false);
}

int ReductionTreeDecomposition::solveDPExact() {
	// Algorithm from https://arxiv.org/pdf/1806.01667
	OGDF_ASSERT(decomposition != nullptr);
	htd::PostOrderTreeTraversal traversal;

	traversal.traverse(*decomposition,
			[&](htd::vertex_t curbag, htd::vertex_t bagparent, std::size_t depth) {
				OGDF_ASSERT(decomposition->isVertex(curbag));
				if (DP.size() <= curbag) {
					DP.resize(curbag + 1);
					// DPCNT.resize(curbag + 1);
					bag_nodes.resize(curbag + 1);
					minkappa.resize(curbag + 1);
				}
				minkappa[curbag] = std::numeric_limits<size_t>::max();
				for (auto u : decomposition->bagContent(curbag)) {
					ogdf::node node = idnode[u];
					bag_nodes[curbag].push_back(node);
				}
				std::sort(bag_nodes[curbag].begin(), bag_nodes[curbag].end());

				if (decomposition->isLeaf(curbag)) {
					handleLeaf(curbag);

				} else if (decomposition->isForgetNode(curbag)) {
					handleForgetNode(curbag);
				} else if (decomposition->isIntroduceNode(curbag)) {
					handleIntroduceNode(curbag);

				} else if (decomposition->isJoinNode(curbag)) {
					handleJoinNode(curbag);
				} else {
					// There seem to be nodes with children having the exact same bag
					handleCopyNode(curbag);
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
	// return ans;

	log << "Solved, now backtracking" << std::endl;

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
					for (auto u : decomposition->children(curbag)) {
						if (chosensig.size() <= u) {
							chosensig.resize(u + 1);
						}
						auto t = DP[curbag][chosensig[curbag]];
						chosensig[u] = std::get<1>(t);
					}
				} else if (decomposition->isJoinNode(curbag)) {
					backTrackJoinNode(curbag);
				} else {
					for (auto u : decomposition->children(curbag)) {
						if (chosensig.size() <= u) {
							chosensig.resize(u + 1);
						}
						auto t = DP[curbag][chosensig[curbag]];
						chosensig[u] = std::get<1>(t);
					}
				}
			});

	int cntadded = I.DS.size() - sizebefore;
	OGDF_ASSERT(cntadded == ans);
	return ans;
}
