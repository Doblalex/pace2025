#include "ogdf_instance.hpp"

#include "hopcroftcarp.hpp"
#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/basic/GraphSets.h"
#include "ogdf/basic/simple_graph_alg.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/decomposition/StaticSPQRTree.h"
#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/fileformats/GraphIO.h"
#include "ogdf/layered/SugiyamaLayout.h"
#include "ogdf_solver.hpp"
#include "ogdf_subsetrefine.hpp"

void Instance::dumpBCTree() {
	std::string stamp = std::to_string(std::chrono::system_clock::now().time_since_epoch().count())
			+ "-" + std::to_string(G.numberOfNodes());
	log << stamp << ".svg" << std::endl;
	std::filesystem::create_directory("out");

	ogdf::BCTree BC(G);
	ogdf::GraphAttributes BCA(BC.bcTree(), ogdf::GraphAttributes::all);
	ogdf::GraphAttributes BA;
	ogdf::NodeSet nodes(BC.auxiliaryGraph());
	ogdf::NodeArray<ogdf::node> nMap(BC.auxiliaryGraph(), nullptr);
	ogdf::EdgeArray<ogdf::edge> eMap(BC.auxiliaryGraph(), nullptr);

	for (auto node : BC.bcTree().nodes) {
		if (BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp) {
			BCA.shape(node) = ogdf::Shape::Ellipse;
			BCA.label(node) = std::to_string(BC.numberOfNodes(node));
			BCA.width(node) = 15 + 5 * BCA.label(node).size();
		} else {
			BCA.shape(node) = ogdf::Shape::Rhomb;
			BCA.label(node) = "";
			continue;
		}

		if ((float)BC.numberOfEdges(node) <= (float)G.numberOfNodes() * BLOCK_FRACTION) {
			continue;
		}

		nodes.clear();
		for (auto e : BC.hEdges(node)) {
			nodes.insert(e->source());
			nodes.insert(e->target());
		}
		nMap.fillWithDefault();
		eMap.fillWithDefault();
		ogdf::Graph B;
		BA.init(B, ogdf::GraphAttributes::nodeLabel);
		B.insert(nodes, BC.hEdges(node), nMap, eMap);
		for (auto n : nodes) {
			BA.label(nMap[n]) = node2ID[BC.original(n)];
		}
		auto file = "out/block-" + stamp + "-b" + std::to_string(node->index()) + "-"
				+ std::to_string(BC.numberOfNodes(node));
		// ogdf::GraphIO::write(B, file + ".gml");
		ogdf::GraphIO::write(B, file + ".s6"); //

		// SPQR-Tree of B
		{
			ogdf::makeParallelFreeUndirected(B);
			// // Smooth out subdivided edges
			// ogdf::safeForEach(B.nodes,
			//                   [&B](ogdf::node n) {
			//                       if (n->degree() == 2) {
			//                           if (n->outdeg()==2 || n->indeg()==2) {
			//                               B.reverseEdge(n->adjEntries.head()->theEdge());
			//                           }
			//                           B.unsplit(n);
			//                       }
			//                   });
			// // And remove now-parallel edges
			// ogdf::makeParallelFreeUndirected(B);
			ogdf::StaticSPQRTree T(B);
			ogdf::GraphAttributes TA(T.tree(), ogdf::GraphAttributes::all);
			for (auto n : T.tree().nodes) {
				switch (T.typeOf(n)) {
				case ogdf::SPQRTree::NodeType::SNode:
					TA.label(n) = "S " + std::to_string(T.skeleton(n).getGraph().numberOfEdges());
					TA.shape(n) = ogdf::Shape::Ellipse;
					TA.strokeColor(n) = ogdf::Color(ogdf::Color::Name::Darkgreen);
					break;
				case ogdf::SPQRTree::NodeType::PNode:
					TA.label(n) = "P " + std::to_string(T.skeleton(n).getGraph().numberOfEdges());
					TA.shape(n) = ogdf::Shape::Rhomb;
					TA.strokeColor(n) = ogdf::Color(ogdf::Color::Name::Darkblue);
					break;
				case ogdf::SPQRTree::NodeType::RNode:
					TA.label(n) = "R " + std::to_string(T.skeleton(n).getGraph().numberOfEdges())
							+ "/" + std::to_string(T.skeleton(n).getGraph().numberOfNodes());
					TA.shape(n) = ogdf::Shape::Rect;
					TA.strokeColor(n) = ogdf::Color(ogdf::Color::Name::Darkred);
					break;
				}
				TA.strokeWidth(n) = 2;
				TA.width(n) = 15 + TA.label(n).size() * 5;
				TA.height(n) = 20 + TA.label(n).size() * 3;
			}
			ogdf::FMMMLayout().call(TA);
			ogdf::GraphIO::write(TA, file + "-spqr-fmmm.svg");
			ogdf::SugiyamaLayout().call(TA);
			ogdf::GraphIO::write(TA, file + "-spqr-sugi.svg");
		}
	}

	ogdf::FMMMLayout().call(BCA);
	ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".svg");
	// ogdf::GraphIO::write(BCA, "out/fmmm-" + stamp + ".gml");
	ogdf::SugiyamaLayout().call(BCA);
	ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".svg");
	// ogdf::GraphIO::write(BCA, "out/sugi-" + stamp + ".gml");
}

bool Instance::reductionExtremeDegrees() {
	int orig_N = G.numberOfNodes(), orig_DS = DS.size();
	bool reduced = false;
	int i = G.numberOfNodes();
	for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
		OGDF_ASSERT(i >= 0); // protect against endless loops
		auto n = *it;
		OGDF_ASSERT(checkNode(n));
		++it; // increment now so that we can safely delete n

		// isolated
		if (n->indeg() == 0 && !is_dominated[n]) {
			addToDominatingSet(n, it);
			reduced = true;
			continue;
		}
		if (n->outdeg() == 0 && is_dominated[n]) {
			safeDelete(n, it);
			reduced = true;
			continue;
		}

		// antenna
		if (n->indeg() == 1
				&& (n->outdeg() == 0
						|| (n->outdeg() == 1
								&& n->adjEntries.head()->twinNode()
										== n->adjEntries.tail()->twinNode()))) {
			auto u = n->adjEntries.head()->twinNode();
			OGDF_ASSERT(u != n);
			addToDominatingSet(u, it);
			reduced = true;
			continue;
		}
		if (n->outdeg() == 1 && is_dominated[n]) {
			OGDF_ASSERT(n->adjEntries.head()->isSource());
			auto u = n->adjEntries.head()->twinNode();
			if (!is_subsumed[u] || u->indeg() >= 2) {
				safeDelete(n, it);
			} else {
				addToDominatingSet(n, it);
			}
			reduced = true;
			continue;
		}
	}
	if (reduced) {
		log << "Simple reduction removed " << (orig_N - G.numberOfNodes()) << " vertices, added "
			<< (DS.size() - orig_DS) << " to DS." << std::endl;
	}
	return reduced;
}

std::list<Instance> Instance::decomposeConnectedComponents() {
	ogdf::Graph::CCsInfo CC(G);
	std::list<Instance> instances;
	if (CC.numberOfCCs() > 1) {
		ogdf::NodeArray<ogdf::node> nMap(G, nullptr);
		ogdf::EdgeArray<ogdf::edge> eMap(G, nullptr);
		for (int i = 0; i < CC.numberOfCCs(); ++i) {
			if (i > 0) {
				nMap.fillWithDefault();
				// eMap can safely be reused for different CCs
			}
			instances.emplace_back();
			Instance& I = instances.back();
			I.G.insert(CC, i, nMap, eMap);
			I.initFrom(*this, CC.nodes(i), CC.edges(i), nMap, eMap);
		}
	}
	return instances;
}

static int max_BC_changes = -1;

bool Instance::reductionBCTree(int depth) {
	enum class Replaced { Unchanged, MarkDominated, AddToDS };

	if (max_BC_changes == 0) {
		return false;
	}
	int orig_N = G.numberOfNodes(), orig_DS = DS.size();
	ogdf::BCTree BC(G);
	ogdf::NodeArray<Replaced> replaced(BC.bcTree(), Replaced::Unchanged);
	const auto original_n = [&BC](ogdf::node n) {
		OGDF_ASSERT(n->graphOf() == &BC.auxiliaryGraph());
		auto on = BC.original(n);
		OGDF_ASSERT(on->graphOf() == &BC.originalGraph());
		return on;
	};
	const auto original_e = [&BC](ogdf::edge n) {
		OGDF_ASSERT(n->graphOf() == &BC.auxiliaryGraph());
		auto on = BC.original(n);
		OGDF_ASSERT(on->graphOf() == &BC.originalGraph());
		return on;
	};
	const auto copy_e = [&BC](ogdf::edge n) {
		OGDF_ASSERT(n->graphOf() == &BC.originalGraph());
		auto cn = BC.rep(n);
		OGDF_ASSERT(cn == nullptr || cn->graphOf() == &BC.auxiliaryGraph());
		return cn;
	};

	ogdf::NodeSet nodes(BC.auxiliaryGraph());
	ogdf::NodeArray<ogdf::node> nMap(BC.auxiliaryGraph(), nullptr);
	ogdf::EdgeArray<ogdf::edge> eMap(BC.auxiliaryGraph(), nullptr);
	bool changed = false;
	for (auto node : BC.bcTree().nodes) {
		if (node->degree() == 1 && BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp
				&& (float)BC.numberOfNodes(node) <= (float)G.numberOfNodes() * BLOCK_FRACTION
				&& replaced[node] == Replaced::Unchanged) {
			ogdf::node h_cv;
			ogdf::node parent = BC.parent(node);
			bool is_root = false;
			if (parent == nullptr) {
				// deg-1 root
				h_cv = BC.hParNode(node->adjEntries.head()->twinNode());
				parent = node->adjEntries.head()->twinNode();
				is_root = true;
			} else {
				h_cv = BC.hRefNode(node);
			}
			OGDF_ASSERT(parent != nullptr);
			OGDF_ASSERT(BC.typeOfBNode(parent) == ogdf::BCTree::BNodeType::CComp);
			OGDF_ASSERT(h_cv != nullptr);
			ogdf::node cv = BC.original(h_cv);
			OGDF_ASSERT(BC.bcproper(cv) == parent);

			if (replaced[parent] != Replaced::Unchanged) {
				continue;
			}
			changed = true;
#ifdef OGDF_DEBUG
			log << "Processing " << (is_root ? "root " : "") << "leaf block " << node << " with "
				<< BC.numberOfNodes(node) << " nodes and "
				<< (is_root ? "single child " : "parent ") << parent << "." << std::endl;
			auto& l = logger.lout()
					<< "CV " << node2ID[cv] << " is "
					<< (is_dominated[cv] ? "already dominated"
										 : (is_subsumed[cv] ? "subsumed" : "unmarked"))
					<< " and has BC " << parent << " with parent " << BC.parent(parent) << "("
					<< BC.numberOfEdges(BC.parent(parent)) << ")" << " and children";
			for (auto adj : parent->adjEntries) {
				if (!adj->isSource()) {
					l << " " << adj->twinNode() << "(" << BC.numberOfEdges(adj->twinNode()) << ")";
				}
			}
			l << std::endl;
#endif

			nodes.clear();
			for (auto e : BC.hEdges(node)) {
				nodes.insert(e->source());
				nodes.insert(e->target());
			}
#ifdef OGDF_DEBUG
			OGDF_ASSERT(nodes.isMember(h_cv));
			for (auto n : nodes) {
				checkNode(BC.original(n));
			}
#endif

			bool en_case_1 = true, en_case_2A = true, en_case_2B = true;
			if (is_dominated[cv]) {
				OGDF_ASSERT(cv->indeg() == 0); // B no in, R no in
				// 2A or 2B
				log << "CV is dominated, disabling case 1" << std::endl;
				en_case_1 = false;
			} else if (is_subsumed[cv]) {
				OGDF_ASSERT(cv->outdeg() == 0); // B no out, R no out
				// 1 or 2B
				log << "CV is subsumed, disabling case 2A" << std::endl;
				en_case_2A = false;
			} else if (h_cv->indeg() == 0) { // B no in
				// 1 or 2A
				log << "No edges from B into CV, disabling case 2B" << std::endl;
				en_case_2B = false;
			} else if (h_cv->outdeg() == 0) { // B no out
				// 1 or 2B
				log << "No edges from CV out to B, disabling case 2A" << std::endl;
				en_case_2A = false;
			}

#define LAZY_INSTANCE(INST, COM, INIT, D)                                                          \
	std::unique_ptr<Instance> INST;                                                                \
	auto get_##INST = [&]() -> Instance& {                                                         \
		if (INST) {                                                                                \
			return *INST;                                                                          \
		}                                                                                          \
		log << COM << std::endl;                                                                   \
		INST = std::make_unique<Instance>();                                                       \
		nMap.fillWithDefault();                                                                    \
		eMap.fillWithDefault();                                                                    \
		INST->G.insert(nodes, BC.hEdges(node), nMap, eMap);                                        \
		INST->initFrom(*this, nodes, BC.hEdges(node), nMap, eMap, original_n, original_e, copy_e); \
		ogdf::Logger::Indent _(logger);                                                            \
		INIT;                                                                                      \
		reduceAndSolve(*INST, depth * 100 + D);                                                    \
		return *INST;                                                                              \
	};

			LAZY_INSTANCE(I1, "I1: Computing ds(X_B - v) with cut-vertex v already dominated.",
					I1->markDominated(nMap[h_cv]), 10);
			LAZY_INSTANCE(I2, "I2: Computing normal ds(X_B).", , 20);
			LAZY_INSTANCE(I3, "I3: Computing ds(X_B) containing v.",
					I3->addToDominatingSet(nMap[h_cv]), 30);
			auto smaller_DS_no_CV = [&] { return get_I1().DS.size() < get_I2().DS.size(); };
			std::unordered_set<int>* opt_DS_with_CV = nullptr;
			auto opt_DS_has_CV = [&] {
				if (opt_DS_with_CV) {
					return true;
				}
				if (get_I2().DS.find(node2ID[cv]) != get_I2().DS.end()) {
					log << "RR-BC Case 2A(sc): The ds(X_B) from I2 already contains v." << std::endl;
					opt_DS_with_CV = &get_I2().DS;
					return true;
				}
				if (get_I2().DS.size() == get_I3().DS.size()) {
					opt_DS_with_CV = &get_I3().DS;
					return true;
				} else {
					return false;
				}
			};

			if (en_case_1 && smaller_DS_no_CV()) {
				log << "RR-BC Case 1: ds(X_B - v) < ds(X_B). Removing block, but leaving cut-vertex unchanged."
					<< std::endl;
				OGDF_ASSERT(get_I2().DS.size() - get_I1().DS.size() == 1);
				addToDominatingSet(get_I1().DS.begin(), get_I1().DS.end(), "ds(X_B - v)");
			} else {
				log << "RR-BC Case 2: ds(X_B - v) = ds(X_B). Looking for optimal ds(X_B) containing cut-vertex v..."
					<< std::endl;
				OGDF_ASSERT(get_I1().DS.size() == get_I2().DS.size());

				if (en_case_2A && (!en_case_2B || opt_DS_has_CV())) {
					log << "RR-BC Case 2A: Found optimal ds(X_B) containing v. "
						<< "Adding v to DS and removing block." << std::endl;
					replaced[parent] = Replaced::AddToDS;
					opt_DS_has_CV(); // ensure opt_DS_with_CV is set (this has no overhead in all cases)
					OGDF_ASSERT(opt_DS_with_CV != nullptr);
					addToDominatingSet(opt_DS_with_CV->begin(), opt_DS_with_CV->end(),
							"ds(B_X) containing v");
				} else {
					OGDF_ASSERT(en_case_2B);
					log << "RR-BC Case 2B: No optimal ds(X_B) contains v. "
						<< "Marking v as dominated (by B) and removing block B." << std::endl;
					OGDF_ASSERT(is_subsumed[cv] || get_I3().DS.size() - get_I2().DS.size() == 1);
					replaced[parent] = Replaced::MarkDominated;
					addToDominatingSet(get_I2().DS.begin(), get_I2().DS.end(),
							"ds(B_X) not containing v");
				}
			}

			// delete all vertices in B now, handle the CV later
			for (auto hn : nodes) {
				auto n = BC.original(hn);
				if (n != cv) {
					safeDelete(n);
				}
			}

			if (max_BC_changes > 0) {
				max_BC_changes--;
			}
			if (max_BC_changes == 0) {
				break;
			}
		}
	}
	int d = 0, r = 0;
	if (changed) {
		log << "Now processing marked cut-vertices..." << std::endl;
		for (auto node : BC.bcTree().nodes) {
			if (replaced[node] == Replaced::Unchanged) {
				continue;
			}
			if (BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp) {
				continue;
			}
			auto cv = BC.original(BC.hRefNode(node));

			if (replaced[node] == Replaced::MarkDominated) {
				markDominated(cv);
				++d;
			} else {
				OGDF_ASSERT(replaced[node] == Replaced::AddToDS);
				auto cv_id = node2ID[cv];
				addToDominatingSet(cv);
				// OGDF_ASSERT(DS.back() == cv_id);
				// DS.erase(cv_id);
				++r;
			}
		}
		log << "BCTree reduction removed " << r << " cut-vertices and marked " << d
			<< " as dominated." << std::endl;
		log << "In total, removed " << (orig_N - G.numberOfNodes()) << " vertices and added "
			<< (DS.size() - orig_DS) << " to DS." << std::endl;
	}
	return changed;
}

bool Instance::subsumptionCondition1(const ogdf::node& u, const ogdf::node& v,
		const ogdf::NodeArray<bool>& adju, const ogdf::NodeArray<u_int64_t>& outadjMask) {
	if (is_subsumed[v]) {
		return true;
	}
	if ((outadjMask[u] | outadjMask[v]) != outadjMask[u]) {
		return false;
	}
	bool cond = true;
	forAllCanDominate(v, [&](ogdf::node adj) {
		if (!adju[adj]) {
			cond = false;
			return false;
		}
		return true;
	});
	return cond;
}

bool Instance::subsumptionCondition2(const ogdf::node& u, const ogdf::node& v,
		ogdf::NodeArray<bool>& inadjv, const ogdf::NodeArray<u_int64_t>& inadjMask) {
	if (is_dominated[u] || is_dominated[v]) {
		return true;
	}
	if ((inadjMask[u] | inadjMask[v]) != inadjMask[v]) {
		return false;
	}
	bool cond = true;
	forAllCanBeDominatedBy(v, [&](ogdf::node adj) {
		inadjv[adj] = true;
		return true;
	});

	forAllCanBeDominatedBy(u, [&](ogdf::node adj) {
		if (!inadjv[adj]) {
			cond = false;
			return false;
		}
		return true;
	});

	forAllCanBeDominatedBy(v, [&](ogdf::node adj) {
		inadjv[adj] = false;
		return true;
	});
	return cond;
}

bool Instance::reductionStrongSubsumption() {
	int cnt_removed = 0;
	int i = G.numberOfNodes();
	ogdf::NodeArray<bool> outadju(G, false);
	ogdf::NodeArray<bool> inadjv(G, false);
	ogdf::NodeArray<ogdf::edge> theedge(G, nullptr);
	auto outAdjMask = computeOutadjMask();
	auto inAdjMask = computeInadjMask();
	for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
		OGDF_ASSERT(i >= 0); // protect against endless loops
		auto u = *it;
		OGDF_ASSERT(checkNode(u));
		it++;

		if (is_dominated[u] && is_subsumed[u]) {
			log << "u is already dominated and subsumed, deleting " << node2ID[u] << std::endl;
			safeDelete(u, it);
			cnt_removed++;
			continue;
		}

		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			outadju[adj->twinNode()] = true;
			theedge[adj->twinNode()] = adj->theEdge();
			return true;
		});
		if (!is_dominated[u] && !is_subsumed[u]) {
			outadju[u] = true;
		}
		ogdf::NodeSet couldBeStronglySubsumed(G);
		// compute all vertices that share an in-neighbor or out-neighbor with u
		forAllInAdj(u, [&](ogdf::adjEntry adj) {
			forAllOutAdj(adj->twinNode(), [&](ogdf::adjEntry adj2) {
				if (adj2->twinNode() == u) {
					return true;
				}
				couldBeStronglySubsumed.insert(adj2->twinNode());
				return true;
			});
			return true;
		});
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			forAllInAdj(adj->twinNode(), [&](ogdf::adjEntry adj2) {
				if (adj2->twinNode() == u) {
					return true;
				}
				couldBeStronglySubsumed.insert(adj2->twinNode());
				return true;
			});
			return true;
		});
		for (auto it2 = couldBeStronglySubsumed.begin(); it2 != couldBeStronglySubsumed.end();) {
			auto v = *it2;
			++it2;
			if (is_dominated[v] && is_subsumed[v]) {
				log << "v is already dominated and subsumed, deleting " << node2ID[v] << std::endl;
				safeDelete(v, it);
				cnt_removed++;
				continue;
			} else if (subsumptionCondition1(u, v, outadju, outAdjMask)
					&& subsumptionCondition2(u, v, inadjv, inAdjMask)) {
				bool udominated = is_dominated[u];
				bool vdominated = is_dominated[v];

				if (udominated && !vdominated) {
					log << "strong subsumption with u dominated (contraction), deleting "
						<< node2ID[v] << std::endl;
					log << countCanBeDominatedBy(u) << std::endl;
					log << countCanBeDominatedBy(v) << std::endl;
					forAllInAdj(v, [&](ogdf::adjEntry adj) {
						if (adj->twinNode() != u) {
							auto e = G.newEdge(adj->twinNode(), ogdf::Direction::before, u,
									ogdf::Direction::after);
							if (outadju[adj->twinNode()]) {
								reverse_edge[e] = theedge[adj->twinNode()];
								reverse_edge[theedge[adj->twinNode()]] = e;
							} else {
								reverse_edge[e] = nullptr;
							}
						}
						return true;
					});
				} else {
					log << "strong subsumption, deleting " << node2ID[v] << std::endl;
					log << "u: (" << is_dominated[u] << ", " << is_subsumed[u] << ") v: ("
						<< is_dominated[v] << ", " << is_subsumed[v] << ")" << std::endl;
				}
				is_dominated[u] = is_dominated[u] && is_dominated[v];
				is_subsumed[u] = is_subsumed[u] && is_subsumed[v];
				safeDelete(v, it);
				cnt_removed++;
			}
		}
		// does not matter if some v are deleted and not reset here because they will never be accessed again
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			outadju[adj->twinNode()] = false;
			return true;
		});
		outadju[u] = false;
	}
	if (cnt_removed > 0) {
		log << "Strong subsumption removed " << cnt_removed << " vertices" << std::endl;
	}
	return cnt_removed > 0;
}

bool Instance::reductionSubsumption() {
	int cnt_subsumed = 0;
	int i = G.numberOfNodes();
	ogdf::NodeArray<bool> outadju(G, false);
	auto outAdjMask = computeOutadjMask();
	ogdf::NodeSet couldBesubsumed(G);
	for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
		OGDF_ASSERT(i >= 0); // protect against endless loops
		auto u = *it;
		OGDF_ASSERT(checkNode(u));
		it++;

		if (is_subsumed[u]) { // if u is subsumed, why should it subsume someone else?
			continue;
		}

		forAllCanDominate(u, [&](ogdf::node adj) {
			outadju[adj] = true;
			return true;
		});
		// compute all vertices that share an out-neighbor with u
		couldBesubsumed.clear();
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			forAllInAdj(adj->twinNode(), [&](ogdf::adjEntry adj2) {
				if (adj2->twinNode() == u) {
					return true;
				}
				couldBesubsumed.insert(adj2->twinNode());
				return true;
			});
			return true;
		});
		for (auto it2 = couldBesubsumed.begin(); it2 != couldBesubsumed.end();) {
			auto v = *it2;
			++it2;
			if (is_subsumed[v]) {
				continue; // v is already subsumed
			}
			// neither u nor v are subsumed here
			if (subsumptionCondition1(u, v, outadju, outAdjMask)) {
				// log<<"Marking vertex "<<node2ID[v]<<" as subsumed"<<std::endl;
				markSubsumed(v);
				cnt_subsumed++;
			}
		}
		forAllCanDominate(u, [&](ogdf::node adj) {
			outadju[adj] = false;
			return true;
		});
	}
	if (cnt_subsumed > 0) {
		log << "Subsumption subsumed " << cnt_subsumed << " vertices" << std::endl;
	}
	return cnt_subsumed > 0;
}

ogdf::NodeArray<u_int64_t> Instance::computeOutadjMask() {
	ogdf::NodeArray<u_int64_t> outadjMask(G, 0);
	for (auto v : G.nodes) {
		forAllCanDominate(v, [&](ogdf::node adj) {
			outadjMask[v] |= (1ull << (nodehash(adj) & 63ull));
			return true;
		});
	}

	return outadjMask;
}

ogdf::NodeArray<u_int64_t> Instance::computeInadjMask() {
	ogdf::NodeArray<u_int64_t> inadjMask(G, 0);
	for (auto v : G.nodes) {
		forAllCanBeDominatedBy(v, [&](ogdf::node adj) {
			inadjMask[v] |= (1ull << (nodehash(adj) & 63ull));
			return true;
		});
	}
	return inadjMask;
}

bool Instance::reductionNeighborhoodSubsets() {
	SubsetRefine refineSubsume(*this, RefineType::Subsume);
	refineSubsume.init();
	size_t cntsubsumed = refineSubsume.doRefinementReduction();
	SubsetRefine refineDominate(*this, RefineType::Dominate);
	refineDominate.init();
	size_t cntdominated = refineDominate.doRefinementReduction();
	if (cntdominated + cntsubsumed > 0) {
		log << "Dominated " << cntdominated << " nodes" << std::endl;
		log << "Subsumed " << cntsubsumed << " nodes" << std::endl;
		size_t cnt_removed = 0;
		for (auto it = G.nodes.begin(); it != G.nodes.end();) {
			auto v = *it;
			it++;
			if (is_subsumed[v] && is_dominated[v]) {
				safeDelete(v, it);
				cnt_removed++;
			}
		}
		log << "Removed " << cnt_removed << " nodes" << std::endl;
		return true;
	}
	return false;
}

bool Instance::reductionContraction() {
	size_t cnt_removed = 0;
	ogdf::NodeArray<ogdf::edge> theedge(G, nullptr);
	for (auto it = G.nodes.begin(); it != G.nodes.end();) {
		auto u = *it;
		it++;

		if (is_dominated[u] && !is_subsumed[u]) {
			bool contracted = false;

			forAllOutAdj(u, [&](ogdf::adjEntry adj) {
				theedge[adj->twinNode()] = adj->theEdge();
				return true;
			});

			ogdf::node takev = nullptr;
			forAllOutAdj(u, [&](ogdf::adjEntry adj) {
				auto v = adj->twinNode();
				if (!is_dominated[v] && is_subsumed[v]) {
					takev = v;
					return false;
				}
				return true;
			});
			if (takev != nullptr) {
				std::list<ogdf::node> in_neighbors;
				forAllInAdj(takev, [&](ogdf::adjEntry adj2) {
					if (adj2->twinNode() != u) {
						in_neighbors.push_back(adj2->twinNode());
					}
					return true;
				});
				for (auto w : in_neighbors) {
					auto e = G.newEdge(w, ogdf::Direction::before, u, ogdf::Direction::after);
					if (theedge[w] != nullptr) {
						reverse_edge[e] = theedge[w];
						reverse_edge[theedge[w]] = e;
					} else {
						reverse_edge[e] = nullptr;
					}
				}
				safeDelete(takev, it);
				cnt_removed++;
				is_dominated[u] = false;
				removeHiddenIncomingEdges(u);
			}

			forAllOutAdj(u, [&](ogdf::adjEntry adj) {
				theedge[adj->twinNode()] = nullptr;
				return true;
			});
		}
	}
	if (cnt_removed > 0) {
		log << "Contraction removed " << cnt_removed << " vertices" << std::endl;
	}
	return cnt_removed > 0;
}

bool Instance::reductionSpecial1() {
	// generalisation of RR 6 of https://www.sciencedirect.com/science/article/pii/S0166218X11002393
	size_t applications = 0;
	ogdf::NodeArray<u_int32_t> nodeinBipGraph(G, 0);
	for (auto it = G.nodes.begin(); it != G.nodes.end();) {
		auto Sv = *it;
		it++;
		std::vector<ogdf::node> freq2neighs;
		forAllCanDominate(Sv, [&](ogdf::node adj) {
			if (countCanBeDominatedBy(adj) == 2) {
				freq2neighs.push_back(adj);
			}
			return true;
		});
		if (freq2neighs.size() > 1) {
			u_int32_t max = 0;
			u_int32_t indexf = 0;
			std::vector<std::pair<u_int32_t, u_int32_t>> edges;
			for (auto n : freq2neighs) {
				indexf++;
				forAllCanBeDominatedBy(n, [&](ogdf::node Siv) {
					if (Siv != Sv) {
						forAllCanDominate(Siv, [&](ogdf::node adj) {
							if (adj != n) {
								if (nodeinBipGraph[adj] == 0) {
									nodeinBipGraph[adj] = ++max;
								}
								edges.push_back({indexf, nodeinBipGraph[adj]});
							}
							return true;
						});
					}
					return true;
				});
			}
			for (auto n : freq2neighs) {
				forAllCanBeDominatedBy(n, [&](ogdf::node Siv) {
					if (Siv != Sv) {
						forAllCanDominate(Siv, [&](ogdf::node adj) {
							nodeinBipGraph[adj] = 0;
							return true;
						});
					}
					return true;
				});
			}
			BipartiteGraph B(freq2neighs.size(), max);
			for (auto& e : edges) {
				B.add_edge(e.first, e.second);
			}
			auto matching_size = B.hopcroftKarp_algorithm();
			if (matching_size < freq2neighs.size()) {
				applications++;
				addToDominatingSet(Sv);
			} else {
				OGDF_ASSERT(matching_size == freq2neighs.size());
			}
		}
	}
	if (applications > 0) {
		log << "Special reduction rule 1 can dominate " << applications << " vertices" << std::endl;
	}

	return applications > 0;
}

bool Instance::reductionSpecial2(int d) {
	for (auto it = G.nodes.begin(); it != G.nodes.end();) {
		auto Rv = *it;
		it++;
		if (countCanDominate(Rv) == 2) {
			std::vector<ogdf::node> freq2neighs;
			forAllCanDominate(Rv, [&](ogdf::node adj) {
				if (countCanBeDominatedBy(adj) == 2) {
					freq2neighs.push_back(adj);
				}
				return true;
			});
			if (freq2neighs.size() == 2) {
				std::unordered_set<ogdf::node> Q;
				std::vector<ogdf::node> Rivs;
				for (auto n : freq2neighs) {
					forAllCanBeDominatedBy(n, [&](ogdf::node Riv) {
						if (Riv != Rv) {
							Rivs.push_back(Riv);
							forAllCanDominate(Riv, [&](ogdf::node adj) {
								if (adj != n) {
									Q.insert(adj);
								}
								return true;
							});
						}
						return true;
					});
				}
				log << "** Applying special reduction rule 2" << std::endl;
				markSubsumed(Rv);
				markSubsumed(Rivs[0]);
				markSubsumed(Rivs[1]);
				markDominated(freq2neighs[0]);
				markDominated(freq2neighs[1]);
				OGDF_ASSERT(Rivs.size() == 2);
				auto Rvid = node2ID[Rv];
				auto Rv1id = node2ID[Rivs[0]];
				auto Rv2id = node2ID[Rivs[1]];

				auto qnode = G.newNode();
				node2ID[qnode] = ++maxid;
				auto qid = maxid;
				markDominated(qnode);
				for (auto q : Q) {
					auto e = G.newEdge(qnode, ogdf::Direction::before, q, ogdf::Direction::after);
					reverse_edge[e] = nullptr;
				}
				Q.clear();
				reduceAndSolve(*this, d);
				if (DS.find(qid) != DS.end()) {
					DS.insert(Rv1id);
					DS.insert(Rv2id);
					DS.erase(qid);
				} else {
					DS.insert(Rvid);
				}
				return true;
			}
		}
	}
	return false;
}