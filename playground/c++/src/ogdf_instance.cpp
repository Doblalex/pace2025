#include "ogdf_instance.hpp"

#include "ogdf/basic/GraphAttributes.h"
#include "ogdf/basic/GraphSets.h"
#include "ogdf/basic/simple_graph_alg.h"
#include "ogdf/decomposition/BCTree.h"
#include "ogdf/decomposition/StaticSPQRTree.h"
#include "ogdf/energybased/FMMMLayout.h"
#include "ogdf/fileformats/GraphIO.h"
#include "ogdf/layered/SugiyamaLayout.h"
#include "ogdf_solver.hpp"

void Instance::dumpBCTree() {
	std::string stamp = std::to_string(std::chrono::system_clock::now().time_since_epoch().count())
			+ "-" + std::to_string(G.numberOfNodes());
	log << stamp << ".svg" << std::endl;
	std::filesystem::create_directory("out");

	ogdf::BCTree BC(G);
	ogdf::GraphAttributes BCA(BC.bcTree(), ogdf::GraphAttributes::all);
	ogdf::GraphAttributes BA;
	ogdf::NodeSet<> nodes(BC.auxiliaryGraph());
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

		if (BC.numberOfEdges(node) < G.numberOfNodes() / 4) {
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
		ogdf::NodeArray<ogdf::node> nMap(G, nullptr); // can safely be reused for different CCs
		ogdf::EdgeArray<ogdf::edge> eMap(G, nullptr);
		for (int i = 0; i < CC.numberOfCCs(); ++i) {
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
		OGDF_ASSERT(cn->graphOf() == &BC.auxiliaryGraph());
		return cn;
	};

	ogdf::NodeSet<> nodes(BC.auxiliaryGraph());
	ogdf::NodeArray<ogdf::node> nMap(BC.auxiliaryGraph(), nullptr);
	ogdf::EdgeArray<ogdf::edge> eMap(BC.auxiliaryGraph(), nullptr);
	bool changed = false;
	for (auto node : BC.bcTree().nodes) {
		if (node->degree() == 1 && BC.typeOfBNode(node) == ogdf::BCTree::BNodeType::BComp
				&& BC.numberOfNodes(node) < G.numberOfNodes() / 4 && replaced[node] == Replaced::Unchanged) {
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

			// TODO what if cv already dominated / subsumed / replaced[parent] == Replaced::Dominated?
			if (replaced[parent] != Replaced::Unchanged || is_dominated[cv] || is_subsumed[cv]) {
				continue;
			}
			changed = true;
#ifdef OGDF_DEBUG
			log << "Processing " << (is_root ? "root " : "") << "leaf block " << node << " with "
				<< BC.numberOfNodes(node) << " nodes and "
				<< (is_root ? "single child " : "parent ") << parent << "." << std::endl;
			auto& l = logger.lout() << "CV " << node2ID[cv] << " has BC " << parent
									<< " with parent " << BC.parent(parent) << "("
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

			log << "I1: Computing ds(X_B - v) with cut-vertex v already dominated." << std::endl;
			Instance I1; // X_B - v
			nMap.fillWithDefault();
			eMap.fillWithDefault();
			I1.G.insert(nodes, BC.hEdges(node), nMap, eMap);
			I1.initFrom(*this, nodes, BC.hEdges(node), nMap, eMap, original_n, original_e, copy_e);
			{
				ogdf::Logger::Indent _(logger);
				I1.markDominated(nMap[h_cv]);
				reduceAndSolve(I1, depth * 100 + 10);
			}

			log << "I2: Computing normal ds(X_B)." << std::endl;
			Instance I2; // X_B
			nMap.fillWithDefault();
			eMap.fillWithDefault();
			I2.G.insert(nodes, BC.hEdges(node), nMap, eMap);
			I2.initFrom(*this, nodes, BC.hEdges(node), nMap, eMap, original_n, original_e, copy_e);
			{
				ogdf::Logger::Indent _(logger);
				reduceAndSolve(I2, depth * 100 + 20);
			}

			if (I1.DS.size() < I2.DS.size()) {
				log << "RR-BC Case 1: ds(X_B - v) < ds(X_B). Removing block, but leaving cut-vertex unchanged."
					<< std::endl;
				OGDF_ASSERT(I2.DS.size() - I1.DS.size() == 1);
				// replaced[parent] = Replaced::Unchanged;
				addToDominatingSet(I1.DS.begin(), I1.DS.end(), "ds(X_B - v)");

			} else {
				log << "RR-BC Case 2: ds(X_B - v) = ds(X_B). Looking for optimal ds(X_B) containing cut-vertex v..."
					<< std::endl;
				OGDF_ASSERT(I1.DS.size() == I2.DS.size());

				if (std::find(I2.DS.begin(), I2.DS.end(), node2ID[cv]) != I2.DS.end()) {
					log << "RR-BC Case 2.1(sc): The ds(X_B) from I2 contains v. "
						<< "Adding v to DS and removing block." << std::endl;
					replaced[parent] = Replaced::AddToDS;
					addToDominatingSet(I2.DS.begin(), I2.DS.end(), "ds(B_X) containing v");

				} else {
					log << "I3: Computing ds(X_B) containing v." << std::endl;
					Instance I3; // X_B with v forced
					nMap.fillWithDefault();
					eMap.fillWithDefault();
					I3.G.insert(nodes, BC.hEdges(node), nMap, eMap);
					I3.initFrom(*this, nodes, BC.hEdges(node), nMap, eMap, original_n, original_e,
							copy_e);
					{
						ogdf::Logger::Indent _(logger);
						I3.addToDominatingSet(nMap[h_cv]);
						reduceAndSolve(I3, depth * 100 + 30);
					}

					if (I2.DS.size() == I3.DS.size()) {
						log << "RR-BC Case 2.1: Found optimal ds(X_B) containing v. "
							<< "Adding v to DS and removing block." << std::endl;
						replaced[parent] = Replaced::AddToDS;
						addToDominatingSet(I3.DS.begin(), I3.DS.end(), "ds(B_X) containing v");
					} else {
						log << "RR-BC Case 2.2: No optimal ds(X_B) contains v. "
							<< "Marking v as dominated (by B) and removing block B." << std::endl;
						OGDF_ASSERT(I3.DS.size() - I2.DS.size() == 1);
						replaced[parent] = Replaced::MarkDominated;
						addToDominatingSet(I2.DS.begin(), I2.DS.end(), "ds(B_X) not containing v");
					}
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
				OGDF_ASSERT(DS.back() == cv_id);
				DS.pop_back();
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
	auto uwithu = outadjMask[u] | (1ull << (nodehash(u) & 63ull));
	if ((uwithu | outadjMask[v]) != uwithu) {
		return false;
	}
	bool cond = true;
	forAllOutAdj(v, [&](ogdf::adjEntry adj) {
		if (adj->twinNode() != u && !adju[adj->twinNode()]) {
			cond = false;
			return false;
		}
		return true;
	});
	return cond;
}

bool Instance::subsumptionCondition2(const ogdf::node& u, const ogdf::node& v,
		const ogdf::NodeArray<bool>& adju) {
	return is_dominated[v] || adju[v];
}

bool Instance::subsumptionCondition3(const ogdf::node& u, const ogdf::node& v,
		ogdf::NodeArray<bool>& inadjv, const ogdf::NodeArray<u_int64_t>& inadjMask) {
	if (is_dominated[v]) {
		return true;
	}
	auto vwithv = inadjMask[v] | (1ull << (nodehash(v) & 63ull));
	bool filter = true;
	if ((inadjMask[u] | vwithv) != vwithv) {
		filter = false; // TODO needed/unused?
	}
	bool cond = true;
	forAllInAdj(v, [&](ogdf::adjEntry adj) {
		inadjv[adj->twinNode()] = true;
		return true;
	});

	forAllInAdj(u, [&](ogdf::adjEntry adj) {
		if (adj->twinNode() != v && !inadjv[adj->twinNode()]) {
			cond = false;
			return false;
		}
		return true;
	});

	forAllInAdj(v, [&](ogdf::adjEntry adj) {
		inadjv[adj->twinNode()] = false;
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
	ogdf::NodeSet<> couldBeStronglySubsumed(G);
	auto outAdjMask = computeOutadjMask();
	auto inAdjMask = computeInadjMask();
	for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
		OGDF_ASSERT(i >= 0); // protect against endless loops
		auto u = *it;
		OGDF_ASSERT(checkNode(u));
		it++;

		if (is_dominated[u] && is_subsumed[u]) {
			log << "condition 1 applies, deleting " << node2ID[u] << std::endl;
			safeDelete(u, it);
			cnt_removed++;
			continue;
		}

		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			outadju[adj->twinNode()] = true;
			theedge[adj->twinNode()] = adj->theEdge();
			return true;
		});
		couldBeStronglySubsumed.clear();
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
				log << "condition 2 applies, deleting " << node2ID[v] << std::endl;
				safeDelete(v, it);
				cnt_removed++;
				continue;
			}
			if (is_subsumed[u]) {
				if (!is_subsumed[v]) {
					if (subsumptionCondition1(u, v, outadju, outAdjMask)
							&& subsumptionCondition3(u, v, inadjv, inAdjMask)) {
						log << "condition 3 applies, deleting " << node2ID[v] << std::endl;
						safeDelete(v, it);
						cnt_removed++;
					}
				} else {
					if (subsumptionCondition3(u, v, inadjv, inAdjMask)) {
						safeDelete(v, it);
						log << "condition 4 applies, deleting" << node2ID[v] << std::endl;
						cnt_removed++;
					}
				}
			} else {
				if (is_dominated[u] && !is_dominated[v]) {
					if (subsumptionCondition2(u, v, outadju)
							&& subsumptionCondition1(u, v, outadju, outAdjMask)) {
						log << "condition 5 applies, deleting " << node2ID[v]
							<< " and contracting edges" << std::endl;
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
								safeDelete(adj->theEdge());
							}
							return true;
						});
						is_dominated[u] = false;
						safeDelete(v, it);
						cnt_removed++;
					}
				} else if (is_dominated[v] && !is_subsumed[v]) {
					if (subsumptionCondition1(u, v, outadju, outAdjMask)) {
						log << "condition 6 applies, deleting " << node2ID[v] << std::endl;
						safeDelete(v, it);
						cnt_removed++;
					}
				} else if (is_subsumed[v]) {
					if (subsumptionCondition2(u, v, outadju)
							&& subsumptionCondition3(u, v, inadjv, inAdjMask)) {
						log << "condition 7 applies, deleting " << node2ID[v] << std::endl;
						safeDelete(v, it);
						cnt_removed++;
					}
				} else {
					if (subsumptionCondition2(u, v, outadju)
							&& subsumptionCondition1(u, v, outadju, outAdjMask)
							&& subsumptionCondition3(u, v, inadjv, inAdjMask)) {
						log << "condition 8 applies, deleting " << node2ID[v] << std::endl;
						safeDelete(v, it);
						cnt_removed++;
					}
				}
			}
		}
		// does not matter if some v are deleted and not reset here because they will never be accessed again
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			outadju[adj->twinNode()] = false;
			return true;
		});
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
	ogdf::NodeSet<> couldBesubsumed(G);
	for (auto it = G.nodes.begin(); it != G.nodes.end(); --i) {
		OGDF_ASSERT(i >= 0); // protect against endless loops
		auto u = *it;
		OGDF_ASSERT(checkNode(u));
		it++;

		if (is_subsumed[u]) { // if u is subsumed, why should it subsume someone else?
			continue;
		}

		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			outadju[adj->twinNode()] = true;
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
			if (!(is_dominated[u] && !is_dominated[v])) {
				// easy case, u can be taken into dominating set instead of v
				if (subsumptionCondition2(u, v, outadju)
						&& subsumptionCondition1(u, v, outadju, outAdjMask)) {
					// log<<"Marking vertex "<<node2ID[v]<<" as subsumed"<<std::endl;
					markSubsumed(v);
					cnt_subsumed++;
				}
			} else {
				// TODO: what do we do here? Maybe again some contraction and marking u as not dominated (before checking condition 1 and 2?
			}
		}
		forAllOutAdj(u, [&](ogdf::adjEntry adj) {
			outadju[adj->twinNode()] = false;
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
		forAllOutAdj(v, [&](ogdf::adjEntry adj) {
			outadjMask[v] |= (1ull << (nodehash(adj->twinNode()) & 63ull));
			return true;
		});
	}

	return outadjMask;
}

ogdf::NodeArray<u_int64_t> Instance::computeInadjMask() {
	ogdf::NodeArray<u_int64_t> inadjMask(G, 0);
	for (auto v : G.nodes) {
		forAllInAdj(v, [&](ogdf::adjEntry adj) {
			inadjMask[v] |= (1ull << (nodehash(adj->twinNode()) & 63ull));
			return true;
		});
	}
	return inadjMask;
}
