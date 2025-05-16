#pragma once

#include <unordered_set>

#include "ogdf_util.hpp"

struct Instance {
private:
	bool subsumptionCondition1(const ogdf::node& u, const ogdf::node& v,
			const ogdf::NodeArray<bool>& adju, const ogdf::NodeArray<u_int64_t>& outadjMask);
	bool subsumptionCondition2(const ogdf::node& u, const ogdf::node& v,
			ogdf::NodeArray<bool>& inadjv, const ogdf::NodeArray<u_int64_t>& inadjMask);
	ogdf::NodeArray<u_int64_t> computeOutadjMask();
	ogdf::NodeArray<u_int64_t> computeInadjMask();
	void read_DS(std::istream& is, std::vector<ogdf::node>& ID2node, unsigned int n, unsigned int m);
	void read_HS(std::istream& is, std::vector<ogdf::node>& ID2node, unsigned int n, unsigned int m);

public:
	ogdf::Graph G;
	// std::vector<ogdf::node> ID2node;
	ogdf::NodeArray<int> node2ID;
	std::unordered_set<int> DS;
	ogdf::NodeArray<bool> is_dominated;
	ogdf::NodeArray<bool> is_subsumed;
	ogdf::NodeArray<bool> is_hidden_loop;
	ogdf::EdgeArray<ogdf::edge> reverse_edge;
	ogdf::Graph::DynamicHiddenEdgeSet hidden_edges;
	std::hash<ogdf::node> nodehash;
	size_t maxid;
	std::string type;

	Instance()
		: node2ID(G, -1)
		, is_dominated(G, false)
		, is_hidden_loop(G, false)
		, is_subsumed(G, false)
		, reverse_edge(G, nullptr)
		, hidden_edges(G) { }

	OGDF_NO_COPY(Instance)
	OGDF_NO_MOVE(Instance)

	void clear() {
		hidden_edges.restore();
		G.clear();
		// DS.clear(); Do not clear DS! This value is still used after computing connected components
		node2ID.init(G, -1);
		is_dominated.init(G, false);
		is_subsumed.init(G, false);
		is_hidden_loop.init(G, false);
		reverse_edge.init(G, nullptr);
	}

	bool checkNode(ogdf::node n) {
		OGDF_ASSERT(!is_dominated[n] || n->indeg() == 0);
		OGDF_ASSERT(!is_subsumed[n] || n->outdeg() == 0);
		OGDF_ASSERT(!is_subsumed[n] || n->indeg() > 0 || is_dominated[n]);
		OGDF_ASSERT(n->outdeg() == 0 || n->adjEntries.head()->isSource());
		OGDF_ASSERT(n->indeg() == 0 || !n->adjEntries.tail()->isSource());
#ifdef OGDF_HEAVY_DEBUG
		size_t c = 0;
		for (auto adj_it = n->adjEntries.begin(); adj_it != n->adjEntries.end(); ++c) {
			auto adj = *adj_it;
			++adj_it;
			OGDF_ASSERT(adj->isSource() == (c < n->outdeg()));
		}
#endif
		return true;
	}

	template<typename NL, typename EL>
	void initFrom(const Instance& other, const NL& nodes, const EL& edges,
			const ogdf::NodeArray<ogdf::node>& nMap, const ogdf::EdgeArray<ogdf::edge>& eMap,
			std::function<ogdf::node(ogdf::node)> originalN = internal::idn,
			std::function<ogdf::edge(ogdf::edge)> originalE = internal::ide,
			std::function<ogdf::edge(ogdf::edge)> copyE = internal::ide) {
		for (auto n : nodes) {
			auto tn = nMap[n];
			auto on = originalN(n);
			// other.hidden_edges.adjEntries(n) TODO
			node2ID[tn] = other.node2ID[on];
			is_dominated[tn] = other.is_dominated[on];
			is_subsumed[tn] = other.is_subsumed[on];
			is_hidden_loop[tn] = other.is_hidden_loop[on];

			// embedding breaks when inserting by edge list, so fix it
			size_t o = 0, i = 0;
			const auto& adjs = tn->adjEntries;
			for (auto adj_it = adjs.begin(); o < tn->outdeg();) {
				OGDF_ASSERT(adj_it != adjs.end());
				auto adj = *adj_it;
				++adj_it;
				if (adj->isSource()) {
					++o;
					if (i > 0) {
						G.moveAdj(adj, ogdf::Direction::before, adjs.head());
					}
				} else {
					++i;
				}
			}

			for (auto adj : other.hidden_edges.adjEntries(on)) {
				if (!adj->isSource()) {
					continue;
				}
				auto ce = copyE(adj->theEdge());
				if (ce == nullptr) {
					continue;
				}
				auto src = nMap[ce->source()];
				auto tgt = nMap[ce->target()];
				if (src != nullptr && tgt != nullptr) {
					hidden_edges.hide(G.newEdge(src, tgt));
				}
			}

#if 0
            if (!(!is_dominated[tn] || tn->indeg() == 0)) {
                {
                    auto& l = logger.lout() << "copied edges:";
                    for (auto e : edges) {
                        l<<" "<<eMap[e]<<"["<<e<<"]";
                    }
                    l<<std::endl;
                }
                {
                    auto& l = logger.lout() << "nodes' edges:";
                    for (auto adj : tn->adjEntries) {
                        l<<" "<<adj->theEdge()<<"["<<adj<<"]";
                    }
                    l<<std::endl;
                }
            }
#endif

			OGDF_ASSERT(checkNode(tn));
		}
		for (auto e : edges) {
			auto ore = other.reverse_edge[originalE(e)];
			if (ore != nullptr) {
				reverse_edge[eMap[e]] = eMap[copyE(ore)];
			} else {
				reverse_edge[eMap[e]] = nullptr;
			}
		}
		this->maxid = other.maxid;
	}

	void read(std::istream& is) {
		std::vector<ogdf::node> ID2node;
		read(is, ID2node);
	}

	void read(std::istream& is, std::vector<ogdf::node>& ID2node);

	void dumpBCTree();

	void safeDelete(ogdf::node n, ogdf::Graph::node_iterator& it) {
		if (*it == n) {
			++it; // no need to check for end, as *end == nullptr
		}
		safeDelete(n);
	}

	void safeDelete(ogdf::node n) {
		// ID2node[node2ID[n]] = nullptr;
		logd << "\tsafe delete " << node2ID[n] << std::endl;
		G.delNode(n);
	}

	void safeDelete(ogdf::edge e) {
		ogdf::edge r = reverse_edge[e];
		if (r != nullptr) {
			OGDF_ASSERT(reverse_edge[r] == e);
			reverse_edge[r] = nullptr;
		}
		G.delEdge(e);
	}

	void addToDominatingSet(ogdf::node v) {
		ogdf::Graph::node_iterator it;
		addToDominatingSet(v, it);
	}

	template<typename IT>
	void addToDominatingSet(IT begin, IT end, std::string comment = "") {
#ifdef OGDF_DEBUG
		int before = DS.size();
		auto& l = logger.lout(ogdf::Logger::Level::Minor)
				<< "Insert " << comment << (comment.empty() ? "" : " ") << "into DS:";
		for (auto it = begin; it != end; ++it) {
			l << " " << *it;
		}
		l << std::endl;
		std::copy(begin, end, std::inserter(DS, DS.end()));
		log << "Updated DS" << (comment.empty() ? "" : " with ") << comment << ": " << before << "+"
			<< (DS.size() - before) << "=" << DS.size() << std::endl;
#else
		std::copy(begin, end, std::inserter(DS, DS.end()));
#endif
	}

	// XXX not using any shortcuts here, as they often hurt in other places (like storing the global universal_in vertex)
	void addToDominatingSet(ogdf::node v, ogdf::Graph::node_iterator& it) {
		OGDF_ASSERT(!is_subsumed[v]);
		logd << "\tadd to DS " << node2ID[v] << std::endl;
		DS.insert(node2ID[v]);
		forAllOutAdj(v, [&](ogdf::adjEntry adj) {
			markDominated(adj->twinNode());
			return true;
		});
		ogdf::safeForEach(hidden_edges.adjEntries(v), [&](ogdf::adjEntry adj) {
			if (adj->isSource()) { // now that the successors are really dominated, they cannot have hidden incoming edges anymore!
				markDominated(adj->twinNode());
			}
		});
		{
			auto gll = logger.globalLogLevel();
			logger.globalLogLevel(ogdf::Logger::Level::Alarm);
			safeDelete(v, it);
			logger.globalLogLevel(gll);
		}
	}

	void markDominated(ogdf::node v, bool byreduction = false) {
		if (!is_dominated[v] && !is_subsumed[v] && byreduction) {
			is_hidden_loop[v] = true;
		} else {
			is_hidden_loop[v] = false;
		}
		is_dominated[v] = true;
		forAllInAdj(v, [&](ogdf::adjEntry adj) {
			ogdf::edge e = adj->theEdge();
			ogdf::edge r = reverse_edge[e];
			if (r != nullptr) {
				OGDF_ASSERT(reverse_edge[r] == e);
				reverse_edge[r] = nullptr;
			}
			if (byreduction) {
				hidden_edges.hide(e);
			} else {
				G.delEdge(e);
			}
			return true;
		});
		if (!byreduction) {
			// v is really dominated, so we need to remove incoming hidden edges
			removeHiddenIncomingEdges(v);
		}
	}

	void markSubsumed(ogdf::node v) {
		is_hidden_loop[v] = false;
		is_subsumed[v] = true;
		forAllOutAdj(v, [&](ogdf::adjEntry adj) {
			safeDelete(adj->theEdge());
			return true;
		});
		removeHiddenOutgoingEdges(v);
	}

	void removeHiddenEdges(ogdf::node n) {
		is_hidden_loop[n] = false;
		ogdf::safeForEach(hidden_edges.adjEntries(n), [&](ogdf::adjEntry adj) {
			hidden_edges.restore(adj->theEdge());
			G.delEdge(adj->theEdge());
		});
	}

	void removeHiddenIncomingEdges(ogdf::node n) {
		is_hidden_loop[n] = false;
		ogdf::safeForEach(hidden_edges.adjEntries(n), [&](ogdf::adjEntry adj) {
			if (!adj->isSource()) {
				hidden_edges.restore(adj->theEdge());
				G.delEdge(adj->theEdge());
			}
		});
	}

	void removeHiddenOutgoingEdges(ogdf::node n) {
		is_hidden_loop[n] = false;
		ogdf::safeForEach(hidden_edges.adjEntries(n), [&](ogdf::adjEntry adj) {
			if (adj->isSource()) {
				hidden_edges.restore(adj->theEdge());
				G.delEdge(adj->theEdge());
			}
		});
	}

	size_t countCanDominate(ogdf::node v) {
		if (!is_subsumed[v] && !is_dominated[v]) {
			return v->outdeg() + 1;
		} else {
			return v->outdeg();
		}
	}

	size_t countCanBeDominatedBy(ogdf::node v) {
		if (!is_subsumed[v] && !is_dominated[v]) {
			return v->indeg() + 1;
		} else {
			return v->indeg();
		}
	}

	inline bool forAllCanDominate(ogdf::node v, std::function<bool(ogdf::node)> f) {
		bool ret = forAllOutAdj(v, [&](ogdf::adjEntry adj) { return f(adj->twinNode()); });
		if (ret && !is_subsumed[v] && !is_dominated[v]) {
			return f(v);
		} else {
			return ret;
		}
	}

	inline bool forAllCanBeDominatedBy(ogdf::node v, std::function<bool(ogdf::node)> f) {
		bool ret = forAllInAdj(v, [&](ogdf::adjEntry adj) { return f(adj->twinNode()); });
		if (ret && !is_subsumed[v] && !is_dominated[v]) {
			return f(v);
		} else {
			return ret;
		}
	}

	bool reductionExtremeDegrees();

	bool reductionStrongSubsumption();

	bool reductionSubsumption();

	bool reductionNeighborhoodSubsets();

	bool reductionContraction();

	std::list<Instance> decomposeConnectedComponents();

	bool reductionBCTree(int depth = 0);

	bool reductionSpecial1();

	bool reductionSpecial2(int d);

	std::pair<size_t, size_t> dominationStats() {
		size_t can = 0, needs = 0;
		for (ogdf::node n : G.nodes) {
			if (!is_dominated[n]) {
				++needs;
			}
			if (!is_subsumed[n]) {
				++can;
			}
		}
		return std::make_pair(can, needs);
	}
};
