#pragma once
#include <queue>
#include <stdexcept>
#include <string>
#include <vector>
#include <cstdint>
#include <limits>

// from https://rosettacode.org/wiki/Hopcroft-Karp_Algorithm

/**
 * Representation of a bipartite graph.
 * Vertices in the left partition, U, are numbered from 1 to m,
 * and vertices in the right partition, V, are numbered 1 to n.
 */
class BipartiteGraph {
public:
	BipartiteGraph(const uint32_t& aM, const uint32_t& aN) {
		m = aM;
		n = aN;

		adjacency_lists = { m + 1, std::vector<uint32_t>() };
		pair_u.assign(m + 1, NIL);
		pair_v.assign(n + 1, NIL);
		levels.assign(m + 1, INF);
	}

	void add_edge(const uint32_t& u, const uint32_t& v) {
		if ( 1 <= u && u <= m && 1 <= v && v <= n ) {
			adjacency_lists[u].emplace_back(v);
		} else {
			throw std::invalid_argument("Attempt to add an edge (" +
				std::to_string(u) + ", " + std::to_string(v) + ") which is out of bounds");
		}
	}

	/**
	 * Return the matching size of the bipartite graph.
	 */
	uint32_t hopcroftKarp_algorithm() {
		pair_u.assign(m + 1, NIL);
		pair_v.assign(n + 1, NIL);
		uint32_t matching_size = 0;

		while ( breadth_first_search() ) {
			for ( uint32_t u = 1; u <= m; ++u ) {
				if ( pair_u[u] == NIL && depth_first_search(u) ) { // vertex u is free and an augmenting path starting
					matching_size++;                               // from u has been found by the depth first search
				}
			}
		}
		return matching_size;
	}

private:
	/**
	 * Determines whether there exists an augmenting path starting from a free vertex in U.
	 *
	 * Return true if an augmenting path could exist, otherwise false.
	 */
	bool breadth_first_search() {
		std::queue<uint32_t> queue;
		for ( uint32_t u = 1; u <= m; ++u ) { // Initialise 'levels' for the vertices in U
			if ( pair_u[u] == NIL ) { // If u is a free vertex, its level is 0 add it is added to the queue
				levels[u] = 0;
				queue.push(u);
			} else { // Otherwise, set 'levels' to infinity
				levels[u] = INF;
			}
		}

		// The 'level' to the NIL node represents the length of the shortest augmenting path
		levels[NIL] = INF;

		while ( ! queue.empty() ) {
			const uint32_t u = queue.front();
			queue.pop();
			if ( levels[u] < levels[NIL] ) { // The path through u could lead to a shorter augmenting path
				for ( const uint32_t& v : adjacency_lists[u] ) { // Explore the neighbours v of u in V
					const uint32_t matched_u = pair_v[v];
					if ( levels[matched_u] == INF ) { // The matched vertex has not been visited yet
					   levels[matched_u] = levels[u] + 1;
						queue.push(matched_u); // Enqueue the matched vertex to explore it further
					}
				}
			}
		}

		// An augmenting path from the initial free vertices was found if levels[NIL] is not INF
		return levels[NIL] != INF;
	}

	/**
	 * Determine whether the shortest path from vertex u in U found by breadth_first_search() can be augmented.
	 *
	 * Return true if an augmenting path was found starting from u, otherwise false.
	 */
	bool depth_first_search(const uint32_t& u) {
		if ( u != NIL ) {
			for ( const uint32_t& v : adjacency_lists[u] ) { // Explore neighbours v of u in V
				const uint32_t matched_u = pair_v[v];
				// Check whether the edge (u, v) leads to a vertex matched_u
				// such that the path u -> v -> matched_u is part of a shortest augmenting path
				if ( levels[matched_u] == levels[u] + 1 ) {
					if ( depth_first_search(matched_u) ) { // An augmenting path is found starting from 'matched_u'
						pair_v[v] = u; // Match v with u,
						pair_u[u] = v; // and u with v
						return true;
					}
				}
			}

			// No augmenting path was found starting from vertex u through any of its neighbours v,
			// so remove u from the depth first search phase of the algorithm
			levels[u] = INF;
			return false;
		}

		return true;
	}

	std::vector<std::vector<uint32_t>> adjacency_lists; // adjacency_lists(u) stores a list of neighbours of u in V
	std::vector<uint32_t> pair_u; // pair_u(u) stores the vertex v in V matched with u in U, or NIL if unmatched
	std::vector<uint32_t> pair_v; // pair_v(v) stores the vertex u in U matched with v in V, or NIL if unmatched
	std::vector<uint32_t> levels; // levels(u) stores the level of vertex u in U during a breadth first search

	uint32_t m; // Index of the vertices in the left partition
	uint32_t n; // Index of the vertices in the right partition

	const uint32_t NIL = 0;
	const uint32_t INF = std::numeric_limits<u_int32_t>::max();
};

struct Edge {
	uint32_t from;
	uint32_t to;
};