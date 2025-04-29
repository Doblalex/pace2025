#pragma once

#include <chrono>
#include <filesystem>

#include "ogdf/basic/Graph.h"

extern ogdf::Logger logger;

#if defined(OGDF_DEBUG) || defined(PACE_LOG)
#	define log   \
		if (true) \
		logger.lout()
#	define logd  \
		if (true) \
		logger.lout(ogdf::Logger::Level::Minor)
#else
#	define log    \
		if (false) \
		logger.lout()
#	define logd   \
		if (false) \
		logger.lout(ogdf::Logger::Level::Minor)
#endif

#ifdef OGDF_DEBUG
inline bool forAllOutAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
	OGDF_ASSERT(v->outdeg() == 0 || v->adjEntries.head()->isSource());
	OGDF_ASSERT(v->indeg() == 0 || !v->adjEntries.tail()->isSource());
	size_t c = 0;
	int outdeg = v->outdeg();
	int degree = v->degree();
	bool call = true;
	auto adj_it = v->adjEntries.begin();
	while (c < outdeg) {
		OGDF_ASSERT(adj_it != v->adjEntries.end());
		auto adj = *adj_it;
		OGDF_ASSERT(adj->isSource());
		++adj_it;
		if (call && !f(adj)) {
			call = false;
		}
		++c;
	}
	while (c < degree) {
		OGDF_ASSERT(adj_it != v->adjEntries.end());
		OGDF_ASSERT(!(*adj_it)->isSource());
		++adj_it;
		++c;
	}
	// FIXME reenable
	// OGDF_ASSERT(adj_it == v->adjEntries.end());
	return call;
}

inline bool forAllInAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
	OGDF_ASSERT(v->outdeg() == 0 || v->adjEntries.head()->isSource());
	OGDF_ASSERT(v->indeg() == 0 || !v->adjEntries.tail()->isSource());
	size_t c = 0;
	int indeg = v->indeg();
	int degree = v->degree();
	bool call = true;
	auto adj_it = v->adjEntries.rbegin();
	while (c < indeg) {
		OGDF_ASSERT(adj_it != v->adjEntries.rend());
		auto adj = *adj_it;
		OGDF_ASSERT(!adj->isSource());
		++adj_it;
		if (call && !f(adj)) {
			call = false;
		}
		++c;
	}
	while (c < degree) {
		OGDF_ASSERT(adj_it != v->adjEntries.rend());
		OGDF_ASSERT((*adj_it)->isSource());
		++adj_it;
		++c;
	}
	OGDF_ASSERT(adj_it == v->adjEntries.rend());
	return call;
}
#else
inline bool forAllOutAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
	for (auto adj_it = (v)->adjEntries.begin(); adj_it != (v)->adjEntries.end();) {
		auto adj = *adj_it;
		++adj_it;
		if (!adj->isSource()) {
			return true;
		}
		if (!f(adj)) {
			return false;
		}
	}
	return true;
}

inline bool forAllInAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
	for (auto adj_it = (v)->adjEntries.rbegin(); adj_it != (v)->adjEntries.rend();) {
		auto adj = *adj_it;
		++adj_it;
		if (adj->isSource()) {
			return true;
		}
		if (!f(adj)) {
			return false;
		}
	}
	return true;
}
#endif

namespace internal {
ogdf::node idn(ogdf::node n);
ogdf::edge ide(ogdf::edge n);
}

#define SMALL_BLOCK 100
#define BLOCK_FRACTION 0.25f

constexpr uint64_t FNV1a_64_SEED = 0xcbf29ce484222325UL;

// https://stackoverflow.com/a/77342581
inline void FNV1a_64_update(uint64_t& h, uint64_t v) {
	h ^= v;
	h *= 0x00000100000001B3UL;
}

inline uint64_t FNV1a_64_one(uint64_t v) {
	uint64_t ret = FNV1a_64_SEED;
	FNV1a_64_update(ret, v);
	return ret;
}

#define EMS_CACHE
// #define USE_ORTOOLS
