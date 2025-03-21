#pragma once

#include "ogdf/basic/Graph.h"
#include <chrono>
#include <filesystem>

extern ogdf::Logger logger;

#ifdef OGDF_DEBUG
#define log if (true) logger.lout()
#else
#define log if (false) logger.lout()
#endif

#ifdef OGDF_DEBUG
inline void forAllOutAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
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
        if (call && !f(adj)) call = false;
        ++c;
    }
    while (c < degree) {
        OGDF_ASSERT(adj_it != v->adjEntries.end());
        OGDF_ASSERT(!(*adj_it)->isSource());
        ++adj_it;
        ++c;
    }
    OGDF_ASSERT(adj_it == v->adjEntries.end());
}

inline void forAllInAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
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
        if (call && !f(adj)) call = false;
        ++c;
    }
    while (c < degree) {
        OGDF_ASSERT(adj_it != v->adjEntries.rend());
        OGDF_ASSERT((*adj_it)->isSource());
        ++adj_it;
        ++c;
    }
    OGDF_ASSERT(adj_it == v->adjEntries.rend());
}
#else
inline void forAllOutAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
    for (auto adj_it = (v)->adjEntries.begin(); adj_it != (v)->adjEntries.end();) {
        auto adj = *adj_it;
        ++adj_it;
        if (!adj->isSource() || !f(adj)) break;
    }
}
inline void forAllInAdj(ogdf::node v, std::function<bool(ogdf::adjEntry)> f) {
    for (auto adj_it = (v)->adjEntries.rbegin(); adj_it != (v)->adjEntries.rend();) {
        auto adj = *adj_it;
        ++adj_it;
        if (adj->isSource() || !f(adj)) break;
    }
}
#endif

namespace internal {
ogdf::node idn(ogdf::node n);
ogdf::edge ide(ogdf::edge n);
}

#define SMALL_BLOCK 100