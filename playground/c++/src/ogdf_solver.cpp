#include "ogdf.hpp"

#include "EvalMaxSAT.h"

void solveGreedy(Instance &I) {
    log << "Solving Greedy with " << I.G.numberOfNodes() << " nodes" << std::endl;
    log << "Old DS size: " << I.DS.size() << std::endl;
    while (!I.G.empty()) {
        ogdf::node bestVertex = nullptr;
        int bestValue = -1;
        for (auto v : I.G.nodes) {
            int wouldbecovered = v->outdeg();
            if (I.is_dominated[v])
                wouldbecovered++;
            if (wouldbecovered > bestValue) {
                bestValue = wouldbecovered;
                bestVertex = v;
            }
        }
        if (bestVertex == nullptr) {
            break;
        }
        I.addToDominatingSet(bestVertex);
    }
    log << "New DS size: " << I.DS.size() << std::endl;
}

void solveEvalMaxSat(Instance &I) {
    EvalMaxSAT solver;
    const int n = I.G.numberOfNodes();
    log << "Solving ILP with " << n << " nodes" << std::endl;
    ogdf::NodeArray<int> varmap(I.G, -1);
    for (auto v : I.G.nodes) {
        if (I.is_subsumed[v]) {
            continue;
        }
        auto var = solver.newVar();
        varmap[v] = var;
        solver.addClause({-var}, 1); // soft clause
    }
    for (auto v : I.G.nodes) {
        if (I.is_dominated[v]) {
            continue;
        }
        std::vector<int> clause;
        if (!I.is_subsumed[v]) {
            clause.push_back(varmap[v]);
        }
        forAllInAdj(v,
                    [&](ogdf::adjEntry adj) {
                        auto w = adj->twinNode();
                        if (!I.is_subsumed[w]) {
                            clause.push_back(varmap[w]);
                        }
                        return true;
                    });
        solver.addClause(clause); //hard clause
    }
    solver.setTargetComputationTime(10 * 60);
    std::cout.setstate(std::ios::failbit); // https://stackoverflow.com/a/8246430
    bool solved = solver.solve();
    std::cout.clear();
    log << "Solving result: " << solved << std::endl;
    log << "Old DS size: " << I.DS.size() << std::endl;
    for (auto v : I.G.nodes) {
        if (!I.is_subsumed[v]) {
            if (solver.getValue(varmap[v])) {
                I.DS.push_back(I.node2ID[v]);
            }
        }
    }
    log << "New DS size: " << I.DS.size() << std::endl;
}
