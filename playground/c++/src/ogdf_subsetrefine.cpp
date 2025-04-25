#include "ogdf_subsetrefine.hpp"
#include <queue>

// void bfsorder(ogdf::Graph G, std::list<ogdf::node>& order) {
//     ogdf::NodeArray<bool> mark(G, false);
//     std::queue<ogdf::node> bfs;
//     if (G.numberOfNodes() == 0) {
//         return;
//     }
//     ogdf::node s = G.chooseNode();
//     bfs.push(s);
//     // mark s and set distance to itself 0
//     mark[s] = true;
//     order.push_back(s);
//     while (!bfs.empty()) {
//         log <<bfs.size() <<std::endl;
//         ogdf::node w = bfs.front();
//         bfs.pop();
//         for (ogdf::adjEntry adj : w->adjEntries) {
//             ogdf::node v = adj->twinNode();
//             if (!mark[v]) {
//                 mark[v] = true;
//                 bfs.push(v);
//                 order.push_back(v);
//             }
//         }
//     }
// }

size_t SubsetRefine::doRefinementReduction() {
    std::vector<ogdf::node> order;

    for (auto u : G.nodes) {
        if ((!instance.is_subsumed(u) && type == RefineType::Dominate) || (!instance.is_dominated(u) && type == RefineType::Subsume)) {
            order.push_back(u);
        }
    }
    std::sort(order.begin(), order.end(), [&](ogdf::node a, ogdf::node b) {
        if (type == RefineType::Dominate) {
            return instance.countCanDominate(a) < instance.countCanDominate(b); 
        }
        else {
            return instance.countCanBeDominatedBy(a) < instance.countCanBeDominatedBy(b);
        }
    });
    // std::list<ogdf::node> order;
    // bfsorder(G, order);
    for (auto u : order) {
        if ((!instance.is_subsumed(u) && type == RefineType::Dominate) || (!instance.is_dominated(u) && type == RefineType::Subsume)) {
            refineByNode(u);
        }
    }
    log << "Partition refinement structure has " << refineG.numberOfNodes() << " nodes and " << refineG.numberOfEdges() << " edges"<<std::endl;
    log << "In total, the structure had " << cntedgesadded << " edges added."<<std::endl;
    for (auto bag: refineG.nodes) {
        if (bag->indeg() > 0) {
            auto vec = bagNodeVec[bag];
            for (auto v: vec) {
                doReduce(v);
            }
        }
        else if (bagNodeVec[bag].size() > 1) {            
            // reduce all but one
            auto vec = bagNodeVec[bag];
            bool someonestay = false;
            for (size_t i = 0; i < vec.size(); i++) {
                if (i < vec.size() - 1 || someonestay) {
                   someonestay |= !doReduce(vec[i]);
                }
            }
        }
    }
    return cntreduced;
}

void SubsetRefine::refineByNode(const ogdf::node& u) {
    std::vector<ogdf::node> nodesToTouch;

    if (type == RefineType::Subsume) {
        instance.forAllCanBeDominatedBy(u, [&](ogdf::node adj) {
            if (bagof[adj] == nullptr) return true;
            nodesToTouch.push_back(adj);
            return true;
        });
    }
    else {
        instance.forAllCanDominate(u, [&](ogdf::node adj) {
            if (bagof[adj] == nullptr) return true;
            nodesToTouch.push_back(adj);
            return true;
        });
    }

    // move vertices into new bags
    std::vector<ogdf::node> touchedbags;
    for (auto v: nodesToTouch) {
        cnttouched[v]++;
        auto bag = bagof[v];
        ogdf::node newbag;
        if (refinedBag[bag] == nullptr) 
        {
            newbag = refineG.newNode();
            touchedbags.push_back(bag);  
            refinedBag[bag] = newbag;          
        }
        else {
            newbag = refinedBag[bag];
        }
        // swap and pop from old bag
        auto& vec = bagNodeVec[bag];
        auto oldindex = vecIndex[v];
        vecIndex[vec[vec.size()-1]] = oldindex;
        std::swap(vec[oldindex], vec[vec.size()-1]);
        vec.pop_back();
        // insert into new bag
        auto& newvec = bagNodeVec[newbag];
        bagof[v] = newbag;
        vecIndex[v] = newvec.size();
        newvec.push_back(v);
    }

    // update bag adjacencies
    for (auto bag: touchedbags) {
        auto newbag = refinedBag[bag];
        if (type == RefineType::Subsume) {
            // the refined bag has vertices with superset outedges
            refineG.newEdge(newbag, ogdf::Direction::before, bag, ogdf::Direction::after);
            // everything that was subsumed by old bag is also subsumed by new bag, but also the refined bags
            cntedgesadded++;

            forAllOutAdj(bag, [&](ogdf::adjEntry adj) {
                refineG.newEdge(newbag, ogdf::Direction::before, adj->twinNode(), ogdf::Direction::after);
                cntedgesadded++;
                if (refinedBag[adj->twinNode()] != nullptr) {
                    refineG.newEdge(newbag, ogdf::Direction::before, refinedBag[adj->twinNode()], ogdf::Direction::after);
                    cntedgesadded++;
                }
                return true;
            });
        }
        else {
            // the refined bag has vertices with superset inedges
            refineG.newEdge(bag, ogdf::Direction::before, newbag, ogdf::Direction::after);
            cntedgesadded++;
            forAllInAdj(bag, [&](ogdf::adjEntry adj) {
                refineG.newEdge(adj->twinNode(), ogdf::Direction::before, newbag, ogdf::Direction::after);
                if (refinedBag[adj->twinNode()] != nullptr) {
                    refineG.newEdge(refinedBag[adj->twinNode()], ogdf::Direction::before, newbag, ogdf::Direction::after);
                }
                return true;
            });
        }
    }

    // remove bags without nodes in it and reset refinedBag
    for (auto it = touchedbags.begin(); it != touchedbags.end(); ) {
        auto bag = *it;
        
        ++it;
        if (bagNodeVec[bag].size() == 0) {
            refineG.delNode(bag);
        }
        refinedBag[bag] = nullptr;
    }

    for (auto v: nodesToTouch) {
        size_t deg = type == RefineType::Subsume ? instance.countCanDominate(v) : instance.countCanBeDominatedBy(v);

        if (cnttouched[v] == deg) {
            if (type == RefineType::Subsume && (bagof[v]->indeg() > 0 || bagNodeVec[bagof[v]].size() > 1)) {
                // v will be subsumed in the end either way by incoming edges or by someone in the bag
                auto& vec = bagNodeVec[bagof[v]];
                auto bag = bagof[v];                
                doReduce(v);

                if (vec.size() == 0) {
                    refineG.delNode(bag);
                }
            }
            else if (type == RefineType::Dominate) {
                // the subset relations for outedges will not change anymore so we can already dominate all outedges and nodes in the bag
                forAllOutAdj(bagof[v], [&](ogdf::adjEntry adj) {
                    for (auto w: std::vector<ogdf::node>(bagNodeVec[adj->twinNode()].begin(), bagNodeVec[adj->twinNode()].end())) {
                        doReduce(w);
                    }
                    if (bagNodeVec[adj->twinNode()].size() == 0) {
                        refineG.delNode(adj->twinNode());
                    }
                    return true;
                });
                for (auto w: std::vector<ogdf::node>(bagNodeVec[bagof[v]].begin(), bagNodeVec[bagof[v]].end())) {
                    if (w != v) {
                        doReduce(w);
                    }
                }
                if (bagNodeVec[bagof[v]].size() == 0) {
                    refineG.delNode(bagof[v]);
                }
            }
        }
    }
}