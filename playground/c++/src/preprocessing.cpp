#include "preprocessing.hpp"


bool reductionIsolated(Instance *instance, VertexList& dominatingSet)
{
    Graph::vertex_iterator v, vend;
    bool reduced = false;
    for (boost::tie(v, vend) = boost::vertices((*instance->G)); v != vend; )
    {   
        auto vd = *v;
        ++v;
        if (boost::out_degree(vd, (*instance->G)) == 0)
        {
            reduced = true;
            if (!(*instance->G)[vd].is_dominated)
            {
                dominatingSet.push_back((*instance->G)[vd].id);
            }
            instance->deleteVertex(vd);
        }
    }
    return reduced;
}

bool reductionCovered(Instance *instance)
{
    Graph::vertex_iterator v, vend;
    bool reduced = false;
    for (boost::tie(v, vend) = boost::vertices((*instance->G)); v != vend; )
    {
        auto vd = *v;
        ++v;
        if ((*instance->G)[vd].is_dominated && (boost::out_degree(vd, (*instance->G)) <= 1 || (*instance->G)[vd].cnt_undominated_neighbors == 0))
        {
            instance->deleteVertex(vd);
            reduced = true;
        }
    }
    return reduced;
}

bool reductionDegree1(Instance *instance, VertexList& dominatingSet)
{
    // TODO: this can be done faster by keeping track of degree 1 vertices
    bool reduced = false;
    Graph::vertex_iterator vi, vend;
    for (boost::tie(vi, vend) = boost::vertices((*instance->G)); vi != vend; )
    {
        auto v = *vi;
        if (boost::out_degree(v, (*instance->G)) == 1)
        {            
            auto w = *boost::adjacent_vertices(v, (*instance->G)).first;   
            if ((*instance->G)[v].is_dominated && ((*instance->G)[w].can_be_dominating_set || (*instance->G)[w].is_dominated))
            {
                vi++;                
                instance->deleteVertex(v);
                reduced = true;
            }
            else if ((*instance->G)[w].can_be_dominating_set)
            {              
                reduced = true;                               
                unordered_set<VD> toDelete;        
                dominatingSet.push_back((*instance->G)[w].id);
                instance->toDominatingSet(w, toDelete);
                instance->deleteVerticesWithIterator(toDelete, vi);
            }
        }
        vi++;
    }

    // reductionCovered(instance); // for sanity. TODO: can probably be removed
    return reduced;
}

list<Instance*> decomposeConnectedComponents(Instance *instance)
{
    boost::property_map<Graph, boost::vertex_index_t>::type index_map = boost::get(boost::vertex_index, (*instance->G));
    int cnt = 0;
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        boost::put(index_map, v, cnt++);
    }
    int number_connected_components = boost::connected_components((*instance->G), index_map);
    vector<Graph*> subgraphs(number_connected_components);
    for (auto &g: subgraphs) {
        g = new Graph();
    }
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {   
        auto idx = boost::get(index_map, v);
        VD new_v = boost::add_vertex(*subgraphs[idx]);
        (*subgraphs[idx])[new_v] = (*instance->G)[v];
        instance->props->id_to_vertex[(*instance->G)[v].id] = new_v;
        (*instance->G)[v].forward = new_v;
    }
    
    BGL_FORALL_EDGES(e, *instance->G, Graph) {
        auto u = boost::source(e, (*instance->G));
        auto v = boost::target(e, (*instance->G));
        auto idx_u = boost::get(index_map, u);
        boost::add_edge((*instance->G)[u].forward, (*instance->G)[v].forward, *subgraphs[idx_u]);
    }
    list<Instance*> instances;
    for (auto &g: subgraphs) {
        Instance* inst = new Instance(instance->props);
        instances.emplace_back(inst);
        inst->G = g;
        inst->n = boost::num_vertices(*g);
        inst->m = boost::num_edges(*g);
    }
    delete instance;
    return instances;
}

bool reductionDomination(Instance* instance) {
    bool reduced = false;
    static vector<bool> is_neighbor(instance->props->n+1, false);
    
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            is_neighbor[(*instance->G)[w].id] = true;
        }        
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            bool dominated = true;
            BGL_FORALL_ADJ_T(w, u, *instance->G, Graph) {
                if (u != v && !is_neighbor[(*instance->G)[u].id] && !(*instance->G)[u].is_dominated) {
                    dominated = false;
                    break;
                }
            }
            if (dominated && (*instance->G)[w].can_be_dominating_set == true && (*instance->G)[v].can_be_dominating_set == true) {
                (*instance->G)[w].can_be_dominating_set = false;
                reduced = true;
            }
        }
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            is_neighbor[(*instance->G)[w].id] = false;
        }
    }
    return reduced;
}

bool reductionTwins(Instance* instance) {
    PartitionRefinement* partition = new PartitionRefinement(instance->G, instance->props);
    partition->dorefine();
    auto partitionElement = partition->start;
    bool reduced = false;
    unordered_set<VD> toDelete;
    while (partitionElement != NULL) {
        if (partitionElement->elements.size() > 1) {
            size_t keep = -1;
            for (size_t i = 0; i < partitionElement->elements.size(); i++) {
                VD v1 = instance->props->id_to_vertex[partitionElement->elements[i]];
                if ((*instance->G)[v1].can_be_dominating_set && !(*instance->G)[v1].is_dominated) { // TODO: keep this safe for now, maybe can be improved. But without this we get wrong results
                    keep = i;
                    break;
                }
            }
            if (keep != -1) {
                for (size_t i = 0; i < partitionElement->elements.size(); i++) {
                    if (i != keep) {
                        VD v2 = instance->props->id_to_vertex[partitionElement->elements[0]];   
                        if (!(*instance->G)[v2].is_dominated) {
                            BGL_FORALL_ADJ_T(v2, w, *instance->G, Graph) {
                                (*instance->G)[w].cnt_undominated_neighbors--;
                            }
                        }       
                        toDelete.insert(v2);
                        reduced = true;
                    }
                }
                
            }
        }
        
        partitionElement = partitionElement->next;
    }
    instance->deleteVertices(toDelete);
    delete partition;
    return reduced;
}

bool reductionByCanBeDominatingSet(Instance* instance, VertexList& dominatingSet) {
    unordered_set<VD> toDelete;
    list<VD> toDominatingSet;
    bool reduced = false;
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        auto vi = (*instance->G)[v];
        if (vi.can_be_dominating_set == false && vi.is_dominated == true) {
            reduced = true;
            toDelete.insert(v);
        }
        else if (vi.can_be_dominating_set && vi.is_dominated == false) {
            bool vshouldbedominatingset = true;
            BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
                auto viw = (*instance->G)[w];
                if (viw.can_be_dominating_set) {
                    vshouldbedominatingset = false;
                    break;
                }
            }
            if (vshouldbedominatingset) {
                reduced = true;
                dominatingSet.push_back(vi.id);
                toDominatingSet.push_back(v);
            }
        }
    }
    for (auto v: toDominatingSet) {
        instance->toDominatingSet(v, toDelete);
    }
    instance->deleteVertices(toDelete);

    return reduced;
    
}

bool reduceUniversal(Instance* instance, VertexList& dominatingSet) {
    bool reduced = false;
    bool someoneneedsdomination = false;
    BGL_FORALL_VERTICES(v, *instance->G, Graph) {
        if (!(*instance->G)[v].is_dominated) {
            someoneneedsdomination = true;
            break;
        }
    }
    if (someoneneedsdomination)
    {
        BGL_FORALL_VERTICES(v, *instance->G, Graph) {
            if (boost::out_degree(v, (*instance->G)) == instance->n-1) {
                dominatingSet.push_back((*instance->G)[v].id);
                unordered_set<VD> toDelete;
                instance->toDominatingSet(v, toDelete);
                instance->deleteVertices(toDelete);
                return true;
            }
        }
    }
    return false;
}

bool reductionDominationPaper(Instance* instance, VertexList& dominatingSet) {
    static vector<int> visited(instance->props->n+1, 0);
    static vector<bool> N2(instance->props->n+1, false);
    
    bool reduced = false;
    Graph::vertex_iterator vi, vend;
    for (boost::tie(vi, vend) = boost::vertices((*instance->G)); vi != vend; )
    {
        auto v = *vi;
        if ((*instance->G)[v].can_be_dominating_set == false) {
            vi++;
            continue;
        }
        visited[(*instance->G)[v].id] = 1;
        bool someonedominated = false;
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            visited[(*instance->G)[w].id] = 1;
            if (!(*instance->G)[w].can_be_dominating_set) {
                someonedominated = true;
            }
        }
        if (!someonedominated) {
            visited[(*instance->G)[v].id] = 0;
            BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
                visited[(*instance->G)[w].id] = 0;
            }
            vi++;
            continue;
        }
        list<VD> N1v;
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            BGL_FORALL_ADJ_T(w, u, *instance->G, Graph) {
                if (visited[(*instance->G)[u].id] == 0) {
                    N1v.push_back(w);
                    visited[(*instance->G)[w].id] = 2;
                    break;
                }
            }
        }

        list<VD> N2v;
        for (auto w: N1v) {
            BGL_FORALL_ADJ_T(w, u, *instance->G, Graph) {
                if (visited[(*instance->G)[u].id] == 1 && u != v && !N2[(*instance->G)[u].id]) {
                    N2[(*instance->G)[u].id] = true;
                    N2v.push_back(u);
                    visited[(*instance->G)[u].id] = 2;
                }
            }
        }

        list<VD> N3v;
        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            if (visited[(*instance->G)[w].id] == 1) {
                N3v.push_back(w);
            }           
        }      
        
        bool someoneneedsdomination = false;
        for (auto w: N3v) {
            if (!(*instance->G)[w].is_dominated) {
                someoneneedsdomination = true;
                break;
            }
        }

        BGL_FORALL_ADJ_T(v, w, *instance->G, Graph) {
            visited[(*instance->G)[w].id] = 0;
            N2[(*instance->G)[w].id] = false;
        }
        visited[(*instance->G)[v].id] = 0;

        if (someoneneedsdomination) {
            // debug("Reduction paper");
            reduced = true;
            dominatingSet.push_back((*instance->G)[v].id);
            unordered_set<VD> toDelete;
            instance->toDominatingSet(v, toDelete);
            if (instance->deleteVerticesWithIterator(toDelete, vi)) {
                continue;
            }
        }
        vi++;
    }
    
    return reduced;
}