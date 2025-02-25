#include "preprocessing.hpp"


bool reductionIsolated(Instance *instance, VertexList& dominatingSet)
{
    Graph::vertex_iterator v, vend;
    bool reduced = false;
    for (boost::tie(v, vend) = IteratorVertices((*instance->G)); v != vend; )
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
    for (boost::tie(v, vend) = IteratorVertices((*instance->G)); v != vend; )
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
    bool progress = true;
    bool reduced = false;

    while (progress)
    {
        progress = false;
        unordered_set<VD> toDelete;

        for (auto v : MyVertices((*instance->G)))
        {
            if (boost::out_degree(v, (*instance->G)) == 1)
            {
                progress = true;
                reduced = true;
                auto w = *boost::adjacent_vertices(v, (*instance->G)).first;   
                if ((*instance->G)[v].is_dominated)
                {
                    // toDelete.insert(v); THIS IS WRONG
                }
                else if ((*instance->G)[w].can_be_dominating_set)
                {                                                         
                    dominatingSet.push_back((*instance->G)[w].id);
                    instance->toDominatingSet(w, toDelete);
                }
            }
        }
        instance->deleteVertices(toDelete);
    }

    // reductionCovered(instance); // for sanity. TODO: can probably be removed
    return reduced;
}

list<Instance*> decomposeConnectedComponents(Instance *instance)
{
    boost::property_map<Graph, boost::vertex_index_t>::type index_map = boost::get(boost::vertex_index, (*instance->G));
    int cnt = 0;
    for (auto v: MyVertices((*instance->G))) {
        boost::put(index_map, v, cnt++);
    }
    int number_connected_components = boost::connected_components((*instance->G), index_map);
    vector<Graph*> subgraphs(number_connected_components);
    for (auto &g: subgraphs) {
        g = new Graph();
    }
    for (VD v: MyVertices((*instance->G))) {        
        auto idx = boost::get(index_map, v);
        VD new_v = boost::add_vertex(*subgraphs[idx]);
        (*subgraphs[idx])[new_v] = (*instance->G)[v];
        instance->props->id_to_vertex[(*instance->G)[v].id] = new_v;
        (*instance->G)[v].forward = new_v;
    }
    for (auto e: Edges((*instance->G))) {
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
    
    for (auto v: MyVertices((*instance->G))) {
        for (auto w: Neighbors((*instance->G), v)) {
            is_neighbor[(*instance->G)[w].id] = true;
        }        
        for (auto w: Neighbors((*instance->G), v)) {
            bool dominated = true;
            for (auto u: Neighbors((*instance->G), w)) {
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
        for (auto w: Neighbors((*instance->G), v)) {
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
    while (partitionElement != NULL) {
        // delete all but one twin
        for (size_t i = 1; i < partitionElement->elements.size(); i++) {
            instance->deleteVertex(instance->props->id_to_vertex[partitionElement->elements[i]]);
            reduced = true;
        }
        partitionElement = partitionElement->next;
    }
    delete partition;
    return reduced;
}

bool reductionByCanBeDominatingSet(Instance* instance, VertexList& dominatingSet) {
    unordered_set<VD> toDelete;
    list<VD> toDominatingSet;
    bool reduced = false;
    for (VD v: MyVertices((*instance->G))) {
        auto vi = (*instance->G)[v];
        if (vi.can_be_dominating_set == false && vi.is_dominated == true) {
            reduced = true;
            toDelete.insert(v);
        }
        else if (vi.can_be_dominating_set && vi.is_dominated == false) {
            bool vshouldbedominatingset = true;
            for (VD w: Neighbors((*instance->G), v)) {
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
    for (VD v: MyVertices((*instance->G)) ) {
        if (!(*instance->G)[v].is_dominated) {
            someoneneedsdomination = true;
            break;
        }
    }
    if (someoneneedsdomination)
    {
        for (VD v: MyVertices((*instance->G)) ) {
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
    for (boost::tie(vi, vend) = IteratorVertices((*instance->G)); vi != vend; )
    {
        auto v = *vi;
        if ((*instance->G)[v].can_be_dominating_set == false) {
            vi++;
            continue;
        }
        visited[(*instance->G)[v].id] = 1;
        bool someonedominated = false;
        for (auto w: Neighbors((*instance->G), v)) {
            visited[(*instance->G)[w].id] = 1;
            if (!(*instance->G)[w].can_be_dominating_set) {
                someonedominated = true;
            }
        }
        if (!someonedominated) {
            visited[(*instance->G)[v].id] = 0;
            for (auto w: Neighbors((*instance->G), v)) {
                visited[(*instance->G)[w].id] = 0;
            }
            vi++;
            continue;
        }
        list<VD> N1v;
        for (auto w: Neighbors((*instance->G), v)) {
            for (auto u: Neighbors((*instance->G), w)) {
                if (visited[(*instance->G)[u].id] == 0) {
                    N1v.push_back(w);
                    visited[(*instance->G)[w].id] = 2;
                    break;
                }
            }
        }

        list<VD> N2v;
        for (auto w: N1v) {
            for (auto u: Neighbors((*instance->G), w)) {
                if (visited[(*instance->G)[u].id] == 1 && u != v && !N2[(*instance->G)[u].id]) {
                    N2[(*instance->G)[u].id] = true;
                    N2v.push_back(u);
                    visited[(*instance->G)[u].id] = 2;
                }
            }
        }

        list<VD> N3v;
        for (auto w: Neighbors((*instance->G), v)) {
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

        for (auto w: Neighbors((*instance->G), v)) {
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