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
            dominatingSet.push_back((*instance->G)[vd].id);
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

        for (auto v : Vertices((*instance->G)))
        {
            if (boost::out_degree(v, (*instance->G)) == 1)
            {
                progress = true;
                reduced = true;
                if ((*instance->G)[v].is_dominated)
                {
                    toDelete.insert(v);
                }
                else
                {
                    auto w = *boost::adjacent_vertices(v, (*instance->G)).first;
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
    for (auto v: Vertices((*instance->G))) {
        boost::put(index_map, v, cnt++);
    }
    int number_connected_components = boost::connected_components((*instance->G), index_map);
    vector<Graph*> subgraphs(number_connected_components);
    for (auto &g: subgraphs) {
        g = new Graph();
    }
    for (VD v: Vertices((*instance->G))) {        
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
    
    for (auto v: Vertices((*instance->G))) {
        for (auto w: Neighbors((*instance->G), v)) {
            is_neighbor[(*instance->G)[w].id] = true;
        }        
        for (auto w: Neighbors((*instance->G), v)) {
            bool dominated = true;
            for (auto u: Neighbors((*instance->G), w)) {
                if (u != v && !is_neighbor[(*instance->G)[u].id]) {
                    dominated = false;
                    break;
                }
            }
            if (dominated && (*instance->G)[v].can_be_dominating_set == true) {
                (*instance->G)[v].can_be_dominating_set = false;
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
    for (VD v: Vertices((*instance->G))) {
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