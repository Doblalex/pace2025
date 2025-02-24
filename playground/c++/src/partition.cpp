#include "partition.hpp"
vector<PartitionElement*> partitionelementAt;
vector<size_t> vectorElementAt; 
// for pop and swap

PartitionRefinement::PartitionRefinement(Graph* graph, globalprops* props) : graph(graph), props(props) {
    partitionelementAt.resize(props->n+1);
    vectorElementAt.resize(props->n+1);
    this->start = new PartitionElement();
    this->end = this->start;
    size_t i = 0;
    for (auto v: Vertices(*graph)) {        
        partitionelementAt[(*graph)[v].id] = this->start;            
        this->start->elements.push_back((*graph)[v].id);
        vectorElementAt[(*graph)[v].id] = i++;
    }
}

bool PartitionRefinement::dorefine()
{
    for (VD v: Vertices(*graph)) {
        vector<Vertex> ids;
        ids.emplace_back((*graph)[v].id);
        for (VD neighbor: Neighbors(*graph, v)) {
            ids.emplace_back((*graph)[neighbor].id);
        }
        refine(ids);
        bool found = !(find(start->elements.begin(), start->elements.end(), 7) == start->elements.end());
    }
    return true;
}

void PartitionRefinement::refine(vector<Vertex> arr)
{
    vector<PartitionElement *> refinedEls;
    for (auto x : arr)
    {        
        PartitionElement *el = partitionelementAt[x];
        if (el->refined.size() == 0)
            refinedEls.push_back(el);
        el->refined.push_back(x);
        if (el->elements[vectorElementAt[x]] != x) exit(0);
        swap(el->elements[vectorElementAt[x]], el->elements[el->elements.size() - 1]);
        vectorElementAt[el->elements[vectorElementAt[x]]] = vectorElementAt[x];               
        el->elements.pop_back();

    }
    for (auto el : refinedEls)
    {
        PartitionElement *newel = new PartitionElement();
        size_t i = 0;
        for (auto x : el->refined)
        {
            newel->elements.push_back(x);
            vectorElementAt[x] = i++;
            partitionelementAt[x] = newel;
        }
        el->refined.clear();
        newel->last = this->end;
        newel->last->next = newel;
        this->end = newel;         

        if (el->elements.size() == 0)
        {
            if (el->next == NULL)
            {
                el->last->next = NULL;
                this->end = el->last;
            }
            else if (el->last == NULL)
            {
                el->next->last = NULL;
                this->start = el->next;
            }
            else
            {
                el->last->next = el->next;
                el->next->last = el->last;
            }
            delete el;
        }
    }
}

void PartitionRefinement::print(int minsize)
{
    auto el = this->start;
    while (el != NULL)
    {
        if (el->elements.size() >= minsize) {
            for (auto x : el->elements)
            {
                cout << x << " ";
            }
            cout << endl;
        }
        
        
        el = el->next;
    }
}