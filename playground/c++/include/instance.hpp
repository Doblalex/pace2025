#ifndef INSTANCE_H
#define INSTANCE_H

#include "util.hpp"

class Instance
{
    public:
    globalprops* props;
    VertexCount n,m;
    Graph* G;

    Instance(globalprops* props): props(props) {}

    Instance(VertexCount n, VertexCount m, globalprops* props) : n(n), m(m), props(props) {}

    ~Instance()
    {
        delete G;
    }

    void deleteVertex(VD v)
    {
        m -= boost::out_degree(v, *G);
        n--;
        boost::clear_vertex(v, *G);
        boost::remove_vertex(v, *G);
    }


    void deleteVertices(unordered_set<VD> &vertices)
    {
        n -= vertices.size();
        for (auto v : vertices)
        {
            m -= boost::out_degree(v, *G);
            clear_vertex(v, *G);
        }
        for (auto v: vertices)
        {
            boost::remove_vertex(v, *G);
        }
    }

    bool deleteVerticesWithIterator(unordered_set<VD>& vertices, Graph::vertex_iterator& iterator) {
        n -= vertices.size();
        bool increment = false;
        for (auto v : vertices)
        {
            m -= boost::out_degree(v, *G);
            clear_vertex(v, *G);
        }
        for (auto v: vertices)
        {
            if (*iterator == v) {
                iterator++;
                increment = true;
            }
            boost::remove_vertex(v, *G);
        }
        return increment;
    }

    void deleteVertices(list<VD> &vertices)
    {
        n -= vertices.size();
        for (auto v : vertices)
        {
            m -= boost::out_degree(v, *G);
            clear_vertex(v, *G);
        }
        for (auto v: vertices)
        {
            boost::remove_vertex(v, *G);
        }
    }

    void toDominatingSet(VD v, unordered_set<VD>& toDelete)
    {
        toDelete.insert(v);
        BGL_FORALL_ADJ_T(v, vd, *G, Graph)
        {
            if (!(*G)[vd].is_dominated)
            {
                BGL_FORALL_ADJ_T(vd, vd2, *G, Graph)
                {
                    (*G)[vd2].cnt_undominated_neighbors--; // vd is dominated
                    if ((*G)[vd2].cnt_undominated_neighbors == 0 && (*G)[vd2].is_dominated)
                    {
                        toDelete.insert(vd2);
                    }
                }
            }            
            (*G)[vd].is_dominated = true; // i am dominated
            if (!(*G)[v].is_dominated) {
                (*G)[vd].cnt_undominated_neighbors--; // v is dominated
            }            
            if ((*G)[vd].cnt_undominated_neighbors == 0)
            {
                toDelete.insert(vd);
            }
        }
        (*G)[v].is_dominated = true;
    }    

    VertexCount CntCanBeDominatingSet()
    {
        VertexCount cnt = 0;
        BGL_FORALL_VERTICES(v, *G, Graph) {
            if ((*G)[v].can_be_dominating_set) {
                cnt++;
            }
        }
        return cnt;
    }

    VertexCount CntNeedsDomination()
    {
        VertexCount cnt = 0;
        BGL_FORALL_VERTICES(v, *G, Graph) {
            if (!(*G)[v].is_dominated) {
                cnt++;
            }
        }
        return cnt;
    }
};

#endif