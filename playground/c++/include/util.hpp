#ifndef UTIL_H
#define UTIL_H


#include <bits/stdc++.h>
#include <random>
#include <tuple>
#include <string>
#include <stack>
#include <map>
#include <chrono>
#include <iostream>
#include <ctime>
#include <boost/config.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/function.hpp>
#include <boost/graph/iteration_macros.hpp>
#include <boost/graph/graph_utility.hpp>

#define MP make_pair
#define PB push_back
#define st first
#define nd second
#define NOT_PLACED INT_MAX
#define all(x) (x).begin(), (x).end()
#define FORE(itt, x) for (auto itt = (x).begin(); itt != (x).end(); ++itt)

#define MYDEBUG
#define MYLOCAL

// #ifdef MYDEBUG
// #undef NDEBUG
// #endif

using namespace std;
using namespace std::chrono;

static std::random_device rd;
static std::mt19937 rng{rd()};

template <typename TH>
void _dbg(const char *sdbg, TH h) { cout << sdbg << "=" << h << "\n"; }
template <typename TH, typename... TA>
void _dbg(const char *sdbg, TH h, TA... t)
{
    while (*sdbg != ',')
        cout << *sdbg++;
    cout << "=" << h << ",";
    _dbg(sdbg + 1, t...);
}
#ifdef MYLOCAL
#define debug(...) _dbg(#__VA_ARGS__, __VA_ARGS__)
#define debugv(x)             \
    {                         \
        {cout << #x << " = "; \
    FORE(itt, (x))            \
    cout << *itt << ", ";     \
    cout << "\n";             \
    }                         \
    }
#else
#define debug(...)
#define debugv(x)
#endif

#if defined(__GNUC__)

#define MY_LIB_IGNORE_DEPRECATED_BEGIN \
    _Pragma("GCC diagnostic push")     \
        _Pragma("GCC diagnostic ignored \"-Wdeprecated-declarations\"")

#define MY_LIB_IGNORE_DEPRECATED_END \
    _Pragma("GCC diagnostic pop")

#elif defined(_MSC_VER)

#define MY_LIB_IGNORE_DEPRECATED_BEGIN \
    _Pragma("warning(push)")           \
        _Pragma("warning(disable : 4996)")

#define MY_LIB_IGNORE_DEPRECATED_END \
    _Pragma("warning(pop)")

#else

#define MY_LIB_IGNORE_DEPRECATED_BEGIN
#define MY_LIB_IGNORE_DEPRECATED_END

#endif

template <class T>
ostream &operator<<(ostream &out, vector<T> vec)
{
    out << "(";
    for (auto &v : vec)
        out << v << ", ";
    return out << ")";
}
template <class T>
ostream &operator<<(ostream &out, set<T> vec)
{
    out << "(";
    for (auto &v : vec)
        out << v << ", ";
    return out << ")";
}

class NotImplemented : public std::logic_error
{
public:
    NotImplemented() : std::logic_error("Function not yet implemented") {};
};

template <typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator &g)
{
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template <typename Iter>
Iter select_randomly(Iter start, Iter end)
{
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

typedef boost::hash<pair<u_int32_t, u_int32_t>> pair_hash;

class BufferStdout
{
public:
    // the collector string is used for collecting the output to stdout
    BufferStdout(std::string &collector) : m_collector(collector),
                                           fp(std::fopen("output.txt", "w"))
    {
        if (fp == nullptr)
            throw std::runtime_error(std::strerror(errno));
        std::swap(stdout, fp); // swap stdout and the temp file
    }

    ~BufferStdout()
    {
        std::swap(stdout, fp); // swap back
        std::fclose(fp);

        // read the content of the temp file into m_collector
        if (std::ifstream is("output.txt"); is)
        {
            m_collector.append(std::istreambuf_iterator<char>(is),
                               std::istreambuf_iterator<char>{});
        }
        std::remove("output.txt"); // cleanup
    }

private:
    std::string &m_collector;
    std::FILE *fp;
};

typedef u_int32_t Vertex;
typedef u_int32_t VertexCount;
typedef pair<Vertex, Vertex> Edge;
typedef list<Vertex> VertexList;

struct vertex_info;
typedef boost::adjacency_list<boost::listS, boost::listS, boost::undirectedS, boost::property<boost::vertex_index_t, size_t, vertex_info>> Graph;
typedef boost::filtered_graph<Graph, boost::keep_all, boost::keep_all> FilteredGraph;
using VD = boost::graph_traits<Graph>::vertex_descriptor;
// #define IteratorVertices(x) boost::vertices(x)
// #define MyVertices(x) boost::make_iterator_range(boost::vertices(x))
// #define IteratorEdges(x) boost::edges(x)
// #define Neighbors(g,v) boost::make_iterator_range(boost::adjacent_vertices(v, g))

struct vertex_info
{
    Vertex id;
    bool is_dominated = false;
    bool can_be_dominating_set = true;
    VertexCount cnt_undominated_neighbors = 0;
    VD forward;
};

struct globalprops {
    int n;
    int m;
    vector<VD> id_to_vertex;
};


#endif