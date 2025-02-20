#pragma once

#include <bits/stdc++.h>
#include <random>
#include <tuple>
#include <string>
#include <stack>
#include <boost/functional/hash.hpp>

#define MP make_pair
#define PB push_back
#define st first
#define nd second
#define NOT_PLACED INT_MAX

// #define MYDEBUG
// #define MYLOCAL

#ifdef MYDEBUG
#undef NDEBUG
#endif

using namespace std;
using namespace std::chrono;

static std::random_device rd; 
static std::mt19937 rng{rd()};

template <typename TH>
void _dbg(const char *sdbg, TH h) { cerr << sdbg << "=" << h << "\n"; }
template <typename TH, typename... TA>
void _dbg(const char *sdbg, TH h, TA... t)
{
    while (*sdbg != ',')
        cerr << *sdbg++;
    cerr << "=" << h << ",";
    _dbg(sdbg + 1, t...);
}
#ifdef MYLOCAL
#define debug(...) _dbg(#__VA_ARGS__, __VA_ARGS__)
#define debugv(x)                 \
    {                             \
        {                         \
            cerr << #x << " = ";  \
            FORE(itt, (x))        \
            cerr << *itt << ", "; \
            cerr << "\n";         \
        }                         \
    }
#else
#define debug(...) (__VA_ARGS__)
#define debugv(x)
#define cerr \
    if (0)   \
    cout
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
    NotImplemented() : std::logic_error("Function not yet implemented") { };
};

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

typedef boost::hash<pair<u_int32_t, u_int32_t>> pair_hash;

class BufferStdout {
public:
    // the collector string is used for collecting the output to stdout
    BufferStdout (std::string& collector) :
        m_collector(collector),
        fp(std::fopen("output.txt", "w"))
    {
        if(fp == nullptr) throw std::runtime_error(std::strerror(errno));
        std::swap(stdout, fp); // swap stdout and the temp file
    }

    ~BufferStdout () {
        std::swap(stdout, fp); // swap back
        std::fclose(fp);

        // read the content of the temp file into m_collector
        if(std::ifstream is("output.txt"); is) {
            m_collector.append(std::istreambuf_iterator<char>(is),
                               std::istreambuf_iterator<char>{});
        }
        std::remove("output.txt"); // cleanup
    }

private:
    std::string& m_collector;
    std::FILE* fp;
};

typedef unsigned int Vertex;
typedef pair<Vertex, Vertex> Edge;
typedef unordered_set<Vertex> VertexSet;
typedef vector<Vertex> VectorList;
typedef vector<unordered_set<Vertex>> AdjacencyList;
typedef vector<Edge> EdgeList;
typedef unordered_set<Edge, pair_hash> EdgeSet;