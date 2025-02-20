#include "util.hpp"
#include "cxxopts.h"
#include "readinstance.hpp"

int main(int argc, char** argv) {
    auto instance = read_instance();
    auto adjList = instance.first;
    auto edgeSet = instance.second;

    cout<<adjList.size()<<" "<<edgeSet.size()<<endl;
}