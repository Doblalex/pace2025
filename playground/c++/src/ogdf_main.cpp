#include "ogdf_util.hpp"
#include "ogdf_instance.hpp"
#include "ogdf_solver.hpp"

ogdf::Logger logger;
ogdf::node internal::idn(ogdf::node n) { return n; }
ogdf::edge internal::ide(ogdf::edge n) { return n; }



int main(int argc, char **argv) {
    logger.localLogLevel(ogdf::Logger::Level::Default);
    auto start = std::chrono::system_clock::now();
    Instance I;
    I.read(std::cin);
    reduceAndSolve(I, 0);
    log<<"size dominating set: "<<I.DS.size()<<std::endl;
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> runtime = end - start;
}
