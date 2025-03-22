#include "ogdf_instance.hpp"
#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"

ogdf::Logger logger;

ogdf::node internal::idn(ogdf::node n) { return n; }

ogdf::edge internal::ide(ogdf::edge n) { return n; }

int main(int argc, char** argv) {
	logger.localLogLevel(ogdf::Logger::Level::Default);
	Instance I;
	I.read(std::cin);
	auto start = std::chrono::high_resolution_clock::now();
	reduceAndSolve(I, 0);
	auto end = std::chrono::high_resolution_clock::now();

	std::cerr << "c DS solution size:\n" << I.DS.size() << "\nc <DS vertices>:" << std::endl;
	for (auto v : I.DS) {
		std::cerr << v << "\n";
	}
	std::cerr << "c </DS vertices>\nc DS solution size: " << I.DS.size() << "\nc solve time: "
			  << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms"
			  << std::endl;
	return 0;
}
