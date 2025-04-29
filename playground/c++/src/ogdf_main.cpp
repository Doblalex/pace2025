#include "ogdf_instance.hpp"
#include "ogdf_solver.hpp"
#include "ogdf_util.hpp"

ogdf::Logger logger;

ogdf::node internal::idn(ogdf::node n) { return n; }

ogdf::edge internal::ide(ogdf::edge n) { return n; }

#define OGDF_DEBUG

int main(int argc, char** argv) {
	logger.localLogLevel(ogdf::Logger::Level::Default);
	ogdf::Logger::globalLogLevel(ogdf::Logger::Level::Default);
	Instance I;
#ifdef OGDF_DEBUG
	std::vector<ogdf::node> ID2node;
	I.read(std::cin, ID2node);
	Instance I2;
	ogdf::NodeArray<ogdf::node> nMap(I.G, nullptr);
	ogdf::EdgeArray<ogdf::edge> eMap(I.G, nullptr);
	I2.G.insert(I.G, nMap, eMap);
	I2.initFrom(I, I.G.nodes, I.G.edges, nMap, eMap);
	for (auto& e : ID2node) {
		e = e == nullptr ? nullptr : nMap[e];
	}
#else
	I.read(std::cin);
#endif

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

#ifdef OGDF_DEBUG
	ogdf::Logger::globalLogLevel(ogdf::Logger::Level::Alarm);
	for (auto v : I.DS) {
		I2.addToDominatingSet(ID2node.at(v));
	}

	int undom = 0;
	for (auto n : I2.G.nodes) {
		if (!I2.is_dominated[n]) {
			std::cerr << "c Not dominated vertex " << I2.node2ID[n] << ", in-neighs:";
			forAllInAdj(n, [&I2](ogdf::adjEntry adj) {
				std::cerr << " " << I2.node2ID[adj->twinNode()];
				return true;
			});
			std::cerr << std::endl;
			++undom;
		}
	}

	if (undom == 0) {
		std::cerr << "c Valid DS" << std::endl;
		return 0;
	} else {
		std::cerr << "c Invalid DS not dominating " << undom << " vertices!" << std::endl;
		return 2;
	}
#endif
	return 0;
}
