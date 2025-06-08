#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"

ogdf::Logger logger;

int main(int argc, char** argv) {
	logger.localLogLevel(ogdf::Logger::Level::Default);
	Instance I;
	if (argc != 2) {
		std::cerr << "Usage: " << argv[0] << " file.gr" << std::endl;
		return 1;
	}
	std::vector<ogdf::node> ID2node;
	{
		std::ifstream fin(argv[1]);
		if (!fin.good()) {
			std::cerr << "Error opening file " << argv[1] << std::endl;
			return 1;
		}
		I.read(fin, ID2node);
	}
	std::cout << "Instance with " << I.G.numberOfNodes() << " vertices and " << I.G.numberOfEdges()
			  << " edges" << std::endl;

	std::string line;
	while (line.empty() || line[0] == 'c') {
		getline(std::cin, line);
	}
	int cnt;
	std::istringstream iss(line);
	iss >> cnt;
	std::cout << I.type << " with " << cnt << " vertices" << std::endl;

	for (int i = 0; i < cnt; i++) {
		getline(std::cin, line);
		if (line.empty() || line[0] == 'c') {
			continue;
		}
		int u;
		std::istringstream iss(line);
		iss >> u;
		// std::cout << u << " ";
		I.addToDominatingSet(ID2node.at(u));
	}
	// std::cout << std::endl;

	int undom = 0;
	for (auto n : I.G.nodes) {
		if (!I.is_dominated[n]) {
			std::cerr << "Not dominated vertex " << I.node2ID[n] << ", in-neighs:";
			forAllInAdj(n, [&I](ogdf::adjEntry adj) {
				std::cerr << " " << I.node2ID[adj->twinNode()];
				return true;
			});
			std::cerr << std::endl;
			++undom;
		}
	}

	if (undom == 0) {
		std::cout << "Valid " << I.type << std::endl;
		return 0;
	} else {
		std::cerr << "Invalid " << I.type << " not dominating " << undom << " vertices!" << std::endl;
		return 2;
	}
}
