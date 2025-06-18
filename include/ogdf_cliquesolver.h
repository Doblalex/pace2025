#pragma once

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

#include "ogdf_instance.hpp"
#include "ogdf_util.hpp"
#include <boost/asio/io_service.hpp>
#include <boost/process.hpp>
#include <boost/process/async.hpp>

namespace bp = boost::process;
namespace fs = std::filesystem;

std::string get_executable_directory() {
	char result[PATH_MAX];
	ssize_t count = readlink("/proc/self/exe", result, PATH_MAX);
	if (count == -1) {
		throw std::runtime_error("Failed to read /proc/self/exe");
	}
	fs::path exe_path(std::string(result, count));
	return exe_path.parent_path().string();
}

bool solveMISInstanceWithCliqueSolver(Instance& I, long limit_branches, bool docheck = false,
		long limit_seconds = 3600) {
	OGDF_ASSERT(I.isVCInstance());

	int numvertices = 1;
	ogdf::NodeArray<int> nodeToIndex(I.G);
	std::vector<ogdf::node> indexToNode;
	indexToNode.push_back(nullptr);
	for (auto node : I.G.nodes) {
		if (!I.is_subsumed[node]) {
			nodeToIndex[node] = numvertices;
			indexToNode.push_back(node);
			numvertices++;
		}
	}
	int numedges;
	for (int i = 1; i < numvertices; i++) {
		I.forAllCanDominate(indexToNode[i], [&](ogdf::node adj) {
			I.forAllCanBeDominatedBy(adj, [&](ogdf::node adj2) {
				if (adj2 != indexToNode[i]) {
					numedges++;
				}
				return true;
			});
			return true;
		});
	}

	bp::opstream ofs;
	ofs << "p td " << (numvertices - 1) << " " << numedges << std::endl;
	for (int i = 1; i < numvertices; i++) {
		I.forAllCanDominate(indexToNode[i], [&](ogdf::node adj) {
			I.forAllCanBeDominatedBy(adj, [&](ogdf::node adj2) {
				if (adj2 != indexToNode[i]) {
					ofs << i << " " << nodeToIndex[adj2] << std::endl;
				}
				return true;
			});
			return true;
		});
	}

	log << "Trying to solve VC instance with " << numvertices - 1 << " vertices and " << numedges
		<< " edges for " << limit_seconds << "s" << std::endl;
	boost::asio::io_service ios;
	std::future<std::string> data;
	bp::child c("/usr/bin/timeout " + std::to_string(limit_seconds) + " "
					+ get_executable_directory() + "/ext/peaty/solve_vc -q",
			bp::std_out > data, bp::std_err > stderr, bp::std_in < ofs, ios);
	ios.run();
	c.wait();
	if (c.exit_code() != 0) {
		log << "VC solver returned error " << c.exit_code() << std::endl;
		return false;
	}

	std::istringstream is(data.get());
	std::string line;
	int size = 0;
	while (is && std::getline(is, line)) {
		if (line.empty() || line[0] == 'c' || line[0] == 's') {
			continue;
		}
		I.DS.insert(I.node2ID[indexToNode[std::stoi(line)]]);
		size++;
	}
	log << "VC solver found solution of size " << size << std::endl;
	return true;
}
