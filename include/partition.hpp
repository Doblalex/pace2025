#ifndef PARTITION_H
#define PARTITION_H


#include "instance.hpp"
#include "util.hpp"

struct PartitionElement {
	vector<Vertex> elements;
	PartitionElement* next = NULL;
	PartitionElement* last = NULL;
	vector<Vertex> refined;
};

struct PartitionRefinement {
	PartitionElement* start;
	PartitionElement* end;
	Graph* graph;
	globalprops* props;

	PartitionRefinement(Graph* graph, globalprops* props);

	~PartitionRefinement() {
		PartitionElement* el = start;
		while (el != NULL) {
			auto nextel = el->next;
			delete el;
			el = nextel;
		}
	}

	bool dorefine();

	void refine(vector<Vertex> arr);
	void print(int minsize = 1);
};

#endif