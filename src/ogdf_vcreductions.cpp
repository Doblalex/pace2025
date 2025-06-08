#include <ogdf/lpsolver/LPSolver.h>

#include "ogdf_instance.hpp"

bool Instance::reductionVCLP() {
	OGDF_ASSERT(isVCInstance());


	int cnt_cols = 0;
	int cnt_rows = 0;
	int cnt_nonzero = 0;
	ogdf::NodeArray<int> colIndices(G);
	ogdf::NodeArray<int> rowIndices(G);
	std::vector<ogdf::node> nodes;
	for (auto n : G.nodes) {
		if (!is_subsumed[n]) {
			colIndices[n] = nodes.size();
			nodes.push_back(n);
			cnt_cols++;
			cnt_nonzero += countCanDominate(n);
		}
		if (!is_dominated[n]) {
			rowIndices[n] = cnt_rows;
			cnt_rows++;
		}
	}

	ogdf::LPSolver::OptimizationGoal goal =
			ogdf::LPSolver::OptimizationGoal::Minimize; // goal of optimization (minimize or maximize)
	ogdf::Array<double> obj(cnt_cols); // objective function vector
	ogdf::Array<int> matrixBegin(cnt_cols); // matrixBegin[i] = begin of column i
	ogdf::Array<int> matrixCount(cnt_cols); // matrixCount[i] = number of nonzeroes in column i
	ogdf::Array<int> matrixIndex(cnt_nonzero); // matrixIndex[n] = index of matrixValue[n] in its column
	ogdf::Array<double> matrixValue(cnt_nonzero); // matrixValue[n] = non-zero value in matrix
	ogdf::Array<double> rightHandSide(cnt_rows); // right-hand side of LP constraints
	ogdf::Array<char> equationSense(cnt_rows); // 'E' ==   'G' >=   'L' <=
	ogdf::Array<double> lowerBound(cnt_cols); // lower bound of x[i]
	ogdf::Array<double> upperBound(cnt_cols); // upper bound of x[i]
	double optimum; // optimum value of objective function (if result is Optimal)
	ogdf::Array<double> x(cnt_cols);

	int index_nonzero = 0;
	for (auto n : G.nodes) {
		if (!is_subsumed[n]) {
			matrixBegin[colIndices[n]] = index_nonzero;
			matrixCount[colIndices[n]] = countCanDominate(n);
			lowerBound[colIndices[n]] = 0;
			upperBound[colIndices[n]] = 1;
			obj[colIndices[n]] = 1;
			forAllCanDominate(n, [&](ogdf::node adj) {
				matrixIndex[index_nonzero] = rowIndices[adj];
				matrixValue[index_nonzero] = 1;
				index_nonzero++;
				return true;
			});
		}
		if (!is_dominated[n]) {
			rightHandSide[rowIndices[n]] = 1;
			equationSense[rowIndices[n]] = 'G';
		}
	}

	ogdf::LPSolver solver;

	solver.optimize(goal, obj, matrixBegin, matrixCount, matrixIndex, matrixValue, rightHandSide,
			equationSense, lowerBound, upperBound, optimum, x);

	log << optimum << std::endl;

	for (int i = 0; i < cnt_cols; i++) {
		log << x[i] << std::endl;
	}

	return false;
}