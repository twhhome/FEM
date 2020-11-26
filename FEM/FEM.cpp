#include "FEM.h"
#include "readFile.h"
#include "printToScreen.h"

FEM::FEM(const char* filename) {
	readFile(filename, nodeDOF, elementType, materials, sections, offsets, nodes, elements, constraints, loads);
}

void FEM::solve() {
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, elements, K);
	processConstraints(constraints, K);

	solveEqns(nodeDOF, nodes, K, loads);

	calRods(elementType, sections, nodeDOF, nodes, elements);
	calConstraintForce(nodeDOF, nodes, elements, constraints);

	printSolution(elementType, nodeDOF, nodes, elements, constraints);

	printFinishMessage();
}