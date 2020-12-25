#include "FEM.h"
#include "readFile.h"
#include "printToScreen.h"

FEM::FEM(const char* inputFileName, const char* outputFileName) : outputFileName(outputFileName) {
	readFile(inputFileName, nodeDOF, elementType, materials, sections, offsets, nodes, secondNodes, elements, constraints, loads);
}

void FEM::solve() {
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, elements, K);
	processConstraints(constraints, K);

	solveEqns(nodeDOF, nodes, K, loads);

	calElements(elementType, sections, nodeDOF, nodes, elements);
	calConstraintForce(nodeDOF, nodes, elements, constraints);

	outputSolution(outputFileName, elementType, nodeDOF, nodes, elements, constraints);

	printFinishMessage();
}