#include "FEMCore.h"
#include "readFile.h"
#include "printToScreen.h"

int main() {
	int elementType;
	int nodeDOF;
	Array<Material> materials;
	Array<Section> sections;
	Array<Offset> offsets;
	Array<Node> nodes;
	Array<Element> elements;
	Array<Constraint> constraints;
	Array<Load> loads;
	
	readFile("5.txt", nodeDOF, elementType, materials, sections, offsets, nodes, elements, constraints, loads);
	
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, elements, K);
	processConstraints(constraints, K);

	solve(nodeDOF, nodes, K, loads);
	
	calRods(elementType, sections, nodeDOF, nodes, elements);
	calConstraintForce(nodeDOF, nodes, elements, constraints);
	
	printSolution(elementType, nodeDOF, nodes, elements, constraints);
	
	printFinishMessage();
	return 0;
}