#pragma once

#include "FEMCore.h"

class FEM
{
public:
	FEM(const char* inputFileName, const char* outputFileName);
	void solve();

private:
	int elementType;
	int nodeDOF;
	Array<Material> materials;
	Array<Section> sections;
	Array<Offset> offsets;
	Array<Node> nodes;
	Array<Node> secondNodes; // 辅助节点，仅在三维梁单元时使用
	Array<Element> elements;
	Array<Constraint> constraints;
	Array<Load> loads;

	const char* outputFileName;
};