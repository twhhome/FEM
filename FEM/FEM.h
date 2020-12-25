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
	Array<Node> secondNodes; // �����ڵ㣬������ά����Ԫʱʹ��
	Array<Element> elements;
	Array<Constraint> constraints;
	Array<Load> loads;

	const char* outputFileName;
};