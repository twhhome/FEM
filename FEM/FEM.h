#pragma once

#include "FEMCore.h"

class FEM
{
public:
	FEM(const char* filename);
	void solve();

private:
	int elementType;
	int nodeDOF;
	Array<Material> materials;
	Array<Section> sections;
	Array<Offset> offsets;
	Array<Node> nodes;
	Array<Element> elements;
	Array<Constraint> constraints;
	Array<Load> loads;
};

