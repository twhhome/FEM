#pragma once

#include "FEMCore.h"

void printParameters(int elementType, int nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets,
						Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads);

void printSolution(int elementType, int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints);

void printFinishMessage();