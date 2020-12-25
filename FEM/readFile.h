#pragma once

#include "FEMCore.h"

void readFile(const char* filename, int &nodeDOF, int &elementType, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets,
						Array<Node> &nodes, Array<Node> &secondNodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads);

void readRodsFromFile(FILE *fp, int &nodeDOF, Array<Material> &materials, Array<Section> &sections, 
						Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads);

void readBeamsFromFile(FILE *fp, int &nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets, 
						Array<Node> &nodes, Array<Node> &secondNodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads);
