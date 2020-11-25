#pragma once

#include "MatrixIn1D.h"
#include "Array.h"

typedef struct {
	double E;
	double poisson;
} Material;

typedef struct {
	double A;
} Section;

typedef struct {
	Matrix<double> pos; // 整体坐标系下的节点位置
	Matrix<double> displacement; // 整体坐标系下的节点位移
} Node;

typedef struct {
	int nodeNum1;
	int nodeNum2;
	int materialNum; // 材料编号
	int sectionNum; // 截面编号
	double L;
	Matrix<double> TK; // 局部坐标系下的刚度矩阵
	Matrix<double> T; // 坐标转换矩阵
	Matrix<double> de; // 整体坐标系下的节点位移
	Matrix<double> dee; // 局部坐标系下的节点位移
	Matrix<double> fe; // 整体坐标系下的节点力
	Matrix<double> fee; // 局部坐标系下的节点力
	double IF; // 单元内力(Internal Force)
	double stress; // 单元应力(Stress)
} Element;

typedef struct {
	int dDirection; // 约束的位移号
	double force; // 约束反力
} Constraint;

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Matrix<double> &K);
void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, MatrixIn1D &K);
void processConstraints(Array<Constraint> &constraints, Matrix<double> &K);
void processConstraints(Array<Constraint> &constraints, MatrixIn1D &K);
void solve(int nodeDOF, Array<Node> &nodes, MatrixIn1D &K, Array<double> &loads);
void calRods(Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements);
void calConstraintForce(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints);