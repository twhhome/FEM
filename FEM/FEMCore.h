#pragma once

#include "MatrixIn1D.h"

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
} Rod;

typedef struct {
	int dDirection; // 约束的位移号
	double force; // 约束反力
} Constraint;

void calStiffnessMatrix(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<double> &K);
void calStiffnessMatrix(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, MatrixIn1D &K);
void processConstraints(Matrix<Constraint> &constraints, Matrix<double> &K);
void processConstraints(Matrix<Constraint> &constraints, MatrixIn1D &K);
void solve(int nodeDOF, Matrix<Node> &nodes, MatrixIn1D &K, Matrix<double> &loads);
void calRods(Matrix<Section> &sections, int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods);
void calConstraintForce(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<Constraint> &constraints);