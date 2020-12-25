#pragma once

#include "MatrixIn1D.h"
#include "Array.h"

typedef struct {
	double E;
	double poisson;
} Material;

typedef struct {
	double A;
	double Ix, Iy, Iz;
	double Az, Ay, K; // 考虑横向效应引入的系数
} Section;

typedef struct {
	double a1, b1, c1;
	double a2, b2, c2;
} Offset; // 梁的截面的偏心参数

typedef struct {
	Matrix<double> pos; // 整体坐标系下的节点位置
	Matrix<double> displacement; // 整体坐标系下的节点位移
} Node;

typedef struct {
	int nodeNum1;
	int nodeNum2;
	int nodeNum3; // 梁需要用到第3个节点，只是用与辅助（仅在三维情况下使用）
	int materialNum; // 材料编号
	int sectionNum; // 截面编号
	double L;
	int offsetNum; // 偏心
	Matrix<double> Kee; // 局部坐标系（节点坐标系）下的刚度矩阵
	Matrix<double> TT; // 位移转换矩阵（考虑偏心时用到）
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

typedef struct {
	int dDirection; // 载荷施加的位移号
	double force; // 载荷力
 } Load;

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Matrix<double> &K);
void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, MatrixIn1D &K);

void processConstraints(Array<Constraint> &constraints, Matrix<double> &K);
void processConstraints(Array<Constraint> &constraints, MatrixIn1D &K);

template <typename T> Matrix<T> cholesky(Matrix<T> &A, Matrix<T> &B);
Matrix<double> cholesky(MatrixIn1D &A, Matrix<double> &B);

void solveEqns(int nodeDOF, Array<Node> &nodes, MatrixIn1D &K, Array<Load> &loads);

void calElements(int elementType, Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements);

void calConstraintForce(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints);