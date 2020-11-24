#pragma once

#include "MatrixIn1D.h"

typedef struct {
	Matrix<double> pos; // ��������ϵ�µĽڵ�λ��
	Matrix<double> displacement; // ��������ϵ�µĽڵ�λ��
} Node;

typedef struct {
	int nodeNum1;
	int nodeNum2;
	double E;
	double A;
	double L;
	Matrix<double> TK; // �ֲ�����ϵ�µĸնȾ���
	Matrix<double> T; // ����ת������
	Matrix<double> de; // ��������ϵ�µĽڵ�λ��
	Matrix<double> dee; // �ֲ�����ϵ�µĽڵ�λ��
	Matrix<double> fe; // ��������ϵ�µĽڵ���
	Matrix<double> fee; // �ֲ�����ϵ�µĽڵ���
	double IF; // ��Ԫ����(Internal Force)
	double stress; // ��ԪӦ��(Stress)
} Rod;

typedef struct {
	int dDirection; // Լ����λ�ƺ�
	double force; // Լ������
} Constraint;

void calStiffnessMatrix(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, MatrixIn1D &K);
void processConstraints(Matrix<Constraint> &constraints, MatrixIn1D &K);
void solve(int nodeDOF, Matrix<Node> &nodes, MatrixIn1D &K, Matrix<double> &loads);
void calRods(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods);
void calConstraintForce(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<Constraint> &constraints);