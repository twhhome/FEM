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
	Matrix<double> pos; // ��������ϵ�µĽڵ�λ��
	Matrix<double> displacement; // ��������ϵ�µĽڵ�λ��
} Node;

typedef struct {
	int nodeNum1;
	int nodeNum2;
	int materialNum; // ���ϱ��
	int sectionNum; // ������
	double L;
	Matrix<double> TK; // �ֲ�����ϵ�µĸնȾ���
	Matrix<double> T; // ����ת������
	Matrix<double> de; // ��������ϵ�µĽڵ�λ��
	Matrix<double> dee; // �ֲ�����ϵ�µĽڵ�λ��
	Matrix<double> fe; // ��������ϵ�µĽڵ���
	Matrix<double> fee; // �ֲ�����ϵ�µĽڵ���
	double IF; // ��Ԫ����(Internal Force)
	double stress; // ��ԪӦ��(Stress)
} Element;

typedef struct {
	int dDirection; // Լ����λ�ƺ�
	double force; // Լ������
} Constraint;

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Matrix<double> &K);
void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, MatrixIn1D &K);
void processConstraints(Array<Constraint> &constraints, Matrix<double> &K);
void processConstraints(Array<Constraint> &constraints, MatrixIn1D &K);
void solve(int nodeDOF, Array<Node> &nodes, MatrixIn1D &K, Array<double> &loads);
void calRods(Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements);
void calConstraintForce(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints);