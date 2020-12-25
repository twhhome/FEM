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
	double Az, Ay, K; // ���Ǻ���ЧӦ�����ϵ��
} Section;

typedef struct {
	double a1, b1, c1;
	double a2, b2, c2;
} Offset; // ���Ľ����ƫ�Ĳ���

typedef struct {
	Matrix<double> pos; // ��������ϵ�µĽڵ�λ��
	Matrix<double> displacement; // ��������ϵ�µĽڵ�λ��
} Node;

typedef struct {
	int nodeNum1;
	int nodeNum2;
	int nodeNum3; // ����Ҫ�õ���3���ڵ㣬ֻ�����븨����������ά�����ʹ�ã�
	int materialNum; // ���ϱ��
	int sectionNum; // ������
	double L;
	int offsetNum; // ƫ��
	Matrix<double> Kee; // �ֲ�����ϵ���ڵ�����ϵ���µĸնȾ���
	Matrix<double> TT; // λ��ת�����󣨿���ƫ��ʱ�õ���
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

typedef struct {
	int dDirection; // �غ�ʩ�ӵ�λ�ƺ�
	double force; // �غ���
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