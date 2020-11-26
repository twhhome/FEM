#pragma once

#include "Matrix.h"

class MatrixIn1D {
public:
	Matrix<int> pDiag; // 主角元在一维储存中的编号
	Matrix<double> data;

	MatrixIn1D() {}

	MatrixIn1D(Matrix<int> pDiag) {
		this->pDiag = pDiag;
		this->data = Matrix<double>(1, this->pDiag(this->pDiag.m - 1) + 1);
	}

	double operator()(int i, int j);

	void set(int i, int j, double n);
};