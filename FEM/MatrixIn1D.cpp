#include "MatrixIn1D.h"

double MatrixIn1D::operator()(int i, int j) {
	int x = MAX(i, j);
	int y = MIN(i, j);
	if (x == 0 && y == 0)
		return data(0);
	if (x - (pDiag(x) - pDiag(x - 1)) + 1 > y)
		return 0;
	else
		return data(y - x + pDiag(x));
}

void MatrixIn1D::set(int i, int j, double n) {
	int x = MAX(i, j);
	int y = MIN(i, j);
	if (x == 0 && y == 0)
		data(0) = n;
	else if (x - (pDiag(x) - pDiag(x - 1)) + 1 > y)
		return;
	else
		data(y - x + pDiag(x)) = n;
}