#include "FEMCore.h"

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Matrix<double> &K) {
	K = Matrix<double>(nodeDOF * nodes.size());
	for (int i = 0; i < elements.size(); i++) {
		Matrix<double> Ke = elements[i].T.trans() * elements[i].Kee * elements[i].T;

		for (int j = 0; j < Ke.n; j++) {
			for (int k = 0; k < Ke.m; k++) {
				int Kj, Kk;
				Kj = nodeDOF * (((j < nodeDOF) ? elements[i].nodeNum1 : elements[i].nodeNum2) - 1) + j % nodeDOF;
				Kk = nodeDOF * (((k < nodeDOF) ? elements[i].nodeNum1 : elements[i].nodeNum2) - 1) + k % nodeDOF;
				K(Kj, Kk) += Ke(j, k);
			}
		}
	}
}

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, MatrixIn1D &K) {
	// calculate pDiag
	K.pDiag = Matrix<int>(1, nodeDOF * nodes.size());
	K.pDiag(0) = 0;
	for (int i = 0; i < nodes.size(); i++) {
		int minIndex = std::numeric_limits<int>::max();
		for (int j = 0; j < elements.size(); j++) {
			if (i == elements[j].nodeNum1 - 1 || i == elements[j].nodeNum2 - 1) {
				if (MIN(elements[j].nodeNum1, elements[j].nodeNum2) < minIndex)
					minIndex = MIN(elements[j].nodeNum1, elements[j].nodeNum2);
			}
		}
		for (int j = 0; j < nodeDOF; j++) {
			if (i == 0 && j == 0)
				continue;
			int width = (i - (minIndex - 1)) * nodeDOF + j + 1;
			K.pDiag(i * nodeDOF + j) = K.pDiag(i * nodeDOF + j - 1) + width;
		}
	}
	K.data = Matrix<double>(1, K.pDiag(K.pDiag.m - 1) + 1);

	for (int i = 0; i < elements.size(); i++) {
		Matrix<double> Ke = elements[i].T.trans() * elements[i].Kee * elements[i].T;
		for (int j = 0; j < Ke.n; j++) {
			for (int k = 0; k < Ke.m; k++) {
				int Kj, Kk;
				Kj = nodeDOF * (((j < nodeDOF) ? elements[i].nodeNum1 : elements[i].nodeNum2) - 1) + j % nodeDOF;
				Kk = nodeDOF * (((k < nodeDOF) ? elements[i].nodeNum1 : elements[i].nodeNum2) - 1) + k % nodeDOF;
				if (Kj < Kk)
					continue;
				int index = K.pDiag(Kj) - (Kj - Kk);
				K.data(index) += Ke(j, k);
			}
		}
	}
}

void processConstraints(Array<Constraint> &constraints, Matrix<double> &K) {
	for (int i = 0; i < constraints.size(); i++) {
		K(constraints[i].dDirection - 1, constraints[i].dDirection - 1) = std::numeric_limits<double>::max();
	}
}

void processConstraints(Array<Constraint> &constraints, MatrixIn1D &K) {
	for (int i = 0; i < constraints.size(); i++) {
		K.set(constraints[i].dDirection - 1, constraints[i].dDirection - 1, std::numeric_limits<double>::max());
	}
}

template <typename T>
Matrix<T> cholesky(Matrix<T> &A, Matrix<T> &B) {
	assert(A.n == B.n);
	assert(B.m == 1);

	int n = A.n;

	Matrix<T> t(n, n), l(n, n), d(n, n);

	t[0][0] = A[0][0];
	d[0][0] = A[0][0];
	l[0][0] = 1;

	for (int i = 1; i < n; i++) {
		for (int j = 0; j <= i - 1; j++) {
			T temp = 0;
			for (int k = 0; k <= j - 1; k++) {
				temp += t[i][k] * l[j][k];
			}
			t[i][j] = A[i][j] - temp;
			l[i][j] = t[i][j] / d[j][j];
		}
		l[i][i] = 1;
		T temp = 0;
		for (int k = 0; k <= i - 1; k++) {
			temp += t[i][k] * l[i][k];
		}
		d[i][i] = A[i][i] - temp;
	}

	Matrix<T> x(n, 1), y(n, 1);
	y[0][0] = B[0][0];
	for (int i = 1; i < n; i++) {
		T temp = 0;
		for (int k = 0; k <= i - 1; k++) {
			temp += l[i][k] * y[k][0];
		}
		y[i][0] = B[i][0] - temp;
	}

	x[n - 1][0] = y[n - 1][0] / d[n - 1][n - 1];
	for (int i = n - 2; i >= 0; i--) {
		T temp = 0;
		for (int k = i + 1; k <= n - 1; k++) {
			temp += l[k][i] * x[k][0];
		}
		x[i][0] = y[i][0] / d[i][i] - temp;
	}

	return x;
}

Matrix<double> cholesky(MatrixIn1D &A, Matrix<double> &B) {
	assert(A.pDiag.m == B.n);
	assert(B.m == 1);

	MatrixIn1D l(A.pDiag);

	int n = A.pDiag.m;

	Matrix<double> d(n, 1);

	d(0) = A(0, 0);
	l.set(0, 0, 1);

	for (int i = 1; i < n; i++) {
		for (int j = i - (l.pDiag(i) - l.pDiag(i - 1)) + 1; j <= i - 1; j++) {
			double temp = 0;
			if (j != 0) {
				for (int k = MAX(i - (l.pDiag(i) - l.pDiag(i - 1)) + 1, j - (l.pDiag(j) - l.pDiag(j - 1)) + 1); k <= j - 1; k++) {
					temp += l(i, k) * d(k) * l(j, k);
				}
			}
			l.set(i, j, (A(i, j) - temp) / d(j));
		}
		l.set(i, i, 1);
		double temp = 0;
		for (int k = i - (l.pDiag(i) - l.pDiag(i - 1)) + 1; k <= i - 1; k++) {
			temp += l(i, k) * l(i, k) * d(k);
		}
		d(i) = A(i, i) - temp;
	}

	Matrix<double> x(n, 1), y(n, 1);
	y(0) = B(0);
	for (int i = 1; i < n; i++) {
		double temp = 0;
		for (int k = i - (l.pDiag(i) - l.pDiag(i - 1)) + 1; k <= i - 1; k++) {
			temp += l(i, k) * y(k);
		}
		y(i) = B(i) - temp;
	}

	x(n - 1) = y(n - 1) / d(n - 1);
	for (int i = n - 2; i >= 0; i--) {
		double temp = 0;
		for (int k = i + 1; k <= n - 1; k++) {
			temp += l(k, i) * x(k);
		}
		x(i) = y(i) / d(i) - temp;
	}

	return x;
}


void solveEqns(int nodeDOF, Array<Node> &nodes, MatrixIn1D &K, Array<Load> &loads) {
	Matrix<double> F(nodeDOF * nodes.size(), 1);
	for (int i = 0; i < loads.size(); i++) {
		F(loads[i].dDirection - 1) = loads[i].force;
	}

	Matrix<double> d = cholesky(K, F);
	for (int i = 0; i < nodes.size(); i++) {
		for (int j = 0; j < nodeDOF; j++) {
			if (fabs(d(i * nodeDOF + j)) < 1e-300)
				d(i * nodeDOF + j) = 0;
			nodes[i].displacement(j) = d(i * nodeDOF + j);
		}
	}
}

void calElements(int elementType, Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements) {
	for (int i = 0; i < elements.size(); i++) {
		Matrix<double> de(2 * nodeDOF, 1);
		for (int j = 0; j < nodeDOF; j++) {
			de(j) = nodes[elements[i].nodeNum1 - 1].displacement(j);
			de(j + nodeDOF) = nodes[elements[i].nodeNum2 - 1].displacement(j);
		}
		elements[i].de = de;
		elements[i].dee = elements[i].T * de;
		elements[i].fee = elements[i].Kee * elements[i].dee;
		elements[i].fe = elements[i].T.trans() * elements[i].fee;
		if (elementType == 1) {
			elements[i].IF = elements[i].fee(1);
			elements[i].stress = elements[i].IF / sections[elements[i].sectionNum - 1].A;
		}
	}
}

void calConstraintForce(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints) {
	Matrix<double> F(nodes.size() * nodeDOF, 1);
	for (int i = 0; i < elements.size(); i++) {
		for (int j = 0; j < nodeDOF; j++) {
			F((elements[i].nodeNum1 - 1) * nodeDOF + j) += elements[i].fe(j);
			F((elements[i].nodeNum2 - 1) * nodeDOF + j) += elements[i].fe(j + nodeDOF);
		}
	}
	for (int i = 0; i < constraints.size(); i++) {
		constraints[i].force = F(constraints[i].dDirection - 1);
	}
}