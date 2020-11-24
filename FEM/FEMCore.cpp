#include "FEMCore.h"

void calStiffnessMatrix(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, MatrixIn1D &K) {
	// calculate pDiag
	K.pDiag = Matrix<int>(1, nodeDOF * nodes.n);
	K.pDiag(0) = 0;
	for (int i = 0; i < nodes.n; i++) {
		int minIndex = std::numeric_limits<int>::max();
		for (int j = 0; j < rods.n; j++) {
			if (i == rods(j).nodeNum1 - 1 || i == rods(j).nodeNum2 - 1) {
				if (MIN(rods(j).nodeNum1, rods(j).nodeNum2) < minIndex)
					minIndex = MIN(rods(j).nodeNum1, rods(j).nodeNum2);
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

	for (int i = 0; i < rods.n; i++) {
		Matrix<double> temp = rods(i).T.trans() * rods(i).TK;
		Matrix<double> Ke = temp * rods(i).T;
		for (int j = 0; j < Ke.n; j++) {
			for (int k = 0; k < Ke.m; k++) {
				int Kj, Kk;
				Kj = nodeDOF * (((j < nodeDOF) ? rods(i).nodeNum1 : rods(i).nodeNum2) - 1) + j % nodeDOF;
				Kk = nodeDOF * (((k < nodeDOF) ? rods(i).nodeNum1 : rods(i).nodeNum2) - 1) + k % nodeDOF;
				if (Kj < Kk)
					continue;
				int index = K.pDiag(Kj) - (Kj - Kk);
				K.data(index) += Ke(j, k);
			}
		}
	}
}

void processConstraints(Matrix<Constraint> &constraints, MatrixIn1D &K) {
	for (int i = 0; i < constraints.n; i++) {
		K.set(constraints(i).dDirection - 1, constraints(i).dDirection - 1, std::numeric_limits<double>::max());
	}
}

void solve(int nodeDOF, Matrix<Node> &nodes, MatrixIn1D &K, Matrix<double> &loads) {
	Matrix<double> d = cholesky(K, loads);
	for (int i = 0; i < nodes.n; i++) {
		for (int j = 0; j < nodeDOF; j++) {
			nodes(i).displacement(j) = d(i * nodeDOF + j);
		}
	}
}

void calRods(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods) {
	for (int i = 0; i < rods.n; i++) {
		Matrix<double> de(2 * nodeDOF, 1);
		for (int j = 0; j < nodeDOF; j++) {
			de(j) = nodes(rods(i).nodeNum1 - 1).displacement(j);
			de(j + nodeDOF) = nodes(rods(i).nodeNum2 - 1).displacement(j);
		}
		rods(i).de = de;
		rods(i).dee = rods(i).T * de;
		rods(i).fee = rods(i).TK * rods(i).dee;
		rods(i).fe = rods(i).T.trans() * rods(i).fee;
		rods(i).IF = rods(i).fee(1);
		rods(i).stress = rods(i).IF / rods(i).A;
	}
}

void calConstraintForce(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<Constraint> &constraints) {
	Matrix<double> F(nodes.n * nodeDOF, 1);
	for (int i = 0; i < rods.n; i++) {
		for (int j = 0; j < nodeDOF; j++) {
			F((rods(i).nodeNum1 - 1) * nodeDOF + j) += rods(i).fe(j);
			F((rods(i).nodeNum2 - 1) * nodeDOF + j) += rods(i).fe(j + nodeDOF);
		}
	}
	for (int i = 0; i < constraints.n; i++) {
		constraints(i).force = F(constraints(i).dDirection - 1);
	}
}