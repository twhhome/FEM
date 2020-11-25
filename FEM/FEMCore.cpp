#include "FEMCore.h"

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Rod> &rods, Matrix<double> &K) {
	K = Matrix<double>(nodeDOF * nodes.size());
	for (int i = 0; i < rods.size(); i++) {
		Matrix<double> temp = rods[i].T.trans() * rods[i].TK;
		Matrix<double> Ke = temp * rods[i].T;

		for (int j = 0; j < Ke.n; j++) {
			for (int k = 0; k < Ke.m; k++) {
				int Kj, Kk;
				Kj = nodeDOF * (((j < nodeDOF) ? rods[i].nodeNum1 : rods[i].nodeNum2) - 1) + j % nodeDOF;
				Kk = nodeDOF * (((k < nodeDOF) ? rods[i].nodeNum1 : rods[i].nodeNum2) - 1) + k % nodeDOF;
				K(Kj, Kk) += Ke(j, k);
			}
		}
	}
}

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Rod> &rods, MatrixIn1D &K) {
	// calculate pDiag
	K.pDiag = Matrix<int>(1, nodeDOF * nodes.size());
	K.pDiag(0) = 0;
	for (int i = 0; i < nodes.size(); i++) {
		int minIndex = std::numeric_limits<int>::max();
		for (int j = 0; j < rods.size(); j++) {
			if (i == rods[j].nodeNum1 - 1 || i == rods[j].nodeNum2 - 1) {
				if (MIN(rods[j].nodeNum1, rods[j].nodeNum2) < minIndex)
					minIndex = MIN(rods[j].nodeNum1, rods[j].nodeNum2);
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

	for (int i = 0; i < rods.size(); i++) {
		Matrix<double> temp = rods[i].T.trans() * rods[i].TK;
		Matrix<double> Ke = temp * rods[i].T;
		for (int j = 0; j < Ke.n; j++) {
			for (int k = 0; k < Ke.m; k++) {
				int Kj, Kk;
				Kj = nodeDOF * (((j < nodeDOF) ? rods[i].nodeNum1 : rods[i].nodeNum2) - 1) + j % nodeDOF;
				Kk = nodeDOF * (((k < nodeDOF) ? rods[i].nodeNum1 : rods[i].nodeNum2) - 1) + k % nodeDOF;
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
		K(constraints[i].dDirection - 1, constraints[i].dDirection - 1) = 1e25;//std::numeric_limits<double>::max();
	}
}

void processConstraints(Array<Constraint> &constraints, MatrixIn1D &K) {
	for (int i = 0; i < constraints.size(); i++) {
		K.set(constraints[i].dDirection - 1, constraints[i].dDirection - 1, std::numeric_limits<double>::max());
	}
}

void solve(int nodeDOF, Array<Node> &nodes, MatrixIn1D &K, Array<double> &loads) {
	Matrix<double> F(loads);
	F = F.trans();
	Matrix<double> d = cholesky(K, F);
	for (int i = 0; i < nodes.size(); i++) {
		for (int j = 0; j < nodeDOF; j++) {
			nodes[i].displacement(j) = d(i * nodeDOF + j);
		}
	}
}

void calRods(Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Rod> &rods) {
	for (int i = 0; i < rods.size(); i++) {
		Matrix<double> de(2 * nodeDOF, 1);
		for (int j = 0; j < nodeDOF; j++) {
			de(j) = nodes[rods[i].nodeNum1 - 1].displacement(j);
			de(j + nodeDOF) = nodes[rods[i].nodeNum2 - 1].displacement(j);
		}
		rods[i].de = de;
		rods[i].dee = rods[i].T * de;
		rods[i].fee = rods[i].TK * rods[i].dee;
		rods[i].fe = rods[i].T.trans() * rods[i].fee;
		rods[i].IF = rods[i].fee(1);
		rods[i].stress = rods[i].IF / sections[rods[i].sectionNum - 1].A;
	}
}

void calConstraintForce(int nodeDOF, Array<Node> &nodes, Array<Rod> &rods, Array<Constraint> &constraints) {
	Matrix<double> F(nodes.size() * nodeDOF, 1);
	for (int i = 0; i < rods.size(); i++) {
		for (int j = 0; j < nodeDOF; j++) {
			F((rods[i].nodeNum1 - 1) * nodeDOF + j) += rods[i].fe(j);
			F((rods[i].nodeNum2 - 1) * nodeDOF + j) += rods[i].fe(j + nodeDOF);
		}
	}
	for (int i = 0; i < constraints.size(); i++) {
		constraints[i].force = F(constraints[i].dDirection - 1);
	}
}