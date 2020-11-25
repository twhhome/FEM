#include "FEMCore.h"

void calStiffnessMatrix(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Matrix<double> &K) {
	K = Matrix<double>(nodeDOF * nodes.size());
	for (int i = 0; i < elements.size(); i++) {
		Matrix<double> temp = elements[i].T.trans() * elements[i].TK;
		Matrix<double> Ke = temp * elements[i].T;

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
		Matrix<double> temp = elements[i].T.trans() * elements[i].TK;
		Matrix<double> Ke = temp * elements[i].T;
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

void calRods(Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements) {
	for (int i = 0; i < elements.size(); i++) {
		Matrix<double> de(2 * nodeDOF, 1);
		for (int j = 0; j < nodeDOF; j++) {
			de(j) = nodes[elements[i].nodeNum1 - 1].displacement(j);
			de(j + nodeDOF) = nodes[elements[i].nodeNum2 - 1].displacement(j);
		}
		elements[i].de = de;
		elements[i].dee = elements[i].T * de;
		elements[i].fee = elements[i].TK * elements[i].dee;
		elements[i].fe = elements[i].T.trans() * elements[i].fee;
		elements[i].IF = elements[i].fee(1);
		elements[i].stress = elements[i].IF / sections[elements[i].sectionNum - 1].A;
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