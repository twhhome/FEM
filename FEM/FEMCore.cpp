#include "FEMCore.h"

double calLength(Node &node1, Node &node2) {
	Matrix<double> delta = node1.pos - node2.pos;
	return sqrt((delta.trans() * delta)(0));
}

void readRodsFromFile(FILE *fp, int &nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads) {
	fscanf(fp, "%d", &nodeDOF);
	
	int materialNum;
	fscanf(fp, "%d", &materialNum);
	materials = Array<Material>(materialNum);
	for (int i = 0; i < materialNum; i++) {
		double E;
		fscanf(fp, "%lf", &E);
		materials[i].E = E;
	}

	int sectionNum;
	fscanf(fp, "%d", &sectionNum);
	sections = Array<Section>(sectionNum);
	for (int i = 0; i < sectionNum; i++) {
		double A;
		fscanf(fp, "%lf", &A);
		sections[i].A = A;
	}

	int nodeNum;
	fscanf(fp, "%d", &nodeNum);
	nodes = Array<Node>(nodeNum);
	for (int i = 0; i < nodeNum; i++) {
		nodes[i].pos = Matrix<double>(nodeDOF, 1);
		nodes[i].displacement = Matrix<double>(nodeDOF, 1);
		for (int j = 0; j < nodeDOF; j++) {
			fscanf(fp, "%lf", &(nodes[i].pos(j)));
		}
	}

	int elementNum;
	fscanf(fp, "%d", &elementNum);
	elements = Array<Element>(elementNum);
	int nodeNum1, nodeNum2, material, section;
	double E, A, L;
	Matrix<double> TK(2, 2), T(2, 2 * nodeDOF);
	for (int i = 0; i < elementNum; i++) {
		fscanf(fp, "%d%d%d%d", &nodeNum1, &nodeNum2, &material, &section);

		E = materials[material - 1].E;
		A = sections[section - 1].A;

		Node node1 = nodes[nodeNum1 - 1];
		Node node2 = nodes[nodeNum2 - 1];
		L = calLength(node1, node2);

		TK(0, 0) = TK(1, 1) = E * A / L;
		TK(0, 1) = TK(1, 0) = -E * A / L;

		Matrix<double> delta = node2.pos - node1.pos;
		delta = delta / L;
		for (int j = 0; j < nodeDOF; j++) {
			T(0, j) = T(1, j + nodeDOF) = delta(j);
		}

		elements[i].nodeNum1 = nodeNum1;
		elements[i].nodeNum2 = nodeNum2;
		elements[i].materialNum = material;
		elements[i].sectionNum = section;
		elements[i].L = L;
		elements[i].Kee = TK;
		elements[i].T = T;
	}

	int constraintNum;
	fscanf(fp, "%d", &constraintNum);
	constraints = Array<Constraint>(constraintNum);
	for (int i = 0; i < constraintNum; i++) {
		fscanf(fp, "%d", &(constraints[i].dDirection));
	}

	int loadNum;
	fscanf(fp, "%d", &loadNum);
	loads = Array<Load>(loadNum);
	for (int i = 0; i < loadNum; i++) {
		fscanf(fp, "%d%lf", &(loads[i].dDirection), &(loads[i].force));
	}
}


void readBeamsFromFile(FILE *fp, int &nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads) {
	fscanf(fp, "%d", &nodeDOF);
	
	int materialNum;
	fscanf(fp, "%d", &materialNum);
	materials = Array<Material>(materialNum);
	for (int i = 0; i < materialNum; i++) {
		double E, poisson;
		fscanf(fp, "%lf%lf", &E, &poisson);
		materials[i].E = E;
		materials[i].poisson = poisson;
	}

	int sectionNum;
	fscanf(fp, "%d", &sectionNum);
	sections = Array<Section>(sectionNum);
	for (int i = 0; i < sectionNum; i++) {
		double A, Ix, Iy, Iz, Ay, Az, K;
		if (nodeDOF == 3) { // 二维
			fscanf(fp, "%lf%lf%lf%lf", &A, &Iz, &Ay, &K);
			sections[i].A = A;
			sections[i].Iz = Iz;
			sections[i].Ay = Ay;
			sections[i].K = K;
		}
		else if (nodeDOF == 6) { // 三维
			fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf", &A, &Ix, &Iy, &Iz, &Ay, &Az, &K);
			sections[i].A = A;
			sections[i].Ix = Ix;
			sections[i].Iy = Iy;
			sections[i].Iz = Iz;
			sections[i].Ay = Ay;
			sections[i].Az = Az;
			sections[i].K = K;
		}
	}

	int offsetNum;
	fscanf(fp, "%d", &offsetNum);
	offsets = Array<Offset>(offsetNum);
	for (int i = 0; i < offsetNum; i++) {
		double a1, b1, c1, a2, b2, c2;
		if (nodeDOF == 3) { // 二维
			fscanf(fp, "%lf%lf%lf%lf", &a1, &b1, &a2, &b2);
			offsets[i].a1 = a1;
			offsets[i].b1 = b1;
			offsets[i].a2 = a2;
			offsets[i].b2 = b2;
		}
		else if (nodeDOF == 6) { // 三维
			fscanf(fp, "%lf%lf%lf%lf%lf%lf", &a1, &b1, &c1, &a2, &b2, &c2);
			offsets[i].a1 = a1;
			offsets[i].b1 = b1;
			offsets[i].c1 = c1;
			offsets[i].a2 = a2;
			offsets[i].b2 = b2;
			offsets[i].c2 = c2;
		}
	}

	int nodeNum;
	fscanf(fp, "%d", &nodeNum);
	nodes = Array<Node>(nodeNum);
	for (int i = 0; i < nodeNum; i++) {
		nodes[i].pos = Matrix<double>((nodeDOF == 3) ? 2 : 3, 1);
		nodes[i].displacement = Matrix<double>(nodeDOF, 1);
		for (int j = 0; j < ((nodeDOF == 3) ? 2 : 3); j++) {
			fscanf(fp, "%lf", &(nodes[i].pos(j)));
		}
	}

	int elementNum;
	fscanf(fp, "%d", &elementNum);
	elements = Array<Element>(elementNum);
	int nodeNum1, nodeNum2, nodeNum3, material, section, offset;
	double E, poisson, G;
	double A, Ix, Iy, Iz, Az, Ay, K;
	double by, bz;
	double a1, b1, c1, a2, b2, c2;
	double L;
	Matrix<double> TK, TT, Kee, T; // TK为主轴坐标系下的刚度矩阵
	for (int i = 0; i < elementNum; i++) {
		if (nodeDOF == 3) { // 二维
			fscanf(fp, "%d%d%d%d%d", &nodeNum1, &nodeNum2, &material, &section, &offset);

			E = materials[material - 1].E;
			poisson = materials[material - 1].poisson;
			G = E / (2 * (1 + poisson));

			A = sections[section - 1].A;
			Iz = sections[section - 1].Iz;
			Ay = sections[section - 1].Ay;
			K = sections[section - 1].K;

			a1 = offsets[offset - 1].a1;
			b1 = offsets[offset - 1].b1;
			a2 = offsets[offset - 1].a2;
			b2 = offsets[offset - 1].b2;

			Node node1 = nodes[nodeNum1 - 1];
			Node node2 = nodes[nodeNum2 - 1];
			L = calLength(node1, node2) + a1 - a2;

			by = 12 * K * E * Iz / G / Ay / L / L;

			TK = Matrix<double>(6, 6);
			TK(0, 0) = TK(3, 3) = E * A / L;
			TK(1, 1) = TK(4, 4) = 12 * E * Iz / (1 + by) / L / L / L;
			TK(2, 2) = TK(5, 5) = (4 + by) * E * Iz / (1 + by) / L;
			TK(2, 1) = TK(1, 2) = TK(5, 1) = TK(1, 5) = 6 * E * Iz / (1 + by) / L / L;
			TK(3, 0) = TK(0, 3) = -E * A / L;
			TK(4, 1) = TK(1, 4) = -12 * E * Iz / (1 + by) / L / L / L;
			TK(4, 2) = TK(2, 4) = TK(5, 4) = TK(4, 5) = -6 * E * Iz / (1 + by) / L / L;
			TK(5, 2) = TK(2, 5) = (2 - by) * E * Iz / (1 + by) / L;

			TT = Matrix<double>(6, 6);
			for (int j = 0; j < 6; j++)
				TT(j, j) = 1;
			TT(1, 2) = -a1;
			TT(0, 2) = b1;
			TT(4, 5) = -a2;
			TT(3, 5) = b2;

			T = Matrix<double>(6, 6);
			Matrix<double> delta1 = (node2.pos - node1.pos) / L;
			for (int j = 0; j < 2; j++) {
				T(0 + 3 * j, 0 + 3 * j) = delta1(0);
				T(0 + 3 * j, 1 + 3 * j) = delta1(1);
				T(1 + 3 * j, 0 + 3 * j) = -delta1(1);
				T(1 + 3 * j, 1 + 3 * j) = delta1(0);
				T(2 + 3 * j, 2 + 3 * j) = 1;
			}

			Kee = TT.trans() * TK * TT;
		}
		else if (nodeDOF == 6) { // 三维
			fscanf(fp, "%d%d%d%d%d%d", &nodeNum1, &nodeNum2, &nodeNum3, &material, &section, &offset);

			E = materials[material - 1].E;
			poisson = materials[material - 1].poisson;
			G = E / (2 * (1 + poisson));

			A = sections[section - 1].A;
			Ix = sections[section - 1].Ix;
			Iy = sections[section - 1].Iy;
			Iz = sections[section - 1].Iz;
			Az = sections[section - 1].Az;
			Ay = sections[section - 1].Ay;
			K = sections[section - 1].K;

			a1 = offsets[offset - 1].a1;
			b1 = offsets[offset - 1].b1;
			c1 = offsets[offset - 1].c1;
			a2 = offsets[offset - 1].a2;
			b2 = offsets[offset - 1].b2;
			c2 = offsets[offset - 1].c2;

			Node node1 = nodes[nodeNum1 - 1];
			Node node2 = nodes[nodeNum2 - 1];
			Node node3 = nodes[nodeNum3 - 1];
			L = calLength(node1, node2) + a1 - a2;

			by = 12 * K * E * Iz / G / Ay / L / L;
			bz = 12 * K * E * Iy / G / Az / L / L;

			TK = Matrix<double>(12, 12);
			TK(0, 0) = TK(6, 6) = E * A / L;
			TK(1, 1) = TK(7, 7) = 12 * E * Iz / (1 + by) / L / L / L;
			TK(2, 2) = TK(8, 8) = 12 * E * Iy / (1 + bz) / L / L / L;
			TK(3, 3) = TK(9, 9) = G * Ix / L;
			TK(4, 4) = TK(10, 10) = (4 + bz) * E * Iy / (1 + bz) / L;
			TK(5, 5) = TK(11, 11) = (4 + by) * E * Iz / (1 + by) / L;
			TK(4, 2) = TK(2, 4) = TK(10, 2) = TK(2, 10) = -6 * E * Iy / (1 + bz) / L / L;
			TK(5, 1) = TK(1, 5) = TK(11, 1) = TK(1, 11) = 6 * E * Iz / (1 + by) / L / L;
			TK(6, 0) = TK(0, 6) = -E * A / L;
			TK(7, 1) = TK(1, 7) = -12 * E * Iz / (1 + by) / L / L / L;
			TK(7, 5) = TK(5, 7) = TK(11, 7) = TK(7, 11) = -6 * E * Iz / (1 + by) / L / L;
			TK(8, 2) = TK(2, 8) = -12 * E * Iy / (1 + bz) / L / L / L;
			TK(8, 4) = TK(4, 8) = TK(10, 8) = TK(8, 10) = 6 * E * Iy / (1 + bz) / L / L;
			TK(9, 3) = TK(3, 9) = -G * Ix / L;
			TK(10, 4) = TK(4, 10) = (2 - bz) * E * Iy / (1 + bz) / L;
			TK(11, 5) = TK(5, 11) = (2 - by) * E * Iz / (1 + by) / L;

			TT = Matrix<double>(12, 12);
			for (int j = 0; j < 12; j++)
				TT(j, j) = 1;
			TT(2, 4) = a1; TT(1, 5) = -a1;
			TT(0, 5) = b1; TT(2, 3) = -b1;
			TT(1, 3) = c1; TT(0, 4) = -c1;
			TT(8, 10) = a2; TT(7, 11) = -a2;
			TT(6, 11) = b2; TT(8, 9) = -b2;
			TT(7, 9) = c2; TT(6, 10) = -c2;

			T = Matrix<double>(12, 12);
			Matrix<double> delta1 = (node2.pos - node1.pos) / L;
			Matrix<double> delta2 = (node3.pos - node1.pos) / calLength(node1, node3);
			for (int j = 0; j < 4; j++) {
				T(0 + 3 * j, 0 + 3 * j) = delta1(0);
				T(0 + 3 * j, 1 + 3 * j) = delta1(1);
				T(0 + 3 * j, 2 + 3 * j) = delta1(2);
				T(1 + 3 * j, 0 + 3 * j) = delta2(0);
				T(1 + 3 * j, 1 + 3 * j) = delta2(1);
				T(1 + 3 * j, 2 + 3 * j) = delta2(2);
				T(2 + 3 * j, 0 + 3 * j) = delta1(1) * delta2(2) - delta1(2) * delta2(1);
				T(2 + 3 * j, 1 + 3 * j) = delta1(2) * delta2(0) - delta1(0) * delta2(2);
				T(2 + 3 * j, 2 + 3 * j) = delta1(0) * delta2(1) - delta1(1) * delta2(0);
			}

			Kee = TT.trans() * TK * TT;
		}

		elements[i].nodeNum1 = nodeNum1;
		elements[i].nodeNum2 = nodeNum2;
		if (nodeDOF == 6)
			elements[i].nodeNum3 = nodeNum3;
		elements[i].materialNum = material;
		elements[i].sectionNum = section;
		elements[i].offsetNum = offset;
		elements[i].L = L;
		elements[i].Kee = Kee;
		elements[i].T = T;
		elements[i].TT = TT;
	}

	int constraintNum;
	fscanf(fp, "%d", &constraintNum);
	constraints = Array<Constraint>(constraintNum);
	for (int i = 0; i < constraintNum; i++) {
		fscanf(fp, "%d", &(constraints[i].dDirection));
	}

	int loadNum;
	fscanf(fp, "%d", &loadNum);
	loads = Array<Load>(loadNum);
	for (int i = 0; i < loadNum; i++) {
		fscanf(fp, "%d%lf", &(loads[i].dDirection), &(loads[i].force));
	}
}

void printParameters(int elementType, int nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads) {
	printf("---------------------------------------输入参数---------------------------------------\n");

	printf("单元种类：");
	if (elementType == 1)
		printf("杆单元\n");
	else if (elementType == 2)
		printf("梁单元\n");
	printf("\n");


	printf("节点自由度数：%d\n", nodeDOF);
	printf("\n");

	printf("材料总数：%d\n", materials.size());
	for (int i = 0; i < materials.size(); i++) {
		printf("第%d种材料：E=%lf", i + 1, materials[i].E);
		if (elementType == 2)
			printf(", 泊松比=%lf", materials[i].poisson);
		printf("\n");
	}
	printf("\n");

	printf("截面总数：%d\n", sections.size());
	for (int i = 0; i < sections.size(); i++) {
		printf("第%d种截面：A=%lf", i + 1, sections[i].A);
		if (elementType == 2) {
			if(nodeDOF == 6)
				printf(", Ix=%lf, Iy=%lf, Iz=%lf, Az=%lf, Ay=%lf, K=%lf", sections[i].Ix, sections[i].Iy, sections[i].Iz, sections[i].Az, sections[i].Ay, sections[i].K);
			else if(nodeDOF == 3)
				printf(", Iz=%lf, Ay=%lf, K=%lf", sections[i].Iz, sections[i].Ay, sections[i].K);
		}
		printf("\n");
	}
	printf("\n");

	if (elementType == 2) {
		printf("偏心总数：%d\n", offsets.size());
		for (int i = 0; i < offsets.size(); i++) {
			printf("第%d种偏心：", i + 1);
			if (nodeDOF == 6)
				printf("a1=%lf, b1=%lf, c1=%lf, a2=%lf, b2=%lf, c2=%lf", offsets[i].a1, offsets[i].b1, offsets[i].c1, offsets[i].a2, offsets[i].b2, offsets[i].c2);
			else if(nodeDOF == 3)
				printf("a1=%lf, b1=%lf, a2=%lf, b2=%lf", offsets[i].a1, offsets[i].b1, offsets[i].a2, offsets[i].b2);
		}
		printf("\n");
	}
	printf("\n");

	printf("节点总数：%d\n", nodes.size());
	printf("节点坐标：\n");
	for (int i = 0; i < nodes.size(); i++) {
		printf("第%d个节点：(", i + 1);
		if (elementType == 1) {
			for (int j = 0; j < nodeDOF; j++) {
				if (j != 0)
					printf(",");
				printf("%lf", nodes[i].pos(j));
			}
		}
		else if (elementType == 2) {
			for (int j = 0; j < ((nodeDOF == 3) ? 2 : 3); j++) {
				if (j != 0)
					printf(",");
				printf("%lf", nodes[i].pos(j));
			}
		}
		printf(")\n");
	}
	printf("\n");

	printf("单元总数：%d\n", elements.size());
	printf("单元信息：\n");
	for (int i = 0; i < elements.size(); i++) {
		printf("第%d个单元：%d号节点与%d号节点，第%d种材料，第%d种截面，L=%lf", i + 1, elements[i].nodeNum1, elements[i].nodeNum2, elements[i].materialNum, elements[i].sectionNum, elements[i].L);
		if (elementType == 2)
			printf("，第%d种偏心", elements[i].offsetNum);
		printf("\n");
	}
	printf("\n");

	printf("受约束的自由度总数：%d\n", constraints.size());
	for (int i = 0; i < constraints.size(); i++) {
		printf("约束的第%d个位移号：%d\n", i + 1, constraints[i].dDirection);
	}
	printf("\n");

	printf("外载荷总数：%d\n", loads.size());
	for (int i = 0; i < loads.size(); i++) {
		printf("第%d个位移方向上的外载荷：F=%lf\n", loads[i].dDirection, loads[i].force);
	}
	printf("\n");
}

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

void solve(int nodeDOF, Array<Node> &nodes, MatrixIn1D &K, Array<Load> &loads) {
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

void calRods(int elementType, Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements) {
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