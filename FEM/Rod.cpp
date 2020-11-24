#define _CRT_SECURE_NO_WARNINGS

#include "FEMCore.h"

double calLength(Node &node1, Node &node2) {
	Matrix<double> delta = node1.pos - node2.pos;
	return sqrt((delta.trans() * delta)(0));
}

void readFile(const char* filename, int &nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<Constraint> &constraints, Matrix<double> &loads) {
	FILE* fp = fopen(filename, "r");

	int nodeNum;
	fscanf(fp, "%d", &nodeNum);
	fscanf(fp, "%d", &nodeDOF);
	nodes = Matrix<Node>(nodeNum, 1);
	for (int i = 0; i < nodeNum; i++) {
		nodes(i).pos = Matrix<double>(nodeDOF, 1);
		nodes(i).displacement = Matrix<double>(nodeDOF, 1);
		for (int j = 0; j < nodeDOF; j++) {
			fscanf(fp, "%lf", &(nodes(i).pos(j)));
		}
	}

	int rodNum;
	fscanf(fp, "%d", &rodNum);
	rods = Matrix<Rod>(rodNum, 1);
	int nodeNum1, nodeNum2;
	double E, A, L;
	Matrix<double> TK(2, 2), T(2, 2 * nodeDOF);
	for (int i = 0; i < rodNum; i++) {
		fscanf(fp, "%d%d%lf%lf", &nodeNum1, &nodeNum2, &E, &A);

		Node node1 = nodes(nodeNum1 - 1);
		Node node2 = nodes(nodeNum2 - 1);
		L = calLength(node1, node2);
		
		TK(0, 0) = TK(1, 1) = E * A / L;
		TK(0, 1) = TK(1, 0) = -E * A / L;

		Matrix<double> delta = node2.pos - node1.pos;
		delta = delta / L;
		for (int j = 0; j < nodeDOF; j++) {
			T(0, j) = T(1, j + nodeDOF) = delta(j);
		}

		rods(i).nodeNum1 = nodeNum1;
		rods(i).nodeNum2 = nodeNum2;
		rods(i).E = E;
		rods(i).A = A;
		rods(i).L = L;
		rods(i).TK = TK;
		rods(i).T = T;
	}

	int constraintNum;
	fscanf(fp, "%d", &constraintNum);
	constraints = Matrix<Constraint>(constraintNum, 1);
	for (int i = 0; i < constraintNum; i++) {
		fscanf(fp, "%d", &(constraints(i).dDirection));
	}

	loads = Matrix<double>(nodeDOF * nodeNum, 1);
	for (int i = 0; i < loads.n; i++) {
		fscanf(fp, "%lf", &(loads(i)));
	}

	fclose(fp);
}

void printParameters(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<Constraint> &constraints, Matrix<double> &loads) {
	printf("---------------------------------------输入参数---------------------------------------\n");

	printf("总节点数：%d\n", nodes.n);
	printf("节点自由度数：%d\n", nodeDOF);
	printf("节点坐标：\n");
	for (int i = 0; i < nodes.n; i++) {
		printf("第%d个节点：(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				printf(",");
			printf("%lf", nodes(i).pos(j));
		}
		printf(")\n");
	}
	printf("\n");
	
	printf("单元总数：%d\n", rods.n);
	printf("单元信息：\n");
	for (int i = 0; i < rods.n; i++) {
		printf("第%d个单元：%d号节点与%d号节点，E=%lf，A=%lf，L=%lf\n", i + 1, rods(i).nodeNum1, rods(i).nodeNum2, rods(i).E, rods(i).A, rods(i).L);
	}
	printf("\n");
	
	printf("受约束的自由度总数：%d\n", constraints.n);
	for (int i = 0; i < constraints.n; i++) {
		printf("约束的第%d个位移号：%d\n", i + 1, constraints(i).dDirection);
	}
	printf("\n");

	for (int i = 0; i < loads.n; i++) {
		printf("第%d个位移方向上的外载荷：%lf\n", i + 1, loads(i));
	}
	printf("\n");
}

void printSolution(int nodeDOF, Matrix<Node> &nodes, Matrix<Rod> &rods, Matrix<Constraint> &constraints) {
	printf("\n");
	printf("---------------------------------------求解结果---------------------------------------\n");
	printf("节点位移：\n");
	for (int i = 0; i < nodes.n; i++) {
		printf("第%d个节点的位移：(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				printf(",");
			printf("%lf", nodes(i).displacement(j));
		}
		printf(")\n");
	}
	printf("\n");

	printf("单元信息：\n");
	for (int i = 0; i < rods.n; i++) {
		printf("第%d个单元：\n", i + 1);
		if(nodeDOF == 2)
			printf("整体坐标系下的节点位移：(%lf,%lf),(%lf,%lf)\n", rods(i).de(0), rods(i).de(1), rods(i).de(2), rods(i).de(3));
		else if(nodeDOF == 3)
			printf("整体坐标系下的节点位移：(%lf,%lf,%lf),(%lf,%lf,%lf)\n", rods(i).de(0), rods(i).de(1), rods(i).de(2), rods(i).de(3), rods(i).de(4), rods(i).de(5));
		printf("局部坐标系下的节点位移：%lf,%lf\n", rods(i).dee(0), rods(i).dee(1));
		if (nodeDOF == 2)
			printf("整体坐标系下的节点力：(%lf,%lf),(%lf,%lf)\n", rods(i).fe(0), rods(i).fe(1), rods(i).fe(2), rods(i).fe(3));
		else if (nodeDOF == 3)
			printf("整体坐标系下的节点力：(%lf,%lf,%lf),(%lf,%lf,%lf)\n", rods(i).fe(0), rods(i).fe(1), rods(i).fe(2), rods(i).fe(3), rods(i).fe(4), rods(i).fe(5));
		printf("局部坐标系下的节点力：%lf,%lf\n", rods(i).fee(0), rods(i).fee(1));
		printf("单元内力：%lf\n", rods(i).IF);
		printf("单元应力：%lf\n", rods(i).stress);
	}
	printf("\n");

	printf("约束反力：\n");
	for (int i = 0; i < constraints.n; i++) {
		printf("位移方向%d的约束反力：%lf\n", constraints(i).dDirection, constraints(i).force);
	}
}

int main() {
	int nodeDOF;
	Matrix<Node> nodes;
	Matrix<Rod> rods;
	Matrix<Constraint> constraints;
	Matrix<double> loads;
	
	readFile("4.txt", nodeDOF, nodes, rods, constraints, loads);
	printParameters(nodeDOF, nodes, rods, constraints, loads);
	
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, rods, K);
	processConstraints(constraints, K);

	solve(nodeDOF, nodes, K, loads);
	
	calRods(nodeDOF, nodes, rods);
	calConstraintForce(nodeDOF, nodes, rods, constraints);
	
	printSolution(nodeDOF, nodes, rods, constraints);
	
	printf("---------------------------------------finished---------------------------------------\n");
	return 0;
}