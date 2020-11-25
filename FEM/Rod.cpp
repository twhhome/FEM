#define _CRT_SECURE_NO_WARNINGS

#include "FEMCore.h"

double calLength(Node &node1, Node &node2) {
	Matrix<double> delta = node1.pos - node2.pos;
	return sqrt((delta.trans() * delta)(0));
}

void readFile(const char* filename, Array<Material> &materials, Array<Section> &sections, int &nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<double> &loads) {
	FILE* fp = fopen(filename, "r");

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
		double A;
		fscanf(fp, "%lf", &A);
		sections[i].A = A;
	}

	int nodeNum;
	fscanf(fp, "%d", &nodeNum);
	fscanf(fp, "%d", &nodeDOF);
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
		elements[i].TK = TK;
		elements[i].T = T;
	}

	int constraintNum;
	fscanf(fp, "%d", &constraintNum);
	constraints = Array<Constraint>(constraintNum);
	for (int i = 0; i < constraintNum; i++) {
		fscanf(fp, "%d", &(constraints[i].dDirection));
	}

	loads = Array<double>(nodeDOF * nodeNum);
	for (int i = 0; i < loads.size(); i++) {
		fscanf(fp, "%lf", &(loads[i]));
	}

	fclose(fp);
}

void printParameters(Array<Material> &materials, Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<double> &loads) {
	printf("---------------------------------------输入参数---------------------------------------\n");

	printf("材料总数：%d\n", materials.size());
	for (int i = 0; i < materials.size(); i++) {
		printf("第%d种材料：E=%lf, 泊松比=%lf\n", i + 1, materials[i].E, materials[i].poisson);
	}
	printf("\n");

	printf("截面总数：%d\n", sections.size());
	for (int i = 0; i < sections.size(); i++) {
		printf("第%d种截面：A=%lf\n", i + 1, sections[i].A);
	}
	printf("\n");

	printf("节点总数：%d\n", nodes.size());
	printf("节点自由度数：%d\n", nodeDOF);
	printf("节点坐标：\n");
	for (int i = 0; i < nodes.size(); i++) {
		printf("第%d个节点：(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				printf(",");
			printf("%lf", nodes[i].pos(j));
		}
		printf(")\n");
	}
	printf("\n");
	
	printf("单元总数：%d\n", elements.size());
	printf("单元信息：\n");
	for (int i = 0; i < elements.size(); i++) {
		printf("第%d个单元：%d号节点与%d号节点，第%d种材料，第%d种截面，L=%lf\n", i + 1, elements[i].nodeNum1, elements[i].nodeNum2, elements[i].materialNum, elements[i].sectionNum, elements[i].L);
	}
	printf("\n");
	
	printf("受约束的自由度总数：%d\n", constraints.size());
	for (int i = 0; i < constraints.size(); i++) {
		printf("约束的第%d个位移号：%d\n", i + 1, constraints[i].dDirection);
	}
	printf("\n");

	for (int i = 0; i < loads.size(); i++) {
		printf("第%d个位移方向上的外载荷：%lf\n", i + 1, loads[i]);
	}
	printf("\n");
}

void printSolution(int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints) {
	printf("\n");
	printf("---------------------------------------求解结果---------------------------------------\n");
	printf("节点位移：\n");
	for (int i = 0; i < nodes.size(); i++) {
		printf("第%d个节点的位移：(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				printf(",");
			printf("%e", nodes[i].displacement(j));
		}
		printf(")\n");
	}
	printf("\n");

	printf("单元信息：\n");
	for (int i = 0; i < elements.size(); i++) {
		printf("第%d个单元：\n", i + 1);
		if(nodeDOF == 2)
			printf("整体坐标系下的节点位移：(%e,%e),(%e,%e)\n", elements[i].de(0), elements[i].de(1), elements[i].de(2), elements[i].de(3));
		else if(nodeDOF == 3)
			printf("整体坐标系下的节点位移：(%e,%e,%e),(%e,%e,%e)\n", elements[i].de(0), elements[i].de(1), elements[i].de(2), elements[i].de(3), elements[i].de(4), elements[i].de(5));
		printf("局部坐标系下的节点位移：%e,%e\n", elements[i].dee(0), elements[i].dee(1));
		if (nodeDOF == 2)
			printf("整体坐标系下的节点力：(%lf,%lf),(%lf,%lf)\n", elements[i].fe(0), elements[i].fe(1), elements[i].fe(2), elements[i].fe(3));
		else if (nodeDOF == 3)
			printf("整体坐标系下的节点力：(%lf,%lf,%lf),(%lf,%lf,%lf)\n", elements[i].fe(0), elements[i].fe(1), elements[i].fe(2), elements[i].fe(3), elements[i].fe(4), elements[i].fe(5));
		printf("局部坐标系下的节点力：%lf,%lf\n", elements[i].fee(0), elements[i].fee(1));
		printf("单元内力：%lf\n", elements[i].IF);
		printf("单元应力：%lf\n", elements[i].stress);
	}
	printf("\n");

	printf("约束反力：\n");
	for (int i = 0; i < constraints.size(); i++) {
		printf("位移方向%d的约束反力：%lf\n", constraints[i].dDirection, constraints[i].force);
	}
}

int main() {
	int nodeDOF;
	Array<Material> materials;
	Array<Section> sections;
	Array<Node> nodes;
	Array<Element> elements;
	Array<Constraint> constraints;
	Array<double> loads;
	
	readFile("3.txt", materials, sections, nodeDOF, nodes, elements, constraints, loads);
	printParameters(materials, sections, nodeDOF, nodes, elements, constraints, loads);
	
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, elements, K);
	processConstraints(constraints, K);

	solve(nodeDOF, nodes, K, loads);
	
	calRods(sections, nodeDOF, nodes, elements);
	calConstraintForce(nodeDOF, nodes, elements, constraints);
	
	printSolution(nodeDOF, nodes, elements, constraints);
	
	printf("---------------------------------------finished---------------------------------------\n");
	return 0;
}