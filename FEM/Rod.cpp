#define _CRT_SECURE_NO_WARNINGS

#include "FEMCore.h"

double calLength(Node &node1, Node &node2) {
	Matrix<double> delta = node1.pos - node2.pos;
	return sqrt((delta.trans() * delta)(0));
}

void readFile(const char* filename, Array<Material> &materials, Array<Section> &sections, int &nodeDOF, Array<Node> &nodes, Array<Rod> &rods, Array<Constraint> &constraints, Array<double> &loads) {
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

	int rodNum;
	fscanf(fp, "%d", &rodNum);
	rods = Array<Rod>(rodNum);
	int nodeNum1, nodeNum2, material, section;
	double E, A, L;
	Matrix<double> TK(2, 2), T(2, 2 * nodeDOF);
	for (int i = 0; i < rodNum; i++) {
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

		rods[i].nodeNum1 = nodeNum1;
		rods[i].nodeNum2 = nodeNum2;
		rods[i].materialNum = material;
		rods[i].sectionNum = section;
		rods[i].L = L;
		rods[i].TK = TK;
		rods[i].T = T;
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

void printParameters(Array<Material> &materials, Array<Section> &sections, int nodeDOF, Array<Node> &nodes, Array<Rod> &rods, Array<Constraint> &constraints, Array<double> &loads) {
	printf("---------------------------------------�������---------------------------------------\n");

	printf("����������%d\n", materials.size());
	for (int i = 0; i < materials.size(); i++) {
		printf("��%d�ֲ��ϣ�E=%lf, ���ɱ�=%lf\n", i + 1, materials[i].E, materials[i].poisson);
	}
	printf("\n");

	printf("����������%d\n", sections.size());
	for (int i = 0; i < sections.size(); i++) {
		printf("��%d�ֽ��棺A=%lf\n", i + 1, sections[i].A);
	}
	printf("\n");

	printf("�ڵ�������%d\n", nodes.size());
	printf("�ڵ����ɶ�����%d\n", nodeDOF);
	printf("�ڵ����꣺\n");
	for (int i = 0; i < nodes.size(); i++) {
		printf("��%d���ڵ㣺(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				printf(",");
			printf("%lf", nodes[i].pos(j));
		}
		printf(")\n");
	}
	printf("\n");
	
	printf("��Ԫ������%d\n", rods.size());
	printf("��Ԫ��Ϣ��\n");
	for (int i = 0; i < rods.size(); i++) {
		printf("��%d����Ԫ��%d�Žڵ���%d�Žڵ㣬��%d�ֲ��ϣ���%d�ֽ��棬L=%lf\n", i + 1, rods[i].nodeNum1, rods[i].nodeNum2, rods[i].materialNum, rods[i].sectionNum, rods[i].L);
	}
	printf("\n");
	
	printf("��Լ�������ɶ�������%d\n", constraints.size());
	for (int i = 0; i < constraints.size(); i++) {
		printf("Լ���ĵ�%d��λ�ƺţ�%d\n", i + 1, constraints[i].dDirection);
	}
	printf("\n");

	for (int i = 0; i < loads.size(); i++) {
		printf("��%d��λ�Ʒ����ϵ����غɣ�%lf\n", i + 1, loads[i]);
	}
	printf("\n");
}

void printSolution(int nodeDOF, Array<Node> &nodes, Array<Rod> &rods, Array<Constraint> &constraints) {
	printf("\n");
	printf("---------------------------------------�����---------------------------------------\n");
	printf("�ڵ�λ�ƣ�\n");
	for (int i = 0; i < nodes.size(); i++) {
		printf("��%d���ڵ��λ�ƣ�(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				printf(",");
			printf("%e", nodes[i].displacement(j));
		}
		printf(")\n");
	}
	printf("\n");

	printf("��Ԫ��Ϣ��\n");
	for (int i = 0; i < rods.size(); i++) {
		printf("��%d����Ԫ��\n", i + 1);
		if(nodeDOF == 2)
			printf("��������ϵ�µĽڵ�λ�ƣ�(%e,%e),(%e,%e)\n", rods[i].de(0), rods[i].de(1), rods[i].de(2), rods[i].de(3));
		else if(nodeDOF == 3)
			printf("��������ϵ�µĽڵ�λ�ƣ�(%e,%e,%e),(%e,%e,%e)\n", rods[i].de(0), rods[i].de(1), rods[i].de(2), rods[i].de(3), rods[i].de(4), rods[i].de(5));
		printf("�ֲ�����ϵ�µĽڵ�λ�ƣ�%e,%e\n", rods[i].dee(0), rods[i].dee(1));
		if (nodeDOF == 2)
			printf("��������ϵ�µĽڵ�����(%lf,%lf),(%lf,%lf)\n", rods[i].fe(0), rods[i].fe(1), rods[i].fe(2), rods[i].fe(3));
		else if (nodeDOF == 3)
			printf("��������ϵ�µĽڵ�����(%lf,%lf,%lf),(%lf,%lf,%lf)\n", rods[i].fe(0), rods[i].fe(1), rods[i].fe(2), rods[i].fe(3), rods[i].fe(4), rods[i].fe(5));
		printf("�ֲ�����ϵ�µĽڵ�����%lf,%lf\n", rods[i].fee(0), rods[i].fee(1));
		printf("��Ԫ������%lf\n", rods[i].IF);
		printf("��ԪӦ����%lf\n", rods[i].stress);
	}
	printf("\n");

	printf("Լ��������\n");
	for (int i = 0; i < constraints.size(); i++) {
		printf("λ�Ʒ���%d��Լ��������%lf\n", constraints[i].dDirection, constraints[i].force);
	}
}

int main() {
	int nodeDOF;
	Array<Material> materials;
	Array<Section> sections;
	Array<Node> nodes;
	Array<Rod> rods;
	Array<Constraint> constraints;
	Array<double> loads;
	
	readFile("3.txt", materials, sections, nodeDOF, nodes, rods, constraints, loads);
	printParameters(materials, sections, nodeDOF, nodes, rods, constraints, loads);
	
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, rods, K);
	processConstraints(constraints, K);

	solve(nodeDOF, nodes, K, loads);
	
	calRods(sections, nodeDOF, nodes, rods);
	calConstraintForce(nodeDOF, nodes, rods, constraints);
	
	printSolution(nodeDOF, nodes, rods, constraints);
	
	printf("---------------------------------------finished---------------------------------------\n");
	return 0;
}