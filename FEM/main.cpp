#include "FEMCore.h"

void readFile(const char* filename, int &nodeDOF, int &elementType, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads) {
	FILE* fp = fopen(filename, "r");

	fscanf(fp, "%d", &elementType);

	if (elementType == 1) // Rod
		readRodsFromFile(fp, nodeDOF, materials, sections, nodes, elements, constraints, loads);
	else if(elementType == 2) // Beam
		readBeamsFromFile(fp, nodeDOF, materials, sections, offsets, nodes, elements, constraints, loads);

	fclose(fp);

	printParameters(elementType, nodeDOF, materials, sections, offsets, nodes, elements, constraints, loads);
}

void printSolution(int elementType, int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints) {
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
	for (int i = 0; i < elements.size(); i++) {
		printf("��%d����Ԫ��\n", i + 1);
		printf("��������ϵ�µĽڵ�λ�ƣ�(");
		for (int j = 0; j < elements[i].de.n; j++) {
			if (j % nodeDOF == 0 && j != 0)
				printf("),(");
			else if (j != 0)
				printf(",");
			printf("%e", elements[i].de(j));
		}
		printf(")\n");

		printf("�ֲ�����ϵ�µĽڵ�λ�ƣ�(");
		for (int j = 0; j < elements[i].dee.n; j++) {
			if (j == elements[i].dee.n / 2)
				printf("),(");
			else if (j != 0)
				printf(",");
			printf("%e", elements[i].dee(j));
		}
		printf(")\n");

		printf("��������ϵ�µĽڵ�����(");
		for (int j = 0; j < elements[i].fe.n; j++) {
			if (j % nodeDOF == 0 && j != 0)
				printf("),(");
			else if (j != 0)
				printf(",");
			printf("%e", elements[i].fe(j));
		}
		printf(")\n");
		
		printf("�ֲ�����ϵ�µĽڵ�����(");
		for (int j = 0; j < elements[i].fee.n; j++) {
			if (j == elements[i].fee.n / 2)
				printf("),(");
			else if (j != 0)
				printf(",");
			printf("%e", elements[i].fee(j));
		}
		printf(")\n");

		if (elementType == 1) {
			printf("��Ԫ������%lf\n", elements[i].IF);
			printf("��ԪӦ����%lf\n", elements[i].stress);
		}
	}
	printf("\n");

	printf("Լ��������\n");
	for (int i = 0; i < constraints.size(); i++) {
		printf("λ�Ʒ���%d��Լ��������%lf\n", constraints[i].dDirection, constraints[i].force);
	}
}

int main() {
	int elementType;
	int nodeDOF;
	Array<Material> materials;
	Array<Section> sections;
	Array<Offset> offsets;
	Array<Node> nodes;
	Array<Element> elements;
	Array<Constraint> constraints;
	Array<Load> loads;
	
	readFile("5.txt", nodeDOF, elementType, materials, sections, offsets, nodes, elements, constraints, loads);
	
	MatrixIn1D K;
	calStiffnessMatrix(nodeDOF, nodes, elements, K);
	processConstraints(constraints, K);

	solve(nodeDOF, nodes, K, loads);
	
	calRods(elementType, sections, nodeDOF, nodes, elements);
	calConstraintForce(nodeDOF, nodes, elements, constraints);
	
	printSolution(elementType, nodeDOF, nodes, elements, constraints);
	
	printf("---------------------------------------finished---------------------------------------\n");
	return 0;
}