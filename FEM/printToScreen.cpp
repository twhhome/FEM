#include "printToScreen.h"

void printParameters(int elementType, int nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads) {
	printf("---------------------------------------�������---------------------------------------\n");

	printf("��Ԫ���ࣺ");
	if (elementType == 1)
		printf("�˵�Ԫ\n");
	else if (elementType == 2)
		printf("����Ԫ\n");
	printf("\n");


	printf("�ڵ����ɶ�����%d\n", nodeDOF);
	printf("\n");

	printf("����������%d\n", materials.size());
	for (int i = 0; i < materials.size(); i++) {
		printf("��%d�ֲ��ϣ�E=%lf", i + 1, materials[i].E);
		if (elementType == 2)
			printf(", ���ɱ�=%lf", materials[i].poisson);
		printf("\n");
	}
	printf("\n");

	printf("����������%d\n", sections.size());
	for (int i = 0; i < sections.size(); i++) {
		printf("��%d�ֽ��棺A=%lf", i + 1, sections[i].A);
		if (elementType == 2) {
			if (nodeDOF == 6)
				printf(", Ix=%lf, Iy=%lf, Iz=%lf, Az=%lf, Ay=%lf, K=%lf", sections[i].Ix, sections[i].Iy, sections[i].Iz, sections[i].Az, sections[i].Ay, sections[i].K);
			else if (nodeDOF == 3)
				printf(", Iz=%lf, Ay=%lf, K=%lf", sections[i].Iz, sections[i].Ay, sections[i].K);
		}
		printf("\n");
	}
	printf("\n");

	if (elementType == 2) {
		printf("ƫ��������%d\n", offsets.size());
		for (int i = 0; i < offsets.size(); i++) {
			printf("��%d��ƫ�ģ�", i + 1);
			if (nodeDOF == 6)
				printf("a1=%lf, b1=%lf, c1=%lf, a2=%lf, b2=%lf, c2=%lf", offsets[i].a1, offsets[i].b1, offsets[i].c1, offsets[i].a2, offsets[i].b2, offsets[i].c2);
			else if (nodeDOF == 3)
				printf("a1=%lf, b1=%lf, a2=%lf, b2=%lf", offsets[i].a1, offsets[i].b1, offsets[i].a2, offsets[i].b2);
		}
		printf("\n");
	}
	printf("\n");

	printf("�ڵ�������%d\n", nodes.size());
	printf("�ڵ����꣺\n");
	for (int i = 0; i < nodes.size(); i++) {
		printf("��%d���ڵ㣺(", i + 1);
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

	printf("��Ԫ������%d\n", elements.size());
	printf("��Ԫ��Ϣ��\n");
	for (int i = 0; i < elements.size(); i++) {
		printf("��%d����Ԫ��%d�Žڵ���%d�Žڵ㣬��%d�ֲ��ϣ���%d�ֽ��棬L=%lf", i + 1, elements[i].nodeNum1, elements[i].nodeNum2, elements[i].materialNum, elements[i].sectionNum, elements[i].L);
		if (elementType == 2)
			printf("����%d��ƫ��", elements[i].offsetNum);
		printf("\n");
	}
	printf("\n");

	printf("��Լ�������ɶ�������%d\n", constraints.size());
	for (int i = 0; i < constraints.size(); i++) {
		printf("Լ���ĵ�%d��λ�ƺţ�%d\n", i + 1, constraints[i].dDirection);
	}
	printf("\n");

	printf("���غ�������%d\n", loads.size());
	for (int i = 0; i < loads.size(); i++) {
		printf("��%d��λ�Ʒ����ϵ����غɣ�F=%lf\n", loads[i].dDirection, loads[i].force);
	}
	printf("\n");
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

void printFinishMessage() {
	printf("---------------------------------------finished---------------------------------------\n");
}