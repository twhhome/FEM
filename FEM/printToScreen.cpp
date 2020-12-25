#include "printToScreen.h"

void printParameters(int elementType, int nodeDOF, Array<Material> &materials, Array<Section> &sections, Array<Offset> &offsets, Array<Node> &nodes, Array<Node> &secondNodes, Array<Element> &elements, Array<Constraint> &constraints, Array<Load> &loads) {
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
			if (nodeDOF == 6)
				printf(", Ix=%lf, Iy=%lf, Iz=%lf, Az=%lf, Ay=%lf, K=%lf", sections[i].Ix, sections[i].Iy, sections[i].Iz, sections[i].Az, sections[i].Ay, sections[i].K);
			else if (nodeDOF == 3)
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
			else if (nodeDOF == 3)
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

	if (elementType == 2 && nodeDOF == 6) {
		printf("辅助节点总数：%d\n", secondNodes.size());
		printf("辅助节点坐标：\n");
		for (int i = 0; i < secondNodes.size(); i++) {
			printf("第%d个辅助节点：(", i + 1);
			for (int j = 0; j < 3; j++) {
				if (j != 0)
					printf(",");
				printf("%lf", secondNodes[i].pos(j));
			}
			printf(")\n");
		}
		printf("\n");
	}

	printf("单元总数：%d\n", elements.size());
	printf("单元信息：\n");
	for (int i = 0; i < elements.size(); i++) {
		printf("第%d个单元：%d号节点与%d号节点，第%d种材料，第%d种截面，L=%lf", i + 1, elements[i].nodeNum1, elements[i].nodeNum2, elements[i].materialNum, elements[i].sectionNum, elements[i].L);
		if (elementType == 2)
			printf("，第%d种偏心", elements[i].offsetNum);
		if(elementType == 2 && nodeDOF == 6)
			printf("，%d号辅助节点", elements[i].nodeNum3);
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

void outputSolution(const char* filename, int elementType, int nodeDOF, Array<Node> &nodes, Array<Element> &elements, Array<Constraint> &constraints) {
	FILE *fp = stdout;
	if (strlen(filename) != 0) {
		fp = fopen(filename, "w");
	}
	fprintf(fp, "\n");
	fprintf(fp, "---------------------------------------求解结果---------------------------------------\n");
	fprintf(fp, "节点位移：\n");
	for (int i = 0; i < nodes.size(); i++) {
		fprintf(fp, "第%d个节点的位移：(", i + 1);
		for (int j = 0; j < nodeDOF; j++) {
			if (j != 0)
				fprintf(fp, ",");
			fprintf(fp, "%e", nodes[i].displacement(j));
		}
		fprintf(fp, ")\n");
	}
	fprintf(fp, "\n");

	fprintf(fp, "单元信息：\n");
	for (int i = 0; i < elements.size(); i++) {
		fprintf(fp, "第%d个单元：\n", i + 1);
		fprintf(fp, "整体坐标系下的节点位移：(");
		for (int j = 0; j < elements[i].de.n; j++) {
			if (j % nodeDOF == 0 && j != 0)
				fprintf(fp, "),(");
			else if (j != 0)
				fprintf(fp, ",");
			fprintf(fp, "%e", elements[i].de(j));
		}
		fprintf(fp, ")\n");

		fprintf(fp, "局部坐标系下的节点位移：(");
		for (int j = 0; j < elements[i].dee.n; j++) {
			if (j == elements[i].dee.n / 2)
				fprintf(fp, "),(");
			else if (j != 0)
				fprintf(fp, ",");
			fprintf(fp, "%e", elements[i].dee(j));
		}
		fprintf(fp, ")\n");

		fprintf(fp, "整体坐标系下的节点力：(");
		for (int j = 0; j < elements[i].fe.n; j++) {
			if (j % nodeDOF == 0 && j != 0)
				fprintf(fp, "),(");
			else if (j != 0)
				fprintf(fp, ",");
			fprintf(fp, "%e", elements[i].fe(j));
		}
		fprintf(fp, ")\n");

		fprintf(fp, "局部坐标系下的节点力：(");
		for (int j = 0; j < elements[i].fee.n; j++) {
			if (j == elements[i].fee.n / 2)
				fprintf(fp, "),(");
			else if (j != 0)
				fprintf(fp, ",");
			fprintf(fp, "%e", elements[i].fee(j));
		}
		fprintf(fp, ")\n");

		if (elementType == 1) {
			fprintf(fp, "单元内力：%lf\n", elements[i].IF);
			fprintf(fp, "单元应力：%lf\n", elements[i].stress);
		}
	}
	fprintf(fp, "\n");

	fprintf(fp, "约束反力：\n");
	for (int i = 0; i < constraints.size(); i++) {
		fprintf(fp, "位移方向%d的约束反力：%lf\n", constraints[i].dDirection, constraints[i].force);
	}

	if (strlen(filename) != 0)
		fclose(fp);
}

void printFinishMessage() {
	printf("---------------------------------------finished---------------------------------------\n");
}