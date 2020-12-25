#include "FEM.h"

int main(int argc, char **argv) {
	char inputFileName[80] = { 0 }, outputFileName[80] = { 0 };
	if (argc == 1) {
		while (strlen(inputFileName) == 0) {
			printf("����������ļ�����");
			gets_s(inputFileName);
		}
		printf("����������ļ��������»س����������Ļ����");
		gets_s(outputFileName);
	}
	else if (argc == 2) {
		strcpy(inputFileName, argv[1]);
	}
	else if (argc == 3) {
		strcpy(inputFileName, argv[1]);
		strcpy(outputFileName, argv[2]);
	}
	else
		exit(0);

	FEM fem(inputFileName, outputFileName);
	fem.solve();

	return 0;
}