#include "FEM.h"

int main(int argc, char **argv) {
	char inputFileName[80] = { 0 }, outputFileName[80] = { 0 };
	if (argc == 1) {
		while (strlen(inputFileName) == 0) {
			printf("请输入读入文件名：");
			gets_s(inputFileName);
		}
		printf("请输入输出文件名（按下回车可输出到屏幕）：");
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