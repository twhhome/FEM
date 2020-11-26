#include "FEM.h"

int main(int argc, char **argv) {
	char filename[80];
	if (argc == 1) {
		printf("请输入文件名：");
		scanf("%s", filename);
	}
	else if (argc == 2) {
		strcpy(filename, argv[1]);
	}
	else
		exit(0);

	FEM fem(filename);
	fem.solve();

	return 0;
}