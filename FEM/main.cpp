#include <stdio.h>
#include <math.h>

int ND; // 每个单元的节点总数
int NF; // 单个节点的自由度数
int NP; // 节点总数
int NE; // 单元总数
int NR; // 受约束的自由度总数
int NM, NMN; // 单元类别总数，单元的类别数
int N; // N=NP*NF
int NN; // 一维存储AK的总容量
double X[50], Y[50], Z[50]; // 各个节点的三维坐标
double X2, X1, Y2, Y1, Z2, Z1, b;
int ME[3][30]; // 每个单元的节点号
int NRR[30]; // 约束的位移号
int LD[50];
int NAE[100]; // 每个单元的类别
double AE[3][100]; // AE[1][i]存放第i种类型的单元的杨氏模量，AE[2][i]存放该种类型单元的横截面积
double P[100], P1[100]; // 节点载荷
double PP[100]; // 整体结构节点力
double A; // 单元横截面
double E; // 单元杨氏模量
double TK[3][3]; // 单元刚度矩阵
double T[3][7]; // 坐标转换矩阵
double TT[7][3]; // 坐标转换矩阵的转置
double AK[100]; // 整体刚度矩阵
double AKEE[7][7]; // 整体坐标系下的单元刚度阵
double s[7][3]; // 作矩阵乘法时的中间矩阵
int IS[7];
double L; // 杆单元的长度
double SG; // 单元应力
double d[100]; // 结构位移矩阵
double ue[7]; // 单元位移矩阵
double dee[50][3]; // 局部坐标系下的单元位移矩阵
double Fee[50][7]; // 局部坐标系下的单元节点力（即单元内力）
double F[50][100]; // 单元整体坐标系中的节点力
double l[100][100], y[100]; // 解方程用到的L，Y矩阵

// 数据输入
void scan() {

}

// 求每个杆的长度（i表示单元号）
void Length(int i) {
	X2 = X[ME[2][i]];
	X1 = X[ME[1][i]];
	Y2 = Y[ME[2][i]];
	Y1 = Y[ME[1][i]];
	Z2 = Z[ME[2][i]];
	Z1 = Z[ME[1][i]];
	L = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1) + (Z2 - Z1) * (Z2 - Z1));
}

// 杆单元的单元刚度阵（单元坐标系下）（i表示单元号）
void StiffnessMatrix_unit(int i) {
	Length(i);
	NMN = NAE[i];
	E = AE[1][NMN];
	A = AE[2][NMN];
	TK[1][1] = E * A / L;
	TK[1][2] = -E * A / L;
	TK[2][1] = -E * A / L;
	TK[2][2] = E * A / L;
}

// 坐标转换矩阵（i表示单元号）
void TransformMatrix(int i) {
	int m, n;
	Length(i);
	T[1][1] = (X2 - X1) / L;
	T[1][2] = (Y2 - Y1) / L;
	T[1][3] = (Z2 - Z1) / L;
	T[2][4] = (X2 - X1) / L;
	T[2][5] = (Y2 - Y1) / L;
	T[2][6] = (Z2 - Z1) / L;
	for (m = 1; m <= 2; m++)
		for (n = 1; n <= (NF * ND); n++)
			TT[n][m] = T[m][n];
}

//总体坐标系下的单元刚度阵（矩阵乘法）（i表示单元号）
void MultiplyMatrix(int i) {
	int j, m, n;
	double b;
	StiffnessMatrix_unit(i);
	TransformMatrix(i);
	for (n = 1; n <= (NF * ND); n++) {
		for (m = 1; m <= 2; m++) {
			b = 0.0;
			for (j = 1; j <= 2; j++)
				b += TT[n][j] * TK[j][m];
			s[n][m] = b;
		}
	}
	for (m = 1; m <= (NF * ND); m++) {
		for (n = 1; n <= (NF * ND); n++) {
			b = 0.0;
			for (j = 1; j <= 2; j++)
				b += s[m][j] * T[j][n];
			AKEE[m][n] = b;
		}
	}
}

// 形成LD数组
void FLd() {
	int k, i, j, L, IG, J, NN, N;
	LD[1] = 1;
	for (k = 1; k <= NP; k++) {
		IG = 100000;
		for (i = 1; i <= NE; i++) {
			for (j = 1; j <= ND; j++) {
				if (ME[j][i] != k)
					continue;
				for (L = 1; L <= ND; L++) {
					if (ME[L][i] >= IG)
						continue;
					IG = ME[L][i];
				}
			}
		}
		for (i = 1; i <= NF; i++) {
			J = NF * (k - 1) + i;
			if (J == 1)
				continue;
			LD[J] = LD[j - 1] + NF * (k - IG) + i;
		}
	}
	N = NP * NF;
	NN = LD[N];
}

// 形成IS数组
void FIS(int i) {
	IS[1] = (ME[1][i] - 1) * NF + 1;
	IS[2] = (ME[1][i] - 1) * NF + 2;
	IS[3] = (ME[1][i] - 1) * NF + 3;
	IS[4] = (ME[2][i] - 1) * NF + 1;
	IS[5] = (ME[2][i] - 1) * NF + 2;
	IS[6] = (ME[2][i] - 1) * NF + 3;
}

// 组集结构刚度阵
void StructMatrix() {
	int i, j, m, ISS, NI, NJ, IJ;
	FLd();
	for (m = 1; m <= NE; m++) {
		MultiplyMatrix(m);
		FIS(m);
		for (i = 1; i <= (NF * ND); i++) {
			for (j = 1; j <= (NF * ND); j++) {
				ISS = IS[i] - IS[j];
				if (ISS >= 0) {
					NI = IS[i];
					IJ = LD[NI] - (NI - IS[j]);
					AK[IJ] += AKEE[i][j];
				}
			}
		}
	}
	for (i = 1; i <= NR; i++) {
		NI = NRR[i];
		NJ = LD[NI];
		AK[NJ] = 1e25;
	}
}

void cholesky(int n, double a[100], double x[100]) {
	int i, j, k, ij, kj, ii, jj, ik, jk, kk, iig, ig, igp, jgp, mi, mj, mij;
	for (i = 1; i <= n; i++) {
		if (i != 1) {
			mi = i - (LD[i] - LD[i - 1]) + 1;
			if (mi != i) {
				iig = LD[i] - i;
				for (j = mi; j <= i - 1; j++) {
					if (j != mi) {
						mj = j - (LD[j] - LD[j - 1]) + 1;
						igp = LD[i] - (i - j);
						if (mj < mi)
							mij = mi;
						else
							mij = mj;
						jgp = LD[j] - j;
						if (mij <= j - 1) {
							for (k = mij; k <= j - 1; k++) {
								ik = iig + k;
								ij = igp + k;
								kk = LD[k];
								a[igp] -= a[ik] * a[kk] * a[jk];
							}
						}
					}
					if (j == mi)
						igp = LD[i - 1] + 1;
					ii = LD[j];
					a[igp] = a[igp] / a[ii];
					x[i] -= a[igp] * a[ii] * x[j];
				}
				ij = LD[i];
				for (k = mi; k <= i - 1; k++) {
					ii = iig + k;
					jj = LD[k];
					a[ij] -= a[ii] * a[ii] * a[ij];
				}
			}
		}
		ij = LD[i];
		x[i] = x[i] / a[ij];
	}
	for (i = n; i >= 2; i--) {
		mi = i - (LD[i] - LD[i - 1]) + 1;
		if (mi == i)
			continue;
		iig = LD[i] - i;
		for (k = mi; k <= i - 1; k++) {
			ij = iig + k;
			x[k] -= a[ij] * x[i];
		}
	}
}

int main() {
	FILE* fp;
	fp = fopen("1.txt", "w");
	fprintf(fp, "%d\n", 123);
	fclose(fp);
	return 0;
}
