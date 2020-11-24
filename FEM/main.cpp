#include <stdio.h>
#include <math.h>

int ND; // ÿ����Ԫ�Ľڵ�����
int NF; // �����ڵ�����ɶ���
int NP; // �ڵ�����
int NE; // ��Ԫ����
int NR; // ��Լ�������ɶ�����
int NM, NMN; // ��Ԫ�����������Ԫ�������
int N; // N=NP*NF
int NN; // һά�洢AK��������
double X[50], Y[50], Z[50]; // �����ڵ����ά����
double X2, X1, Y2, Y1, Z2, Z1, b;
int ME[3][30]; // ÿ����Ԫ�Ľڵ��
int NRR[30]; // Լ����λ�ƺ�
int LD[50];
int NAE[100]; // ÿ����Ԫ�����
double AE[3][100]; // AE[1][i]��ŵ�i�����͵ĵ�Ԫ������ģ����AE[2][i]��Ÿ������͵�Ԫ�ĺ�����
double P[100], P1[100]; // �ڵ��غ�
double PP[100]; // ����ṹ�ڵ���
double A; // ��Ԫ�����
double E; // ��Ԫ����ģ��
double TK[3][3]; // ��Ԫ�նȾ���
double T[3][7]; // ����ת������
double TT[7][3]; // ����ת�������ת��
double AK[100]; // ����նȾ���
double AKEE[7][7]; // ��������ϵ�µĵ�Ԫ�ն���
double s[7][3]; // ������˷�ʱ���м����
int IS[7];
double L; // �˵�Ԫ�ĳ���
double SG; // ��ԪӦ��
double d[100]; // �ṹλ�ƾ���
double ue[7]; // ��Ԫλ�ƾ���
double dee[50][3]; // �ֲ�����ϵ�µĵ�Ԫλ�ƾ���
double Fee[50][7]; // �ֲ�����ϵ�µĵ�Ԫ�ڵ���������Ԫ������
double F[50][100]; // ��Ԫ��������ϵ�еĽڵ���
double l[100][100], y[100]; // �ⷽ���õ���L��Y����

// ��������
void scan() {

}

// ��ÿ���˵ĳ��ȣ�i��ʾ��Ԫ�ţ�
void Length(int i) {
	X2 = X[ME[2][i]];
	X1 = X[ME[1][i]];
	Y2 = Y[ME[2][i]];
	Y1 = Y[ME[1][i]];
	Z2 = Z[ME[2][i]];
	Z1 = Z[ME[1][i]];
	L = sqrt((X2 - X1) * (X2 - X1) + (Y2 - Y1) * (Y2 - Y1) + (Z2 - Z1) * (Z2 - Z1));
}

// �˵�Ԫ�ĵ�Ԫ�ն��󣨵�Ԫ����ϵ�£���i��ʾ��Ԫ�ţ�
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

// ����ת������i��ʾ��Ԫ�ţ�
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

//��������ϵ�µĵ�Ԫ�ն��󣨾���˷�����i��ʾ��Ԫ�ţ�
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

// �γ�LD����
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

// �γ�IS����
void FIS(int i) {
	IS[1] = (ME[1][i] - 1) * NF + 1;
	IS[2] = (ME[1][i] - 1) * NF + 2;
	IS[3] = (ME[1][i] - 1) * NF + 3;
	IS[4] = (ME[2][i] - 1) * NF + 1;
	IS[5] = (ME[2][i] - 1) * NF + 2;
	IS[6] = (ME[2][i] - 1) * NF + 3;
}

// �鼯�ṹ�ն���
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
