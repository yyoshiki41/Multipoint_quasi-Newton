#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define I 5 // 初期点の数
#define N 2
#define K 200

double __2n_minima (double x1, double x2);

int main()
{
	double x[I][N];
	int count = 0;
	char filepath[256];
	sprintf(filepath, "quasiPso_2n/run.txt");
	ofstream fout; // file出力の為の定義
	fout.open(filepath); // fileを開く
	fout << "#m\tx1\tx2" << endl; // 見出し出力

	for (int m = 1; m <= 100; m++) {
		for (int i = 0; i < I; i++) {
			for (int n = 0; n < N; n++) {
				x[i][n] = ((double)rand() / ((double)RAND_MAX + 1)) * 10 - 5;
			}
		}

		double tmp_x, tmp_gbest, x_gbest[N];
		for (int i = 0; i < I; i++) {
			tmp_x = __2n_minima(x[i][0], x[i][1]);
			if (i == 0 || tmp_gbest > tmp_x) {
				for (int n = 0; n < N; n++) {
					x_gbest[n] = x[i][n];
				}
				tmp_gbest = tmp_x;
			}
		}

		double p[N], q[N], dx[N];
		double f[I][N], z[I][N];
		double v[I][N];
		for (int i = 0; i < I; i++) {
			for (int n = 0; n < N; n++) {
				v[i][n] = 0;
			}
		}
		double c1 = 0.0007, c2 = 0.001;
		double λ, r1, r2;
		double λ_max = 1.1, λ_min = 0.9;

		// 近似行列の初期値は、単位行列
		double j11 = 1, j12 = 0, j21 = 0, j22 = 1;
		double g11 = 0, g12 = 0, g21 = 0, g22 = 0;

		for (int k = 0; k < K; k++) {
			// Inertia Weight Approach
			λ = λ_max - (λ_max - λ_min) * k / K;
			r1 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;
			r2 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;

			for (int i = 0; i < I; i++) {
				for (int n = 0; n < N; n++) {
					z[i][n] = f[i][n];
					f[i][n] = 4 * pow(x[i][n], 3) - 32 * x[i][n] + 5;
				}

				if (k != 0) {
					// PSOのアルゴリズム
					tmp_x = __2n_minima(x[i][0], x[i][1]);
					if (tmp_gbest > tmp_x) {
						for (int n = 0; n < N; n++) {
							x_gbest[n] = x[i][n];
						}
						tmp_gbest = tmp_x;
					}

					for (int n = 0; n < N; n++) {
						p[n] = v[i][n];
						q[n] = f[i][n] - z[i][n];
					}
					g11 = p[0]*p[0] / (p[0]*q[0]+p[1]*q[1]) - ((j11*q[0]+j21*q[1]) * (j11*q[0]+j12*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
					g12 = p[0]*p[1] / (p[0]*q[0]+p[1]*q[1]) - ((j12*q[0]+j22*q[1]) * (j11*q[0]+j12*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
					g21 = p[1]*p[0] / (p[0]*q[0]+p[1]*q[1]) - ((j11*q[0]+j21*q[1]) * (j21*q[0]+j22*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
					g22 = p[1]*p[1] / (p[0]*q[0]+p[1]*q[1]) - ((j12*q[0]+j22*q[1]) * (j21*q[0]+j22*q[1])) / (q[0]*(q[0]*j11+q[1]*j21) + q[1]*(q[0]*j12+q[1]*j22));
				}

				if (v[i][0] * f[i][0] + v[i][1] * f[i][1] < 0) {
					// Hesseの近似行列の各成分
					j11 = j11 + g11;
					j12 = j12 + g12;
					j21 = j21 + g21;
					j22 = j22 + g22;
				} else {
					// 単位行列にリセット
					j11 = j22 = 1;
					j12 = j21 = 0;
				}
				dx[0] = j11 * f[i][0] + j12 * f[i][1];
				dx[1] = j21 * f[i][0] + j22 * f[i][1];

				for (int n = 0; n < N; n++) {
					v[i][n] = λ*v[i][n] - c1*r1*dx[n] + c2*r2*(x_gbest[n] - x[i][n]);
					x[i][n] = x[i][n] + v[i][n];
				}
			}
		}

		fout << m << "\t";
		fout << x_gbest[0] << "\t";
		fout << x_gbest[1] << endl;
		if (tmp_gbest <= -150) {
			count++;
		}
	}
	fout << count << endl;

	return 0;
}

// 目的関数
double __2n_minima (double x1, double x2)
{
	double f;
	f = pow(x1, 4) - 16 * pow(x1, 2) + 5 * x1 + pow(x2, 4) - 16 * pow(x2, 2) + 5 * x2;
	return f;
}
