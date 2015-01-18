#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define I 5 // 初期点の数
#define N 2
#define M 1000
#define K 1000
#define ε 0.7

double __2n_minima (double x1, double x2);

int main()
{
	double x[I][N];
	int count = 0;

	double λ, r1, r2;
	double λ_max = 0.9, λ_min = 0.4;
	double c1 = 2, c2 = 2;

	char filepath[256];
	sprintf(filepath, "run/c1=%lf-c2=%lf-reset-εI.txt", c1, c2);
	ofstream fout; // file出力の為の定義
	fout.open(filepath);
	fout << "#m\tBest" << endl; // 見出し出力

	for (int m = 1; m <= M; m++) {
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

		double p[N], q[N];
		double f[I][N], z[I][N], dx[I][N];
		double v[I][N];
		for (int i = 0; i < I; i++) {
			for (int n = 0; n < N; n++) {
				v[i][n] = 0;
			}
		}

		// 近似行列の初期値は、単位行列
		double j11[I], j12[I], j21[I], j22[I];
		double g11[I], g12[I], g21[I], g22[I];
		for (int i = 0; i < I; i++) {
			j11[i] = j22[i] = 1;
			j12[i] = j21[i] = 0;
			g11[i] = g12[i] = g21[i] = g22[i] = 0;
		}

		for (int k = 0; k < K; k++) {
			// Inertia Weight Approach
			λ = λ_max - (λ_max - λ_min) * k / K;

			for (int i = 0; i < I; i++) {
				r1 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;
				r2 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;
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
					g11[i] = p[0]*p[0] / (p[0]*q[0]+p[1]*q[1]) - ((j11[i]*q[0]+j21[i]*q[1]) * (j11[i]*q[0]+j12[i]*q[1])) / (q[0]*(q[0]*j11[i]+q[1]*j21[i]) + q[1]*(q[0]*j12[i]+q[1]*j22[i]));
					g12[i] = p[0]*p[1] / (p[0]*q[0]+p[1]*q[1]) - ((j12[i]*q[0]+j22[i]*q[1]) * (j11[i]*q[0]+j12[i]*q[1])) / (q[0]*(q[0]*j11[i]+q[1]*j21[i]) + q[1]*(q[0]*j12[i]+q[1]*j22[i]));
					g21[i] = p[1]*p[0] / (p[0]*q[0]+p[1]*q[1]) - ((j11[i]*q[0]+j21[i]*q[1]) * (j21[i]*q[0]+j22[i]*q[1])) / (q[0]*(q[0]*j11[i]+q[1]*j21[i]) + q[1]*(q[0]*j12[i]+q[1]*j22[i]));
					g22[i] = p[1]*p[1] / (p[0]*q[0]+p[1]*q[1]) - ((j12[i]*q[0]+j22[i]*q[1]) * (j21[i]*q[0]+j22[i]*q[1])) / (q[0]*(q[0]*j11[i]+q[1]*j21[i]) + q[1]*(q[0]*j12[i]+q[1]*j22[i]));
				}

				// Hesseの近似行列の各成分
				j11[i] = j11[i] + g11[i];
				j12[i] = j12[i] + g12[i];
				j21[i] = j21[i] + g21[i];
				j22[i] = j22[i] + g22[i];

				dx[i][0] = j11[i] * f[i][0] + j12[i] * f[i][1];
				dx[i][1] = j21[i] * f[i][0] + j22[i] * f[i][1];

				if (dx[i][0] * f[i][0] + dx[i][1] * f[i][1] >= 0) {
					// H = εI + H
					j11[i] = j11[i] + ε;
					j22[i] = j22[i] + ε;
					// 移動方向の再計算
					dx[i][0] = j11[i] * f[i][0] + j12[i] * f[i][1];
					dx[i][1] = j21[i] * f[i][0] + j22[i] * f[i][1];
				}

				for (int n = 0; n < N; n++) {
					v[i][n] = λ*v[i][n] - c1*r1*dx[i][n] + c2*r2*(x_gbest[n] - x[i][n]);
					x[i][n] = x[i][n] + v[i][n];
				}
			}

			if (tmp_gbest <= (-78.332 * N) + 5) {
				fout << m << "\t";
				fout << tmp_gbest << endl;
				count++;
				break;
			} else if (k == K - 1) {
				fout << m << "\t";
				fout << tmp_gbest << endl;
			}
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
