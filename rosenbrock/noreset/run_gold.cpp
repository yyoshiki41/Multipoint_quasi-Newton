#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define I 5 // 初期点の数
#define N 2
#define M 1000
#define K 1000

double __rosenbrock (double x1, double x2);
void gold (double *x1, double *x2, double *dx1, double *dx2);

int main()
{
	double x[I][N];
	int count = 0;
	char filepath[256];
	sprintf(filepath, "./gold-c2=2.txt");
	ofstream fout; // file出力の為の定義
	fout.open(filepath); // fileを開く
	fout << "#m\tBest" << endl; // 見出し出力

	double λ, r2;
	double λ_max = 0.9, λ_min = 0.4;
	//double c2 = 0.003;
	double c2 = 2;

	for (int m = 1; m <= M; m++) {
		for (int i = 0; i < I; i++) {
			for (int n = 0; n < N; n++) {
				x[i][n] = ((double)rand() / ((double)RAND_MAX + 1)) * 10 - 5;
			}
		}

		double tmp_x, tmp_gbest, x_gbest[N];
		for (int i = 0; i < I; i++) {
			tmp_x = __rosenbrock(x[i][0], x[i][1]);
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
				for (int n = 0; n < N; n++) {
					z[i][n] = f[i][n];
				}
				f[i][0] = 400*x[i][0]*(x[i][0]*x[i][0] - x[i][1]) + 2*(x[i][0] - 1);
				f[i][1] = 200*(x[i][1] - x[i][0]*x[i][0]);

				if (k != 0) {
					// PSOのアルゴリズム
					tmp_x = __rosenbrock(x[i][0], x[i][1]);
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
				gold(&x[i][0], &x[i][1], &dx[i][0], &dx[i][1]);

				r2 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;
				for (int n = 0; n < N; n++) {
					v[i][n] = λ*v[i][n] - dx[i][n] + c2*r2*(x_gbest[n] - x[i][n]);
					x[i][n] = x[i][n] + v[i][n];
				}
			}

			if (tmp_gbest <= 0 + 0.01) {
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
double __rosenbrock (double x1, double x2)
{
	double f;
	f = pow((1-x1), 2) + 100 * pow((x2-x1*x1), 2);
	return f;
}

// 黄金分割法
void gold (double *x1, double *x2, double *dx1, double *dx2)
{
	double x1_min, x2_min, x1_max, x2_max;
	double a1, a2, b1, b2;
	double t_dx1, t_dx2;
	double f1, f2;

	int count = 0;
	double eps = 0.01;
	double tau = (sqrt(5) - 1) / 2;

	double norm = sqrt((*dx1)*(*dx1) + (*dx2)*(*dx2));

	t_dx1 = *dx1 / norm;
	t_dx2 = *dx2 / norm;
	x1_min = *x1, x2_min = *x2;
	x1_max = x1_min - t_dx1, x2_max = x2_min - t_dx2;

	while (count < 100) {
		count += 1;
		a1 = x1_min - (1-tau) * t_dx1;
		a2 = x2_min - (1-tau) * t_dx2;
		b1 = x1_min - tau * t_dx1;
		b2 = x2_min - tau * t_dx2;

		f1 = __rosenbrock(a1, a2);
		f2 = __rosenbrock(b1, b2);

		if (f1 < f2) {
			x1_max = b1;
			x2_max = b2;
			*dx1 = (*x1) - x1_max;
			*dx2 = (*x2) - x2_max;
		} else {
			x1_min = a1;
			x2_min = a2;
			*dx1 = (*x1) - x1_min;
			*dx2 = (*x2) - x2_min;
		}
		t_dx1 = a1 - b1;
		t_dx2 = a2 - b2;
		norm = sqrt(t_dx1*t_dx1 + t_dx2*t_dx2);
		if (norm < eps) break;
	}
}
