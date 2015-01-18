#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define I 5    // 初期点の数
#define M 1000  // 試行回数
#define N 2  // 変数の次元数
#define K 1000 // イテレーション回数

void __initialPoints (VectorXd *, int *);
double __rosenbrock (VectorXd *);
void gold (VectorXd *, VectorXd *);

int main()
{
	ofstream fout;
	fout.open("./2k-value-gold.txt");
	// 見出し
	fout << "try";
	for (int k = 0; k < K; k++) {
		if (k < 100 || k % 5 == 0) {
			fout << "\t" << k;
		}
	}
	fout << endl;

	// 係数
	double λ, r1, r2;
	double λ_max = 0.9, λ_min = 0.4;
	double c2 = 0.03;

	int g_line = 1;
	for (int m = 1; m <= M; m++) {
		// Initial Points
		vector<VectorXd> x(I);
		// Initialize Matrix
		vector<VectorXd> y(I);
		vector<VectorXd> f(I);
		vector<VectorXd> z(I);
		vector<VectorXd> p(I);
		vector<VectorXd> q(I);
		vector<VectorXd> v(I);
		vector<VectorXd> dx(I);
		// Approximate matrix of the Hessian of the inverse matrix
		vector<MatrixXd> J(I);
		// Update matrix in the DFP
		MatrixXd G = MatrixXd::Zero(N, N);
		// The best point in the group
		VectorXd x_gbest(N);
		double f_gbest, f_temp;

		for (int i = 0; i < I; i++) {
			x[i] = VectorXd(N);
			__initialPoints (&(x[i]), &g_line);

			J[i] = MatrixXd::Identity(N, N);
			f[i] = VectorXd::Zero(N);
			v[i] = VectorXd::Zero(N);
			dx[i] = VectorXd::Zero(N);
		}

		fout << m;
		for (int k = 0; k < K; k++) {
			z = f;
			for (int i = 0; i < I; i++) {
				for (int n = 0; n < N; n++) {
					if (n == 0) {
						f[i](n) = 400*x[i](n)*(x[i](n)*x[i](n) - x[i](n+1)) + 2*(x[i](n) - 1);
					} else if (n + 1 == N) {
						f[i](n) = 200*(x[i](n) - x[i](n-1)*x[i](n-1));
					} else {
						f[i](n) = 400*x[i](n)*(x[i](n)*x[i](n) - x[i](n+1)) + 2*(x[i](n) - 1) + 200*(x[i](n) - x[i](n-1)*x[i](n-1));
					}
				}

				// PSO Algorithm
				f_temp = __rosenbrock(&(x[i]));
				if ((k == 0 && i == 0) || f_gbest > f_temp) {
					x_gbest = x[i];
					f_gbest = f_temp;
				}
			}

			if (k < 100 || k % 5 == 0) {
				fout << "\t" << f_gbest;
			}

			if (k != 0) {
				// quasi-Newton Algorithm
				for (int i = 0; i < I; i++) {
					p[i] = v[i];
					q[i] = f[i] - z[i];

					// DFP
					G = (p[i] * p[i].transpose() / (p[i].transpose() * q[i])(0, 0)) -
						(J[i] * q[i] * q[i].transpose() * J[i] / (q[i].transpose() * J[i] * q[i])(0, 0));
					J[i] = J[i] + G;
					dx[i] = J[i] * f[i];
					gold(&x[i], &dx[i]);
				}
			}

			// Inertia Weight Approach
			λ = λ_max - (λ_max - λ_min) * k / K;
			for (int i = 0; i < I; i++) {
				r2 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;

				v[i] = λ * v[i] - dx[i] + c2*r2*(x_gbest - x[i]);
				x[i] = x[i] + v[i];
			}
		}
		fout << endl;
	}

	return 0;
}

/*
 * Return initial points
 */
void __initialPoints (VectorXd *X, int *g_line)
{
	// 使用する初期点file
	ifstream ifs("../../../initialPoints/rosenbrock.txt");
	string str;
	char *e;
	double r;
	if (ifs.fail()) {
		cout << "初期点Fileの読み込み失敗" << endl;
		exit(1);
	}

	int l_line = 1, n = 0;
	while (getline(ifs, str)) {
		r = strtod(str.c_str(), &e);

		if (strlen(e) != 0) {
			continue;
		} else if (l_line >= (*g_line)) {
			(*X)(n) = r;
			(*g_line)++;
			n++;
			if (n >= N) break;
		}
		l_line++;
	}
}

/*
 * Rosenbrock
 */
double __rosenbrock (VectorXd *X)
{
	double f = 0;
	for (int n = 0; n < N-1; n++) {
		double x0 = (*X)(n), x1 = (*X)(n+1);
		f += 100*(x1 - x0*x0)*(x1 - x0*x0) + (1 - x0) * (1 - x0);
	}
	return f;
}

/*
 * 黄金分割法
 */
void gold (VectorXd *X, VectorXd *dx)
{
	double f1, f2;
	int count = 0;
	double eps = 0.01;
	double tau = (sqrt(5) - 1) / 2;
	double norm = (*dx).lpNorm<N>();
	VectorXd temp1(N), temp2(N);

	VectorXd X_0(N);
	X_0 = (*X);
	VectorXd dX(N);
	dX = (*dx) / norm;
	VectorXd X_max(N);
	X_max = X_0 - dX;

	while (count < 100) {
		count += 1;
		temp1 = X_0 - (1-tau) * dX;
		f1 = __rosenbrock(&temp1);
		temp2 = X_0 - tau * dX;
		f2 = __rosenbrock(&temp2);

		if (f1 < f2) {
			X_max = temp2;
			*dx = (*X) - X_max;
		} else {
			X_0 = temp1;
			*dx = (*X) - X_0;
		}
		dX = (*X) - X_max;
		norm = dX.lpNorm<N>();
		if (norm < eps) break;
	}
}
