#include <iostream>
#include <fstream>
#include <math.h>
#include <vector>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define I 5 // 初期点の数
#define M 100 // 試行回数
#define N 100 // 変数の次元数
#define K 1000 // イテレーション回数

void __initialPoints (VectorXd *, int *);
double __2n_minima (VectorXd *);

int main()
{
	ofstream fout;
	fout.open("n-dimensions/k-reset.txt");
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
	//double c1 = 0.03, c2 = 0.03;
	double c1 = 0.0001, c2 = 0.01;

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
		// Reset count
		int count = 0;

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
				// Differential function
				for (int n = 0; n < N; n++) {
					f[i](n) = 4 * pow(x[i](n), 3) - 32 * x[i](n) + 5;
				}

				// PSO Algorithm
				f_temp = __2n_minima(&(x[i]));
				if ((k == 0 && i == 0) || f_gbest > f_temp) {
					x_gbest = x[i];
					f_gbest = f_temp;
				}
			}

			if (k < 100 || k % 5 == 0) {
				fout << "\t" << count;
			}

			if (k != 0) {
				// quasi-Newton Algorithm
				for (int i = 0; i < I; i++) {
					p[i] = v[i];
					q[i] = f[i] - z[i];

					// DFP 公式
					G = (p[i] * p[i].transpose() / (p[i].transpose() * q[i])(0, 0)) -
						(J[i] * q[i] * q[i].transpose() * J[i] / (q[i].transpose() * J[i] * q[i])(0, 0));
					J[i] = J[i] + G;
					dx[i] = J[i] * f[i];

					if (((dx[i]).transpose() * f[i])(0, 0) >= 0) {
						count++;
						// 単位行列にリセット
						J[i] = MatrixXd::Identity(N, N);
						dx[i] = J[i] * f[i];
					}
				}
			}

			// Inertia Weight Approach
			λ = λ_max - (λ_max - λ_min) * k / K;
			for (int i = 0; i < I; i++) {
				r1 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;
				r2 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;

				v[i] = λ * v[i] - c1*r1*dx[i] + c2*r2*(x_gbest - x[i]);
				x[i] = x[i] + v[i];
			}
		}
		fout << endl;
	}

	return 0;
}

/*
 * 初期点群を返す
 */
void __initialPoints (VectorXd *X, int *g_line)
{
	// 使用する初期点file
	ifstream ifs("../../initialPoints/2n-minima.txt");
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
 * 目的関数  2n-minima
 */
double __2n_minima (VectorXd *X)
{
	double f = 0;
	for (int n = 0; n < N; n++) {
		double x = (*X)(n);
		f += pow(x, 4) - 16 * pow(x, 2) + 5 * x;
	}
	return f;
}
