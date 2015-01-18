#include <iostream>
#include <fstream>
#include <math.h>
#include <Eigen/Dense>

using namespace std;
using Eigen::MatrixXd;

#define I 5 // 初期点の数
#define K 200 // イテレーション回数

double __2n_minima (MatrixXd *);

int main()
{
	// n次元設定
	cout << "Dimention: ";
	int N;
	cin >> N;

	// PSO の初期点生成
	MatrixXd x(N, I);
	for (int i = 0; i < I; i++) {
		for (int n = 0; n < N; n++) {
			x(n, i) = ((double)rand() / ((double)RAND_MAX + 1)) * 10 - 5;
		}
	}
	// 一様乱数は、こっちに変更
	// MatrixXd x_* = MatrixXd::Random(n,n);

	// file 出力
	char filepath[256];
	sprintf(filepath, "n-dimensions/%d-dimensions.txt", N);
	ofstream fout;
	fout.open(filepath);
	fout << "k\tx_best" << endl;

	// 係数
	double λ, r1, r2;
	double λ_max = 1.0, λ_min = 0.6;
	double c1 = 0.0001, c2 = 0.01;

	// 行列の初期化
	MatrixXd y(N, I);
	MatrixXd f(N, I);
	MatrixXd z(N, I);
	MatrixXd p(N, I);
	MatrixXd q(N, I);
	MatrixXd dx(N, 1);
	MatrixXd v = MatrixXd::Zero(N, I);     // 移動方向ベクトル
	MatrixXd J = MatrixXd::Identity(N, N); // Hesse行列の近似行列の初期値は、単位行列
	MatrixXd G = MatrixXd::Zero(N, N);     // 更新行列の初期化

	MatrixXd J1, J2, J3, J4, J5;
	J1 = J2 = J3 = J4 = J5 = J;

	MatrixXd x_gbest(N, 1), x_temp(N, 1);
	MatrixXd p_temp(N, 1), q_temp(N, 1);
	double f_gbest, f_temp;

	for (int k = 0; k < K; k++) {
		// Inertia Weight Approach
		λ = λ_max - (λ_max - λ_min) * k / K;
		r1 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;
		r2 = ((double)rand() / ((double)RAND_MAX + 1)) * 1.0;

		z = f;
		for (int i = 0; i < I; i++) {
			for (int n = 0; n < N; n++) {
				f(n, i) = 4 * pow(x(n, i), 3) - 32 * x(n, i) + 5;
			}

			// PSO項のアルゴリズム
			x_temp = x.col(i);
			f_temp = __2n_minima(&x_temp);
			if ((k == 0 && i == 0) || f_gbest > f_temp) {
				x_gbest = x_temp;
				f_gbest = f_temp;
			}
		}

		if (k != 0) {
			p = v;
			q = f - z;

			for (int i = 0; i < I; i++) {
				// 準ニュートン項のアルゴリズム
				p_temp = p.col(i);
				q_temp = q.col(i);

				// DFP 公式
				switch (i) {
					case 0:
						J = J1;
						break;
					case 1:
						J = J2;
						break;
					case 2:
						J = J3;
						break;
					case 3:
						J = J4;
						break;
					case 4:
						J = J5;
						break;
					default: break;
				}
				G = (p_temp * p_temp.transpose() / (p_temp.transpose() * q_temp)(0, 0)) -
					(J * q_temp * q_temp.transpose() * J / (q_temp.transpose() * J * q_temp)(0, 0));

				if (((v.col(i)).transpose() * f.col(i))(0, 0) < 0) {
					// Hesseの近似行列の各成分
					J = J + G;
				} else {
					// 単位行列にリセット
					J = MatrixXd::Identity(N, N);
				}

				switch (i) {
					case 0:
						J1 = J;
						break;
					case 1:
						J2 = J;
						break;
					case 2:
						J3 = J;
						break;
					case 3:
						J4 = J;
						break;
					case 4:
						J5 = J;
						break;
					default: break;
				}
			}
		}

		for (int i = 0; i < I; i++) {
			dx = J * f.col(i);
			v.col(i) = λ * v.col(i) - c1*r1*dx + c2*r2*(x_gbest - x.col(i));
			x.col(i) = x.col(i) + v.col(i);
		}

		fout << k << "\t" << x_gbest.transpose() << endl;
	}

	return 0;
}

/*
 * 目的関数  2n-minima
 */
double __2n_minima (MatrixXd *A)
{
	double f = 0;
	for (int n = 0; n < A->rows(); n++) {
		double x = (*A)(n, 0);
		f += pow(x, 4) - 16 * pow(x, 2) + 5 * x;
	}
	return f;
}
