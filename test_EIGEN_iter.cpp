#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

double SHRINK(double x, double tau) {
	double sign(0);
	if(fabs(x) < 1e-15) sign = 0;
	else if(x > 0) sign = 1.0;
	else sign = -1.0;
	return sign * std::max(std::abs(x)-tau, 0.0);
}

typedef Eigen::Triplet<double> T;
int main(){
	u_int N(512), m(15);

	double lambda = 0.5;
	double mu = 0.5;
	double delta = 1.0;

	Eigen::SparseMatrix<double> A(N*N, N*N);
	Eigen::VectorXd u(N*N);
	Eigen::VectorXd f(N*N);
	Eigen::VectorXd d1(N*N);
	Eigen::VectorXd b1(N*N);
	Eigen::VectorXd d2(N*N);
	Eigen::VectorXd b2(N*N);
	Eigen::VectorXd RHS(N*N);
	std::vector<std::vector<double> > K;
	K.resize(m);
	for(u_int i = 0; i < m; ++i)
		K[i].resize(m);
	std::ifstream input("f.dat", std::ios::in);
	if(!input) std::cerr << "Opening file failed!!!" << std::endl;
	for(u_int i = 0; i < N*N; ++i) {
		input >> f(i);
	}
	input.close();
	std::ifstream inputK("k.dat", std::ios::in);
	if(!inputK) std::cerr << "Opening file failed!!!" << std::endl;
	for(u_int i = 0; i < m; ++i) {
		for(u_int j = 0; j < m; ++j) {
			inputK >> K[i][j];
		}
	}
	inputK.close();

	std::vector< T > List;
	u_int row, col;
	for(u_int l = 0; l != N; ++l){//col of u
		for(u_int k = 0; k != N; ++k){//row of u
			for(u_int i = 0; i < m; ++i) {
				for(u_int j = 0; j < m; ++j) {
					row = k - (m-1)/2 + i;
					col = l - (m-1)/2 + j;
					if(col >= 0 && col < N && row >= 0 && row < N)
						List.push_back( T(l*N + k, col*N + row, K[i][j]) );
				}
			}
		}
	}
	A.setFromTriplets(List.begin(), List.end());
	std::cout << A.rows() << std::endl;
	std::cout << A.cols() << std::endl;
	std::cout << A.nonZeros() << std::endl;

	Eigen::SparseMatrix<double> B(N*N, N*N);
	Eigen::SparseMatrix<double> laplace(N*N, N*N);
	B = A*A;
	List.clear();
	for(u_int i = 0; i < N; ++i) {
		for(u_int j = 0; j < N; ++j) {
			List.push_back( T(j*N + i, j*N + i, -4) );
			if(j > 0) List.push_back( T(j*N + i, (j-1)*N + i, 1) );
			if(j < N-1) List.push_back( T(j*N + i, (j+1)*N + i, 1) );
			if(i > 0) List.push_back( T(j*N + i, j*N + i+1, 1) );
			if(i < N-1) List.push_back( T(j*N + i, j*N + i-1, 1) );
		}
	}
	laplace.setFromTriplets(List.begin(), List.end());
	B = A + mu * laplace;

	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
	solver.compute(B);
	u_int MaxIte(10);
	for(u_int ite = 0; ite < MaxIte; ++ite) {
		//求解线性方程组
		RHS = A * f;
    for(u_int j = 0; j < N; ++j) {
      for(u_int i = 0; i < N; ++i) {
				if(i != N-1) RHS[j*N+i] += mu * (d1[j*N+i+1] - d1[j*N+i] - b1[j*N+i+1] + b1[j*N+i]);
        else RHS[j*N+i] = RHS[j*N+i-1];
				if(j != N-1) RHS[j*N+i] += mu * (d2[(j+1)*N+i] - d2[j*N+i] - b2[j*N+i+1] + b2[j*N+i]);
        else RHS[j*N+i] = RHS[(j-1)*N+i];
			}
		}
		u = solver.solve(RHS);
    std::cout << "iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;

		//求解d_k+1
		for(u_int j = 0; j < N; ++j) {
			for(u_int i = 0; i < N; ++i) {
				if(i != N-1) d1[i] = SHRINK(u[j*N+i+1] - u[j*N+i] + b1[j*N+i], lambda/mu);
        else d1[i] = SHRINK(u[j*N+i] - u[j*N+i-1] + b1[j*N+i], lambda/mu);
				if(j != N-1) d2[i] = SHRINK(u[(j+1)*N+i] - u[j*N+i] + b2[j*N+i], lambda/mu);
        else d2[i] = SHRINK(u[j*N+i] - u[(j-1)*N+i] + b2[j*N+i], lambda/mu);
			}
		}
		//求解b_k+1
		for(u_int j = 0; j < N; ++j) {
			for(u_int i = 0; i < N; ++i) {
				if(i != N-1) b1[i] += delta * (u[j*N+i+1] - u[j*N+i] - d1[j*N+i]);
        else b1[i] = delta * (u[j*N+i+1] - u[j*N+i] - d1[j*N+i]);
				if(j != N-1) b2[i] += delta * (u[(j+1)*N+i] - u[j*N+i] - d2[j*N+i]);
        else b2[i] = delta * (u[j*N+i] - u[(j-1)*N+i] - d2[j*N+i]);
			}
		}

    std::cout << "ite:" << ite << std::endl;
	}

	return 0;
}

