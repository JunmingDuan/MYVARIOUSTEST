/*
 * =====================================================================================
 *
 *       Filename:  test_eigenv.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016年11月12日 11时04分08秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */


#include <vector>
#include <iostream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include "../Eigen/Dense"

int main(){
	int N(20);
	//EIGEN求Hilbert矩阵特征值
	Eigen::MatrixXd A;
	A.resize(N, N);
	for(u_int i = 0; i < N; ++i) {
		for(u_int j = 0; j < N; ++j) {
			A(i,j) = 1./(i+j+1);
		}
	}
	std::cout << A.eigenvalues() << std::endl;
	//GSL求Hilbert矩阵特征值
	gsl_matrix* m = gsl_matrix_alloc(N, N);
	for(u_int i = 0; i < N; ++i) {
		for(u_int j = 0; j < N; ++j) {
			gsl_matrix_set(m, i, j, 1./(i+j+1));
		}
	}
	gsl_vector * eval = gsl_vector_alloc(N);
	gsl_matrix * evec = gsl_matrix_alloc(N, N);
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(N);
	gsl_eigen_symmv(m, eval, evec, w);
	gsl_eigen_symmv_free(w);
	gsl_vector_fprintf(stdout, eval, "%.16g");

	gsl_vector_free(eval);
	gsl_matrix_free(evec);

	return 0;
}

