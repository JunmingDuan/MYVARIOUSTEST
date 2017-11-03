/*
 * =====================================================================================
 *
 *       Filename:  openmp_test.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016年10月08日 20时15分11秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */

#include <cstdlib>
#include <omp.h>
#include "MAT.h"


int main () {
	const int N(1e3);
	int n_threads(2);
	MAT<int> a(N,N);
	//MAT<int> new_a(N,N);
	double t1 =	omp_get_wtime();
	clock_t t3 = clock();
	for(u_int i = 0; i < N; ++i)
#pragma omp parallel for num_threads(n_threads)
		for(u_int j = 0; j < N; ++j) {
			for(u_int k = 0; k < N; ++k)
				a[i][j] = a[i][j] + 1;
		}

	//for(u_int i = 0; i < N; ++i)
		//for(u_int j = 0; j < N; ++j)
			//if(a[i][j] != N) std::cout << "Wrong!" << i << " " << j << " " << a[i][j] << std::endl;

	double t2 =	omp_get_wtime();
	clock_t t4 = clock();

	std::cout << (t2-t1) << std::endl;
	std::cout << (t4-t3)/CLOCKS_PER_SEC << std::endl;

	return 0;
}

