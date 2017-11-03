/*
 * =====================================================================================
 *
 *       Filename:  test_allga.cpp
 *
 *    Description:  
 *
 *        Version:  1.0
 *        Created:  2016年10月28日 15时09分00秒
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Duan Junming (DSEC), duanjm@pku.edu.cn
 *   Organization:  
 *
 * =====================================================================================
 */


#include <iostream>
#include <cstdlib>
#include <fstream>
#include <mpi.h>
#include "Eigen/SparseCore"

int main(int argc, char* argv[]){
	int n_rank, rank;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &n_rank);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	Eigen::VectorXd x;
	x.resize(n_rank*3);
	int Ns[n_rank], Ne[n_rank];
	Ns[0] = 0; Ne[0] = 3;
	Ns[1] = 3; Ne[1] = 5;
	Ns[2] = 5; Ne[2] = 9;

	for(int i = 0; i < 3*n_rank; ++i) {
		if(i >= Ns[rank] && i < Ne[rank]) x(i) = i;
		else x(i) = 0;
	}

	std::cout << "qian: " << std::endl;
	std::cout << "rank: " << rank << " ,";
	for(int i = 0; i < 3*n_rank; ++i) {
		std::cout << x(i) << " ";
	}
	std::cout << std::endl;

	std::vector<int> Nrec;
	std::vector<int> disrec;
	Nrec.resize(n_rank);
	disrec.resize(n_rank);
	for(int i = 0; i < n_rank; ++i) {
		Nrec[i] = Ne[i] - Ns[i];
		disrec[i] = Ns[i];
	}
	MPI_Allgatherv(&x(Ns[rank]), Nrec[rank], MPI_DOUBLE,
			&x(0), &Nrec[0], &disrec[0], MPI_DOUBLE, MPI_COMM_WORLD);

	std::cout << "hou: " << std::endl;
	std::cout << "rank: " << rank << " ,";
	for(int i = 0; i < 3*n_rank; ++i) {
		std::cout << x(i) << " ";
	}
	std::cout << std::endl;


	MPI_Finalize();
	return 0;
}

