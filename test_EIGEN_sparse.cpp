#include "Eigen/Sparse"
#include "Eigen/SparseQR"
#include <vector>
#include <iostream>

typedef Eigen::Triplet<double> T;
int main(){
	int N(10);
	Eigen::SparseMatrix<double> A(N, N);
	std::vector< T > List;
	for(int i=0; i!=N; ++i){
		List.push_back( T(i, i-1, 1) );
	}
	A.setFromTriplets(List.begin(), List.end());
	std::cout<<"A:"<<"\n"<<A<<std::endl;

	return 0;
}

