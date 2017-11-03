#ifndef MAT_H
#define MAT_H
//矩阵

#include "vvector.h"
#include <cmath>

enum{ NOR=0,SYM=1,SPD=2 };

template<class T>
class MAT : public vvector< vvector<T> >{
	private:
		typedef vvector< vvector<T> > Base;

		typedef vvector<int> SWAP;

	protected:
		int row;//行数	

		int col;//列数

		int flag;
	public:
		MAT(){
		}

		MAT(const int m1,const int n1,const T a=0,int FLAG=NOR):row(m1),col(n1),flag(FLAG){
			vvector<T> v(n1,a);
			Base::assign(m1,v);
		}

		void initial(T a){
			for(int i=0;i!=row;++i){
				for(int j=0;j!=col;++j){
					(*this)[i][j]=a;
				}
			}
		}

		void resize(int m1,int n1,T a=0){
			row=m1,col=n1;
			Base::resize(row);
			for(int i=0;i!=row;++i){
				(*this)[i].resize(col,a);
			}
		}

		int N_row() const{//返回行数
			return row;
		}
		
		int N_col() const{//返回列数
			return col;
		}

		void rowSWAP(int i,int j){//交换两行
			std::swap((*this)[i],(*this)[j]);
		}

		void colSWAP(int i,int j){//交换两列
			for(size_t k=0;k!=row;++k){
				std::swap((*this)[k][i],(*this)[k][j]);
			}
		}

		MAT<T> trans(){
			MAT<T> trans((*this).col,(*this).row);
			for(int i=0;i!=(*this).row;++i){
				for(int j=0;j!=(*this).col;++j){
					trans[j][i]=(*this)[i][j];
				}
			}

			return trans;
		}

		void Regular(){//A^TA
			MAT<T> X(col,col);
			for(int i=0;i!=col;++i){
				for(int j=0;j!=col;++j){
					for(int k=0;k!=row;++k){
						X[i][j]+=(*this)[k][i]*(*this)[k][j];
					}
				}
			}

			resize(col,col);
			(*this)=X;
			flag=SYM;
		}

		int rowMAX(int i,int j,int k){//第k列的第i到第j个元素中最大的，并返回下标
			T max;
			int t=i;
			max=std::abs((*this)[i][k]);
			for(int p=i+1;p!=j+1;++p){
				if(max<std::abs((*this)[p][k])){
					max=std::abs((*this)[p][k]);
					t=p;
				}
			}
			return t;
		}

		void nummultiply(int i,int j,int k,int m,T x)//数乘，将i-j行，k-m列乘以x
		{
			for(int p=i;p!=j+1;++p){
				for(int q=k;q!=m+1;++q){
					(*this)[p][q]*=x;
				}
			}
		}

		T NormInfinite(){
			vvector<T> sum(row);
			for(int i=0;i!=row;++i){
				sum[i]=(*this)[i].Norm1();
			}
			return sum.NormInfinite();
		}
		
		T AbsMax(){
			T a=0;
			for(int i=0;i!=row;++i){
				for(int j=0;j!=col;++j){
					if(std::abs((*this)[i][j])>a){
						a=std::abs((*this)[i][j]);
					}
				}
			}
			return a;
		}

		template<class C>
		vvector<T> arrmultiply(const vvector<C> &b)//将矩阵乘以b，输出
		{
			vvector<T> c(row,0);
			for(int i=0;i!=row;++i){
				for(int j=0;j!=col;++j){
					c[i]+=(*this)[i][j]*b[j];
				}
			}
			return c;
		}

		bool gauss(){//Gauss消去法，L为单位下三角矩阵，主对角为1未保存，其余保存在原矩阵下半部分，U为上三角矩阵，保存在原矩阵上半部分
			int dim;
			if(row>col) dim=col;
			else dim=row;
			for(int k=0;k!=dim-1;++k){
				if((*this)[k][k]==0){
					std::cerr<<"对角为0"<<std::endl;
					return 0;
				}
				
				for(int i=k+1;i!=row;++i){
					(*this)[i][k]/=(*this)[k][k];
					for(int j=k+1;j!=col;++j){
						(*this)[i][j]-=(*this)[i][k]*(*this)[k][j];
					}
				}
			}
			return 1;
		}

		SWAP colgauss(){//列主元Gauss消去法，形参保存行交换信息，其它如上
			int dim;
			if(row>col) dim=col;
			else dim=row;

			SWAP a(dim,-1);

			for(int k=0;k!=dim-1;++k){
				int p=rowMAX(k,row-1,k);
				if((*this)[p][k]==0){
					std::cerr<<"对角为0"<<std::endl;
					return 0;
				}
				
				rowSWAP(p,k);
				a[k]=p;

				for(int i=k+1;i!=row;++i){
					(*this)[i][k]/=(*this)[k][k];
					for(int j=k+1;j!=col;++j){
						(*this)[i][j]-=(*this)[i][k]*(*this)[k][j];
					}
				}
			}
			return a;
		}
		
		bool cholesky(){//分解后下三角阵L保存在原矩阵下三角部分
			if(flag!=SYM){
				std::cerr<<"非对称"<<std::endl;
				return 0;
			}
			for(int i=0;i!=col;++i){
				if((*this)[i][i]<=0){
					std::cerr<<"非正定"<<std::endl;
					return 0;
				}
				(*this)[i][i]=sqrt((*this)[i][i]);
				nummultiply(i+1,col-1,i,i,1/(*this)[i][i]);
				for(int j=i+1;j!=col;++j){
					for(int k=j;k!=col;++k){
						(*this)[k][j]-=(*this)[k][i]*(*this)[j][i];
					}
				}
			}
			return 1;
		}

		bool imcholesky(){//分解后单位下三角阵L，主对角为1未保存，其余部分保存在原矩阵严格下三角部分，D保存在原矩阵主对角
			if(flag!=SYM){
				std::cerr<<"非对称"<<std::endl;
				return 0;
			}
			for(int i=0;i!=col;++i){
				if((*this)[i][i]<=0){
					std::cerr<<"非正定"<<std::endl;
					return 0;
				}
				double a[col];
				for(int j=0;j!=i;++j){
					a[j]=(*this)[i][j]*(*this)[j][j];

					(*this)[i][i]-=(*this)[i][j]*a[j];
					for(int k=i+1;k!=col;++k){
						(*this)[k][i]-=(*this)[k][j]*a[j];
					}
				}
				for(int k=i+1;k!=col;++k){
					(*this)[k][i]/=(*this)[i][i];
				}
			}
			return 1;
		}

		//MAT<T> Householder(const T beta,const vvector<T>& v){
		//	int n=v.size();
		//	MAT<T> H(n,n);
		//	for(int i=0;i!=n;++i){
		//		H[i][i]=1;
		//	}
		//	if(beta!=0){
		//		for(int i=0;i!=n;++i){
		//			for(int j=0;j!=n;++j){
		//				H[i][j]-=beta*v[i]*v[j];
		//			}
		//		}
		//	}
		//	return H;

		//}

		vvector<T> QR(){
			vvector<T> d(std::min(row-1,col));
			for(int j=0;j!=col;++j){
				if(j<row-1){
					vvector<T> ANext(row-j);
					for(int i=j;i!=row;++i){
						ANext[i-j]=(*this)[i][j];
					}
					vvector<T> v(ANext.size());
					d[j]=ANext.Householder(v);
					vvector<T> w(col-j);
					for(int i=0;i!=w.size();++i){
						for(int k=0;k!=v.size();++k){
							w[i]+=(*this)[k+j][i+j]*v[k];
						}
						w[i]*=d[j];
					}
					for(int i=0;i!=v.size();++i){
						for(int k=0;k!=w.size();++k){
							(*this)[i+j][k+j]-=v[i]*w[k];
						}
					}
					//MAT<T> H=Householder(d[j],v);
					//MAT<T> temp(v.size(),v.size());
					//for(int i=j;i!=row;++i){
					//	for(int k=j;k!=col;++k){
					//		for(int p=j;p!=row;++p){
					//			temp[i-j][k-j]+=H[i-j][p-j]*(*this)[p][k];
					//		}
					//	}
					//}
					//for(int i=j;i!=row;++i){
					//	for(int k=j;k!=col;++k){
					//		(*this)[i][k]=temp[i-j][k-j];
					//	}
					//}
					for(int i=j+1;i!=row;++i){
						(*this)[i][j]=v[i-j];
					}
				}
			}
			return d;
		}

		MAT<T> QR_Q(vvector<T>& d0){
			int N=std::min(row-1,col);
			MAT<T> Q(row,row);
			for(int i=0;i!=row;++i){
				Q[i][i]=1;
			}
			for(int j=N-1;j!=-1;--j){//左乘H_j
				vvector<T> v(row-j);
				v[0]=1;
				for(int i=1;i!=v.size();++i){
					v[i]=(*this)[i+j][j];
				}
				vvector<T> w(row-j);
				for(int i=0;i!=w.size();++i){
					for(int k=0;k!=v.size();++k){
						w[i]+=Q[k+j][i+j]*v[k];
					}
					w[i]*=d0[j];
				}
				for(int i=0;i!=v.size();++i){
					for(int k=0;k!=w.size();++k){
						Q[i+j][k+j]-=v[i]*w[k];
					}
				}

				//MAT<T> H_i=Householder(d0[i],v);
				//MAT<T> temp(v.size(),v.size());
				//for(int j=i;j!=row;++j){
				//	for(int k=i;k!=row;++k){
				//		for(int p=i;p!=row;++p){
				//			temp[j-i][k-i]+=H_i[j-i][p-i]*Q[p][k];
				//		}
				//	}
				//}
				//for(int j=i;j!=row;++j){
				//	for(int k=i;k!=row;++k){
				//		Q[j][k]=temp[j-i][k-i];
				//	}
				//}
			}

			return Q;
		}

		T Norm1Optimi(){
			vvector<T> x(col,1./col),w(row),z(col);
			vvector<int> v(row);
			while(1){
				w=arrmultiply(x);
				v=w.sign();
				z=((*this).trans()).arrmultiply(v);
				
				if(z.NormInfinite()<=z.arrmultiply(x)){
					return w.Norm1();
				}
				else {
					x.assign(col,0);
					T Norm=z.NormInfinite();
					for(int i=0;i!=z.size();++i){
						if(std::abs(z[i])==Norm) { x[i]=1; break; }
					}
				}
			}

			std::cerr<<"Cannot stop"<<std::endl;
			return 0;
		}

		template<class type> friend std::ostream& operator<<(std::ostream&,const MAT<type>&);
};

template<class T>
std::ostream&
operator<<(std::ostream& out,const MAT<T>& A)
{
	for(int i=0;i!=A.N_row();++i){
		out<<A[i]<<"\n";
	}
	return out;
}

template<class T>
MAT<T> operator+(const MAT<T>& v1, const MAT<T>& v2)
{
	if(v1.N_row() != v2.N_row()){
		std::cerr<<"different dim"<<std::endl;
		return v1;
	}
    MAT<T> v3(v1);
	for(int i=0;i!=v1.N_row();++i){
		v3[i] += v2[i];
	}
    return v3;
}

template<class T>
MAT<T> operator-(const MAT<T>& v1, const MAT<T>& v2)
{
	if(v1.N_row() != v2.N_row()){
		std::cerr<<"different dim"<<std::endl;
		return v1;
	}
    MAT<T> v3(v1);
	for(int i=0;i!=v1.N_row();++i){
		v3[i] -= v2[i];
	}
    return v3;
}

template<class T>
MAT<T> operator*(const MAT<T>& v1, const double v2)
{
    MAT<T> v3(v1);
	for(int i=0;i!=v1.N_row();++i){
		v3[i] *= v2;
	}
    return v3;
}

template<class T>
MAT<T> operator/(const MAT<T>& v1, const double v2)
{
    MAT<T> v3(v1);
	for(int i=0;i!=v1.N_row();++i){
		v3[i] /= v2;
	}
    return v3;
}

#endif
