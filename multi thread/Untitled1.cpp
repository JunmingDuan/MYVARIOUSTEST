#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

#define DOS 4
#define M 3
#define GAMMA 1.4
#define N  ( (M+1)*DOS )
#define CFL 0.05
#define gpn 3
#define ncell 80

class MeshCell
{
public:
	double solu[N];
	double UGL1[N];
	double UGL2[N];
	double UGL3[N];
	double UGL4[N];
	double AUdP[N];
	double oldU[N];
};




double global_w[20]={0.1527533871307258	,
	0.1527533871307258	,
	0.1491729864726037	,
	0.1491729864726037	,
	0.1420961093183820	,
	0.1420961093183820	,
	0.1316886384491766	,
	0.1316886384491766	,
	0.1181945319615184	,
	0.1181945319615184	,
	0.1019301198172404	,
	0.1019301198172404	,
	0.0832767415767048	,
	0.0832767415767048	,
	0.0626720483341091	,
	0.0626720483341091	,
	0.0406014298003869	,
	0.0406014298003869	,
	0.0176140071391521	,
	0.0176140071391521	};

double global_p[20] = {	-0.0765265211334973,
	0.0765265211334973,
	-0.2277858511416451,
	0.2277858511416451,
	-0.3737060887154195,
	0.3737060887154195,
	-0.5108670019508271,
	0.5108670019508271,
	-0.6360536807265150,
	0.6360536807265150,
	-0.7463319064601508,
	0.7463319064601508,
	-0.8391169718222188,
	0.8391169718222188,
	-0.9122344282513259,
	0.9122344282513259,
	-0.9639719272779138,
	0.9639719272779138,
	-0.9931285991850949,
	0.9931285991850949};


bool RKUpdate1D ( MeshCell * mesh1D, const int Ncell, double lamx );
double  calallAU( double *hU, double  allAU[N][N] , bool flag);
double calallAUM (double *hUL, double *hUR, double  allAUM0[N][N], double  allAUM1[N][N]);
double calallhA ( double *hU, double  allhA0[N][N], double  allhA1[N][N] );
double calhAU( double *hU, double xi, double  hA0[N][N], double  hA1[N][N] );
void calAU( double U[DOS], double  AU0[DOS][DOS], double  AU1[DOS][DOS] , double & maxeig) ;
void SolverAxb(double AA[N][N], double b[N], double  x[N]) ;
double WENOGL41(double v1,double v2,double v3,double v4,double v5);
double WENOGL42(double v1,double v2,double v3,double v4,double v5);
double WENOGL43(double v1,double v2,double v3,double v4,double v5);
double WENOGL44(double v1,double v2,double v3,double v4,double v5);
bool checkallhU( double *hU );
double Phix(int n, double x);


double **Iget2Ddouble(int row,int col)   
{
	int i;
	double **y;
	y=new double *[row];
	if (!y) printf("No enough memory!");
	for (i=0;i<row;i++)
	{
		y[i]=new double [col];	
		if (!y) printf("No enough memory!");
	}
	return y;
}

void Idel2Ddouble(int row,double **p)  
{
	for(int i=0;i<row;i++)
		delete[] p[i];
	delete[] p;
}


double ***Iget3Dmatrix(int row,int col,int col2)   
{
	int i,k;
	double ***y;
	y=new double **[row];
	for (i=0;i<row;i++)
	{
		y[i]=new double *[col];
		for (k=0;k<col;k++)
			y[i][k] = new double [col2];
	}
	for ( i=0;i<row;i++)
		for (k=0; k<col; k++)
			for(int j=0;j<col2;j++)
				*(*(*(y+i)+k)+j)=1.;
	return y;
}



void Idel3Dmatrix(int row,int col,double ***p)  
{
	for(int i=0;i<row;i++)
	{
		for (int k=0;k<col;k++)
			delete[] p[i][k];
		delete[] p[i];
	}
	delete[] p;
}










bool RKUpdate1D ( MeshCell * mesh1D, const int Ncell, double lamx )
{
	double face_AP[Ncell-1][N], face_AN[Ncell-1][N];
	double allAUM0[N][N], allAUM1[N][N], allAUM[N][N];
	double AA[N][N],knowb[N],unknowx[N];
	double allAU1[N][N], allAU2[N][N], allAU3[N][N], allAU4[N][N];
	double dP1[N],dP2[N],dP3[N],dP4[N];
	
	int RK, k, i, j, kk;
	
	const double GLp1 = -0.5;
	const double GLp2 = -0.1*sqrt(5.);
	const double GLp3 = 0.1*sqrt(5.);
	const double GLp4 = 0.5;	
	
	

		// boundary condition 
		for ( i = 0; i < N; i++ )
			for ( k = 0; k < gpn ; k++ )
			{
				mesh1D[k].solu[i] = 	mesh1D[k+ncell].solu[i] ;
	
				mesh1D[Ncell-k-1].solu[i] = 	mesh1D[Ncell-k-1-ncell].solu[i] ;
			}


	
	for ( RK = 1; RK <=3; RK++ )
	{
	  //std::cout << "RK = " << RK << std::endl;
		// WENO reconstruction
		for ( k = 0; k < (Ncell -1); k++ )
		{
			if ( k>=2 && k<=(Ncell-3) )
				for ( i=0; i<N; i++ ) {
					mesh1D[k].UGL3[i] = WENOGL43 ( mesh1D[k-2].solu[i], mesh1D[k-1].solu[i],
						mesh1D[k].solu[i], mesh1D[k+1].solu[i], mesh1D[k+2].solu[i] );
					mesh1D[k].UGL4[i] = WENOGL44 ( mesh1D[k-2].solu[i], mesh1D[k-1].solu[i],
						mesh1D[k].solu[i], mesh1D[k+1].solu[i], mesh1D[k+2].solu[i] );
				}					

			if ( k>=1 && k<=(Ncell-4) )
				for ( i=0; i<N; i++ ) {
					mesh1D[k+1].UGL2[i] = WENOGL42 ( mesh1D[k-1].solu[i], mesh1D[k].solu[i],
						mesh1D[k+1].solu[i], mesh1D[k+2].solu[i], mesh1D[k+3].solu[i] );
					mesh1D[k+1].UGL1[i] = WENOGL41 ( mesh1D[k-1].solu[i], mesh1D[k].solu[i],
						mesh1D[k+1].solu[i], mesh1D[k+2].solu[i], mesh1D[k+3].solu[i] );
				}											
		}
		
		
		
		// Positivity-preserving limiter
		
		double **UU;
		UU = Iget2Ddouble ( 4, N );
		
		for ( k = 3; k < (Ncell -3); k++ )
		{
			if ( ! checkallhU( mesh1D[k].solu ) )
			{
				std::cout << k << ":negative p or rho" << std::endl;
				//exit(1);
				return 1;
			}
				
	
				
				
			for ( i = 0; i < N; i ++ )
			{
				UU[0][i] = mesh1D[k].UGL1[i];
				UU[1][i] = mesh1D[k].UGL2[i];
				UU[2][i] = mesh1D[k].UGL3[i];	
				UU[3][i] = mesh1D[k].UGL4[i];			 
			}
			bool flagU[4],flagU0[4];
			
		
			for (kk=0; kk<4;kk++)
				flagU[kk] = checkallhU( UU[kk] );
				
			
				
			if ( flagU[0] + flagU[1] + flagU[2] + flagU[3] < 4  )
			{
				std:: cout << "PP limiter is used\n";
				double ThetaL = 0.;
				double ThetaR = 1.;
				
				while ( (ThetaR - ThetaL) > 1.0e-13)
				{ 
					double ThetaM = 0.5*(ThetaL+ThetaR);
					for ( kk = 0; kk < 4; kk++ )  
					{
						if ( !flagU[kk] )
						{
							double tpU[N];
							for ( j =0;j<N;j++ )
								tpU[j] = ThetaM*UU[kk][j]+(1.0-ThetaM)*mesh1D[k].solu[j];
							flagU0[kk] = checkallhU(tpU);
						} else 
							flagU0[kk] = true;		
					}
					
					if ( flagU0[0]+flagU0[1]+flagU0[2]+flagU0[3] == 4 )
						ThetaL = ThetaM;
					else
						ThetaR = ThetaM;
				} 
				
				//std:: cout << "theta is: " << ThetaL << std::endl;
				
				printf("theta is: %16.16f\n",ThetaL);
				
				
				for ( i = 0; i < N; i ++ )
				{
					mesh1D[k].UGL1[i] = ThetaL*mesh1D[k].UGL1[i] + (1.-ThetaL)*mesh1D[k].solu[i];
					mesh1D[k].UGL2[i] = ThetaL*mesh1D[k].UGL2[i] + (1.-ThetaL)*mesh1D[k].solu[i];
					mesh1D[k].UGL3[i] = ThetaL*mesh1D[k].UGL3[i] + (1.-ThetaL)*mesh1D[k].solu[i];	
					mesh1D[k].UGL4[i] = ThetaL*mesh1D[k].UGL4[i] + (1.-ThetaL)*mesh1D[k].solu[i];			 
				}								
				
			}
			
		}
		
		Idel2Ddouble(4,UU);
		
	
		
		// boundary condition 
		/*
		for ( i = 0; i < N; i++ )
			for ( k = 0; k < gpn ; k++ )
			{
				mesh1D[k].UGL1[i] = 	mesh1D[gpn].UGL1[i] ;
				mesh1D[k].UGL2[i] = 	mesh1D[gpn].UGL2[i] ;
				mesh1D[k].UGL3[i] = 	mesh1D[gpn].UGL3[i] ;
				mesh1D[k].UGL4[i] = 	mesh1D[gpn].UGL4[i] ;
				
				mesh1D[Ncell-k-1].UGL1[i] = 	mesh1D[Ncell-gpn-1].UGL1[i] ;
				mesh1D[Ncell-k-1].UGL2[i] = 	mesh1D[Ncell-gpn-1].UGL2[i] ;
				mesh1D[Ncell-k-1].UGL3[i] = 	mesh1D[Ncell-gpn-1].UGL3[i] ;
				mesh1D[Ncell-k-1].UGL4[i] = 	mesh1D[Ncell-gpn-1].UGL4[i] ;
			}
		*/


		for ( i = 0; i < N; i++ )
			for ( k = 0; k < gpn ; k++ )
			{
				mesh1D[k].UGL1[i] = 	mesh1D[k+ncell].UGL1[i] ;
				mesh1D[k].UGL2[i] = 	mesh1D[k+ncell].UGL2[i] ;
				mesh1D[k].UGL3[i] = 	mesh1D[k+ncell].UGL3[i] ;
				mesh1D[k].UGL4[i] = 	mesh1D[k+ncell].UGL4[i] ;
				
				mesh1D[Ncell-k-1].UGL1[i] = 	mesh1D[Ncell-k-1-ncell].UGL1[i] ;
				mesh1D[Ncell-k-1].UGL2[i] = 	mesh1D[Ncell-k-1-ncell].UGL2[i] ;
				mesh1D[Ncell-k-1].UGL3[i] = 	mesh1D[Ncell-k-1-ncell].UGL3[i] ;
				mesh1D[Ncell-k-1].UGL4[i] = 	mesh1D[Ncell-k-1-ncell].UGL4[i] ;
			}

		
		// approximate LLF matrix
		for ( k=2; k<(Ncell-3); k++ )
		{
			double temp = calallAUM ( mesh1D[k].UGL4, mesh1D[k+1].UGL1, allAUM0, allAUM1 ) ;
			temp = temp * 1.1;
			
			for ( kk = 0; kk < N; kk ++ )
			{
				
				for ( i=0; i<N; i++ )
				{
					knowb[i] = allAUM1[i][kk];
					for ( j=0; j<N; j++ )
						AA[i][j] = allAUM0[i][j];
				}
				SolverAxb(AA, knowb, unknowx);
				for ( i=0; i<N; i++ )
					allAUM[i][kk] = unknowx[i]; 
				
			}
			
			for ( i = 0; i < N; i++ )
			{
				face_AP[k][i] = 0.0;
				face_AN[k][i] = 0.0;
				for ( j = 0; j < N; j++ )
				{
					double temp2 = temp*(i==j) ;
					face_AP[k][i] += ( temp2 + allAUM[i][j] ) * ( mesh1D[k+1].UGL1[j] - mesh1D[k].UGL4[j] );
					face_AN[k][i] += ( -temp2 + allAUM[i][j] ) * ( mesh1D[k+1].UGL1[j] - mesh1D[k].UGL4[j] );
				}
				face_AP[k][i] *= 0.5;
				face_AN[k][i] *= 0.5;	
			}
			
		}
		

		
		// AUdP


		
		for ( k=3; k<(Ncell-3); k++ )
		{
			double temp = calallAU( mesh1D[k].UGL1, allAU1 , 1 );
			temp = calallAU( mesh1D[k].UGL2, allAU2 , 1 );
			temp = calallAU( mesh1D[k].UGL3, allAU3 , 1 );
			temp = calallAU( mesh1D[k].UGL4, allAU4 , 1 );	
			
			for ( i=0; i<N; i++ )
			{
	            double tp = (GLp1-GLp2)*(GLp1-GLp3)*(GLp1-GLp4);
	            dP1[i] =  (  (GLp1-GLp3)*(GLp1-GLp4) + (GLp1-GLp2)*(GLp1-GLp4) + (GLp1-GLp2)*(GLp1-GLp3) ) / tp * mesh1D[k].UGL1[i];
	            dP2[i] =  (  (GLp2-GLp3)*(GLp2-GLp4) ) / tp * mesh1D[k].UGL1[i];
	            dP3[i] =  (  (GLp3-GLp2)*(GLp3-GLp4) ) / tp * mesh1D[k].UGL1[i];
	            dP4[i] =  (  (GLp4-GLp2)*(GLp4-GLp3) ) / tp * mesh1D[k].UGL1[i];
	            
	            tp = (GLp2-GLp1)*(GLp2-GLp3)*(GLp2-GLp4);
	            dP1[i] = dP1[i] + (  (GLp1-GLp3)*(GLp1-GLp4) ) / tp * mesh1D[k].UGL2[i];
	            dP2[i] = dP2[i] + (  (GLp2-GLp3)*(GLp2-GLp4) + (GLp2-GLp1)*(GLp2-GLp4) + (GLp2-GLp1)*(GLp2-GLp3) ) / tp * mesh1D[k].UGL2[i];
	            dP3[i] = dP3[i] + (  (GLp3-GLp1)*(GLp3-GLp4) ) / tp * mesh1D[k].UGL2[i];
	            dP4[i] = dP4[i] + (  (GLp4-GLp1)*(GLp4-GLp3) ) / tp * mesh1D[k].UGL2[i];
	            
	            tp = (GLp3-GLp1)*(GLp3-GLp2)*(GLp3-GLp4);
	            dP1[i] = dP1[i] + (  (GLp1-GLp2)*(GLp1-GLp4) ) / tp * mesh1D[k].UGL3[i];
	            dP2[i] = dP2[i] + (  (GLp2-GLp1)*(GLp2-GLp4) ) / tp * mesh1D[k].UGL3[i];
	            dP3[i] = dP3[i] + (  (GLp3-GLp1)*(GLp3-GLp4) + (GLp3-GLp1)*(GLp3-GLp2) + (GLp3-GLp2)*(GLp3-GLp4)  ) / tp * mesh1D[k].UGL3[i];
	            dP4[i] = dP4[i] + (  (GLp4-GLp1)*(GLp4-GLp2) ) / tp * mesh1D[k].UGL3[i];
	            
	            tp = (GLp4-GLp1)*(GLp4-GLp2)*(GLp4-GLp3);
	            dP1[i] = dP1[i] + (  (GLp1-GLp2)*(GLp1-GLp3) ) / tp * mesh1D[k].UGL4[i];
	            dP2[i] = dP2[i] + (  (GLp2-GLp1)*(GLp2-GLp3) ) / tp * mesh1D[k].UGL4[i];
	            dP3[i] = dP3[i] + (  (GLp3-GLp1)*(GLp3-GLp2)) / tp * mesh1D[k].UGL4[i];
	            dP4[i] = dP4[i] + (  (GLp4-GLp1)*(GLp4-GLp2) + (GLp4-GLp1)*(GLp4-GLp3) + (GLp4-GLp2)*(GLp4-GLp3) ) / tp * mesh1D[k].UGL4[i];
	            				
			}
			
			for ( i = 0; i < N; i++ )
			{
				mesh1D[k].AUdP[i] = 0.;
				for ( j = 0; j < N; j++ )
					mesh1D[k].AUdP[i] += allAU1[i][j]*dP1[j]+5.*allAU2[i][j]*dP2[j]+5.*allAU3[i][j]*dP3[j]+allAU4[i][j]*dP4[j];
				mesh1D[k].AUdP[i] = mesh1D[k].AUdP[i] / 12.0;	
			}			
					
		}
		

		
		
        // evolution
		if ( RK == 1 )	
		  for ( k=0; k<Ncell; k++ )
		    for ( i = 0; i< N; i++ )
		      mesh1D[k].oldU[i] = mesh1D[k].solu[i];


        for ( k=3; k<(Ncell-3); k++ )
        {
        	for ( i = 0; i< N; i++ )
        	{
        		double temp = mesh1D[k].solu[i] - 
	                    lamx*( mesh1D[k].AUdP[i] + face_AN[k][i]  + face_AP[k-1][i] );
	                    
        		if ( RK == 1 ){
        			mesh1D[k].oldU[i] = mesh1D[k].solu[i];
				mesh1D[k].solu[i] = temp;
        		} else if ( RK == 2 )
        			mesh1D[k].solu[i] = 0.75*mesh1D[k].oldU[i] + 0.25* temp;
        		else
        			mesh1D[k].solu[i] = (1.0/3)*mesh1D[k].oldU[i] + (2.0/3)* temp;
			if (mesh1D[k].oldU[i] != mesh1D[k].oldU[i] ) {
			  std::cout << mesh1D[k].oldU[i] << std::endl;
			  exit(1);
			}
        	}
        }
  
		// boundary condition 
		for ( i = 0; i < N; i++ )
			for ( k = 0; k < gpn ; k++ )
			{
				mesh1D[k].solu[i] = 	mesh1D[k+ncell].solu[i] ;
	
				mesh1D[Ncell-k-1].solu[i] = 	mesh1D[Ncell-k-1-ncell].solu[i] ;
			}		
		

	}

	return 0;
	 
	
	
}


double  calallAU( double *hU, double  allAU[N][N] , bool flag)
{
	double allhA0[N][N], allhA1[N][N];
	double maxeig;
	maxeig = calallhA( hU, allhA0, allhA1 );
	
	if ( !flag )
		return maxeig;

	for ( int kk = 0; kk < N; kk ++ )
	{
		double AA[N][N],knowb[N],unknowx[N];
		for ( int i=0; i<N; i++ )
		{
			knowb[i] = allhA1[i][kk];
			for ( int j=0; j<N; j++ )
				AA[i][j] = allhA0[i][j];
		}
		SolverAxb(AA, knowb, unknowx);
		for ( int i=0; i<N; i++ )
			allAU[i][kk] = unknowx[i]; 
				
	}
	
	return maxeig;
}




double calallAUM (double *hUL, double *hUR, double  allAUM0[N][N], double  allAUM1[N][N])
{
	const double weight[5]={0.5688888888888889,
							0.4786286704993665,
							0.4786286704993665,
							0.2369268850561891,
							0.2369268850561891};
	const double gs_p[5] = {0.0000000000000000,
							-0.5384693101056831,
							0.5384693101056831,
							-0.9061798459386640,
							0.9061798459386640};
	const int ngp = 5;
	int k,i,j;
	double allhA0[N][N], allhA1[N][N];
	double hU[N];
	double maxeig = 0.;	
	
	for ( i=0;i<N;i++ )
		for (j=0;j<N;j++) 
			allAUM0[i][j] = allAUM1[i][j] = 0.;
	
	for ( k=0; k < ngp; k ++ )
	{
		for ( i=0;i<N;i++ )
			hU[i] = hUL[i] + 0.5*(gs_p[k]+1.) * ( hUR[i] - hUL[i] );
		double temp = calallhA ( hU, allhA0, allhA1 );
		
		if ( temp > maxeig )
			maxeig = temp;
			
		for ( i=0;i<N;i++ )
			for (j=0;j<N;j++) {
				allAUM0[i][j] += 0.5*weight[k] * allhA0[i][j];
				allAUM1[i][j] += 0.5*weight[k] * allhA1[i][j];
			}
	}
	
	return maxeig;
}







double calallhA ( double *hU, double  allhA0[N][N], double  allhA1[N][N] )
{
	double maxeig = 0.;
	int k,i,j;
	double hA0[N][N], hA1[N][N];

	for ( i=0;i<N;i++ )
		for (j=0;j<N;j++) 
			allhA0[i][j] = allhA1[i][j] = 0.;

	for (k=0; k<20; k++)
	{
		double temp = calhAU( hU, global_p[k], hA0, hA1 ) ;
		if ( temp > maxeig )
			maxeig = temp;
		for ( i=0;i<N;i++ )
			for (j=0;j<N;j++) {
				allhA0[i][j] += global_w[k] * hA0[i][j];
				allhA1[i][j] += global_w[k] * hA1[i][j];
			}			
	}
	
	return maxeig;
} 







double calhAU( double *hU, double xi, double  hA0[N][N], double  hA1[N][N] )
{
	double allPhi[M+1];
	int k,i;
	for ( k = 0; k<=M; k++ )
		allPhi[k] = Phix(k,xi) * sqrt(k+0.5);
	double uM[DOS];
	for (i =0; i<DOS; i++) {
		uM[i] = 0.;
		for ( k = 0; k<=M; k++ )
			uM[i] += allPhi[k] * hU[k*DOS+i];
	}	
	
	double A0[DOS][DOS], A1[DOS][DOS];
	double maxeig=0;
	calAU( uM, A0, A1 , maxeig );
	
	int kk,ii;
	
	for (k=0; k<=M; k++)
		for (i=0; i<=M; i++)
			for (kk=0; kk<DOS; kk++)
				for (ii=0; ii<DOS; ii++)  {
					hA0[k*DOS+kk][i*DOS+ii] = allPhi[k]*allPhi[i]*A0[kk][ii]; 
					hA1[k*DOS+kk][i*DOS+ii] = allPhi[k]*allPhi[i]*A1[kk][ii]; 
				}
	return maxeig; 
}







void calAU( double U[DOS], double  AU0[DOS][DOS], double  AU1[DOS][DOS] , double & maxeig) 
{
	const double rho = U[0];
	const double u = U[1]/U[0];
	const double v = U[2]/U[0];
	const double p = (GAMMA-1.)*(U[3]-0.5*rho*(u*u+v*v));

	if ( rho < 0 || p < 0 )
	  {
	    std::cout << "negative rho or p:" << rho << " " << p <<"\n";exit(1);
	  }



	const double c2 = GAMMA*p/rho;
	const double c = sqrt(c2);
	const double B1 = (GAMMA-1.)/c2;
	const double B2 = B1*(u*u+v*v)/2.;
	
	double Lx[DOS][DOS];
	Lx[0][0] = 0.5*(B2+u/c);
	Lx[0][1] = -0.5*(B1*u+1./c);
	Lx[0][2] = -0.5*B1*v;
	Lx[0][3] = 0.5*B1;
	
	Lx[1][0] = v;
	Lx[1][1] = 0.;
	Lx[1][2] = -1.;
	Lx[1][3] = 0.;
	
	Lx[2][0] = 1.-B2;
	Lx[2][1] = B1*u;
	Lx[2][2] = B1*v;
	Lx[2][3] = -B1;
	
	Lx[3][0] = 0.5*(B2-u/c);
	Lx[3][1] = -0.5*(B1*u-1./c);
	Lx[3][2] = -0.5*B1*v;
	Lx[3][3] = 0.5*B1;
	
	maxeig = fabs(u)+c;
	
	double Lam[DOS];
	
	Lam[0] = u - c;
	Lam[1] = u;
	Lam[2] = u;
	Lam[3] = u + c;
	
	double tp = 2.0*U[0]*U[3]-(U[1]*U[1]+U[2]*U[2]);
	tp = tp*tp;
	
	int i,j,k;
	
	for ( i=0; i < DOS; i++ )
		for ( j=0; j<DOS; j++ )
		{
			AU0[i][j] = 0.0;
			for ( k=0; k<DOS; k++ )
				AU0[i][j] += Lx[k][i] * Lx[k][j];
			AU0[i][j] *= tp;

			AU1[i][j] = 0.0;
			for ( k=0; k<DOS; k++ )
				AU1[i][j] += Lam[k] * Lx[k][i] * Lx[k][j];
			AU1[i][j] *= tp;			
		}
	
}




















bool checkpositive( double U[DOS] )
{
	const double _eps = 1.0e-12;
	if ( U[0] < _eps || (GAMMA-1.)*( U[3]-0.5*(U[1]*U[1]+U[2]*U[2])/U[0] ) < _eps )
		return 0;
	else
		return 1; 	
}



bool checkhU( double *hU, double xi )
{
	double allPhi[M+1];
	int k,i;
	for ( k = 0; k<=M; k++ )
		allPhi[k] = Phix(k,xi) * sqrt(k+0.5);
	double uM[DOS];
	for (i =0; i<DOS; i++) {
		uM[i] = 0.;
		for ( k = 0; k<=M; k++ )
			uM[i] += allPhi[k] * hU[k*DOS+i];
	}	
		
	return checkpositive(uM);
}

bool checkallhU( double *hU )
{ 
	for ( int k=0; k<20; k++ ) 
		if ( ! checkhU ( hU, global_p[k] ) )
			return 0;
	return 1;	
}




double Phix(int n, double x)
{ 
	if ( n == 0 )
    	return 1.;
	else if ( n == 1 )
		return x;
	else
		return (2.-1.0/n)*x*Phix(n-1,x)-(1.0-1.0/n)*Phix(n-2,x);
} 





/** 
 * 
 * 
 * @param v5 cell average u_{i-2}
 * @param v4 
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v2 
 * @param v1 
 * 
 * @return  point value u_{i-1/2}
 */
double WENOGL41(double v5,double v4,double v3,double v2,double v1){
     double s1,s2,s3;
     double sum;
     double tmp1,tmp2;
     double eps=1e-14;
     double w1,w2,w3;
     
     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s1=0.1/(s1*s1);

     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s2=0.6/(s2*s2);

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s3=0.3/(s3*s3);

     sum=s1+s2+s3;
     w1=s1/sum;
     w2=s2/sum;
     w3=s3/sum;

     return (w1*(v1/3-7*v2/6+11*v3/6)+w2*(-v2/6+5*v3/6+v4/3)+w3*(v3/3+5*v4/6-v5/6));
}

/** 
 * 
 * 
 * @param v5 cell average u_{i-2}
 * @param v4 
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v2 
 * @param v1 
 * 
 * @return  point value u_{i-sqrt(5)/10}
 */
double WENOGL42(double v5,double v4,double v3,double v2,double v1){
     double s1,s2,s3;
     double sum;
     double tmp1,tmp2;
     double w1,w2,w3;
     double value;
     double sq5=sqrt(5.0);
     double eps=1e-14;
     
     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;

     s1=(91+9*sq5)/(440*s1*s1);

     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;

     s2=129./(220*s2*s2);

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s3=(91-9*sq5)/(440*s3*s3);

     sum=s1+s2+s3;
     w1=s1/sum;
     w2=s2/sum;
     w3=s3/sum;

     value=w1*((-1./60+0.05*sq5)*v1+(1./30-sq5*0.2)*v2+(59./60+0.15*sq5)*v3)+w2*((-1./60-0.05*sq5)*v2+31./30*v3+(-1.0/60+0.05*sq5)*v4)+w3*((59./60-0.15*sq5)*v3+(1./30+0.2*sq5)*v4+(-1./60-0.05*sq5)*v5);
     return value;
}

/** 
 * 
 * 
 * @param v1 cell average u_{i-2}
 * @param v2 
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v4 
 * @param v5 
 * 
 * @return  point value u_{i+sqrt(5)/10}
 */
double WENOGL43(double v1,double v2,double v3,double v4,double v5){
     double s1,s2,s3;
     double sum;
     double tmp1,tmp2;
     double w1,w2,w3;
     double value;
     double sq5=sqrt(5.0);
     double eps=1e-14;
     
     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;

     s1=(91+9*sq5)/(440*s1*s1);

     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;

     s2=129./(220*s2*s2);

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s3=(91-9*sq5)/(440*s3*s3);

     sum=s1+s2+s3;
     w1=s1/sum;
     w2=s2/sum;
     w3=s3/sum;

     value=w1*((-1./60+0.05*sq5)*v1+(1./30-sq5*0.2)*v2+(59./60+0.15*sq5)*v3)+w2*((-1./60-0.05*sq5)*v2+31./30*v3+(-1.0/60+0.05*sq5)*v4)+w3*((59./60-0.15*sq5)*v3+(1./30+0.2*sq5)*v4+(-1./60-0.05*sq5)*v5);
     return value;
}


/** 
 * 
 * 
 * @param v1 cell average u_{i-2}
 * @param v2 
 * @param v3  cell average u_i , cell i: [i-1/2,i+1/2]
 * @param v4 
 * @param v5 
 * 
 * @return  point value u_{i+1/2}
 */
double WENOGL44(double v1,double v2,double v3,double v4,double v5){
     double s1,s2,s3;
     double sum;
     double tmp1,tmp2;
     double eps=1e-14;
     double w1,w2,w3;
     
     tmp1=v1-2*v2+v3;
     tmp2=v1-4*v2+3*v3;
     s1=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s1=0.1/(s1*s1);

     tmp1=v2-2*v3+v4;
     tmp2=v2-v4;
     s2=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s2=0.6/(s2*s2);

     tmp1=v3-2*v4+v5;
     tmp2=3*v3-4*v4+v5;
     s3=13.0/12*tmp1*tmp1+0.25*tmp2*tmp2+eps;
     s3=0.3/(s3*s3);

     sum=s1+s2+s3;
     w1=s1/sum;
     w2=s2/sum;
     w3=s3/sum;

     return (w1*(v1/3-7*v2/6+11*v3/6)+w2*(-v2/6+5*v3/6+v4/3)+w3*(v3/3+5*v4/6-v5/6));
}

/** solve Ax=b by Gauss's method **/
void SolverAxb(double AA[N][N], double b[N], double  x[N]) 
{
    int i, j, d, row;
    double temp;
    double known_items;
    double l[N][N],a[N][N];

	
	for ( i = 0; i < N; i++ )	
		for (j = 0; j < N ;j++)
			a[i][j] = AA[i][j];		

    for (d=0; d<N-1;d++)
    {
		row=d;
		for (i=d+1; i<N; i++) {
			if(fabs(a[i][d])>fabs(a[row][d]))
				row=i;
		}
		if (row!=d)	{
			for (j=d;j<N; j++) {
				temp=a[row][j];
				a[row][j]=a[d][j];
				a[d][j]=temp;
			}
			temp=b[row];
			b[row]=b[d];
			b[d]=temp;
		}

		if ( fabs(a[d][d]) < 1.0e-10 ) {
			printf("sigular in solver Ax=b! : %e\n",a[d][d]);
			fflush(0);exit(1);
		}

		for (i=d+1; i<N; i++)
		{		   			
			l[i][d]=-a[i][d]/a[d][d];
			for (j=d;j<N; j++)
				a[i][j]=a[i][j]+a[d][j]*l[i][d];
			b[i]=b[i]+b[d]*l[i][d];
		}
    }

    for (i=N-1; i>=0; i--)
    {
		known_items=0.0;
		for (j=1; j<N-i; j++)
			known_items=known_items+a[i][i+j]*x[i+j];

		if ( fabs(a[i][i]) < 1.0e-10 ) {
			printf("sigular in solver Ax=b! : %e\n",a[i][i]);
			fflush(0);exit(1);
		}

		x[i]=(b[i]-known_items)/a[i][i];
    }

}

void uexact(double x, double y, double  xi, double U[DOS])
{
  double rho,u,v,p,T,_eps,xbar,ybar,r,S;

  xbar = x - 5.;
  ybar = y - 5.;
  r = sqrt(xbar*xbar+ybar*ybar);
  _eps = 5. + xi * 1.;
  
  const double PI = 4. * atan(1.);
  

  u = 1. + _eps*exp(0.5*(1.-r*r))*(-ybar)/(2.*PI);
  v = 1. + _eps*exp(0.5*(1.-r*r))*xbar/(2.*PI);
  T = 1. - (GAMMA-1.)*_eps*_eps*exp(1.-r*r)/(8.*GAMMA*PI*PI);
  S = 1.;
  rho = pow(T/S,1./(GAMMA-1.));
  p = rho * T;
  
  U[0] = rho;
  U[1] = rho*u;
  U[2] = rho*v;
  U[3] = p/(GAMMA-1.) + 0.5*rho*(u*u+v*v);
}


void initial ( double x, double y, double hU[N] )
{
	int i,j,k;
	for ( i=0; i<=M; i++ )
	{
		for ( j=0;j<DOS;j++ )
		{
			hU[i*DOS+j] = 0.0;
			for ( k=0; k<20;k++ )
			{
				double U[DOS];
				uexact(x,y,global_p[k],U);
				hU[i*DOS+j] += global_w[k]*Phix(i,global_p[k])*sqrt(i+0.5)*U[j];
			}
		}
	}
	
} 





int main()
{
	double xa = 0., xb = 10.0;//8.0*atan(1.0);
	//const int ncell = 50;
	const int Ncell = ncell + 2*gpn;
	double dx = (xb - xa)/ncell;
	double ntime = 0;
	double tstop = 0.01;
	int i,j,k,ii,kk;
	
	double xmid[Ncell];
	MeshCell *mesh;
	mesh = new MeshCell [Ncell];
	
	double ***mesh2Dsolu;
	mesh2Dsolu = Iget3Dmatrix(Ncell,Ncell,N);
			
	
	for ( k = 0; k < Ncell; k++ )
		xmid[k] = xa - 2.5*dx + k*dx;

	double tp_g[4],tp_w[4];
	const double _eps = 2.0e-16;
	tp_g[0] = -0.5+2*_eps/dx;
	tp_g[1] = -0.1*sqrt(5.0);
	tp_g[2] = 0.1*sqrt(5.0);
	tp_g[3] = 0.5-2*_eps/dx;
	tp_w[0] = 1./12.;
	tp_w[1] = 5./12.;
	tp_w[2] = 5./12.;
	tp_w[3] = 1./12.;


	
	// initial
	for ( k = 0; k < Ncell; k++ )
	{
	  std::cout << k << std::endl;
	  for ( i = 0; i < Ncell; i++)
	    {	
	      for ( j = 0; j < N; j++ )
		mesh2Dsolu[k][i][j] = 0.;

	      for ( kk = 0; kk < 4; kk++ )
		{
		  for ( ii = 0; ii < 4; ii++ )
		    {
		      initial( xmid[k]+tp_g[kk]*dx,  xmid[i]+tp_g[ii]*dx, mesh[0].UGL1 );
		      for ( j = 0; j < N; j++ )
			mesh2Dsolu[k][i][j] += tp_w[kk]*tp_w[ii]*mesh[0].UGL1[j];

		    }

		}
		
	    }	  
	}
	
	
	std::cout << "initial solution obtained\n";
	
	
	// updata
	double allAU[N][N];
	double tp_U[N];
	int IT = 0;
	while ( ntime < tstop )
	{
		double maxeig = 0;
		for ( k = 0; k < Ncell; k++ )
		{
		  for ( i = 0; i < Ncell; i++ )
		    {
			double temp = calallAU( mesh2Dsolu[k][i], allAU, 0 );
			if (temp > maxeig)
				maxeig = temp;

			for ( j = 0; j < N; j++ )
			  {
			    if ( j%DOS == 1 )
			      tp_U[j] = mesh2Dsolu[k][i][j+1];
			    else if ( j%DOS == 2 )
			      tp_U[j] = mesh2Dsolu[k][i][j-1];
			    else
			      tp_U[j] = mesh2Dsolu[k][i][j];
			  }

			temp = calallAU( tp_U, allAU, 0 );
			if (temp > maxeig)
				maxeig = temp;
		    }
		}
		
						
		double dt = CFL * dx / maxeig;

		dt = std::min(dt, pow(dx,5.0/3.0));

		if ( ntime + dt > tstop )
			dt = tstop - ntime;

		
		std::cout << "time step obtained\n";
	       
			
		//double lamx = dt / dx;
		bool halfdt = 1;



		//time step for 3rd order Strang spliting
		const double taustep1 = 2.0 / ( 5.0-sqrt(13.0) + sqrt(2.0*(1.0+sqrt(13.0))) );
		const double taustep2 = ( 7.0 + sqrt(13.0) - sqrt(2.0*(1.0+sqrt(13.0))) ) / 12.0;
		const double taustep3 = taustep1 * taustep1 / ( taustep2 - taustep1 );
		const double taustep4 = 1.0 - ( taustep1 + taustep2 + taustep3 );


		if ( IT % 2 == 0 )
		  {
		   
		    std::cout << "strang step 1" << std::endl; 

		    // E_x tau_1
		    for ( i = 0; i < Ncell; i++ )
		      {
			
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh[k].solu[j] = mesh2Dsolu[k][i][j];

			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep1*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh2Dsolu[k][i][j] = mesh[k].solu[j];		    
		      }


		    // E_y tau_1+tau_2

		    std::cout << "strang step 2" << std::endl; 

		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j+1];
			    else if ( j%DOS == 2 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j-1];
			    else
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j];
			  }

			if ( halfdt = RKUpdate1D (  mesh,  Ncell, (taustep1+taustep2)*dt/dx ) )
			  exit(1);
		    

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh2Dsolu[i][k][j+1] = mesh[k].solu[j];
			    else if ( j%DOS == 2 )
			      mesh2Dsolu[i][k][j-1] = mesh[k].solu[j];
			    else
			      mesh2Dsolu[i][k][j] = mesh[k].solu[j];
			  }
		    
		      }

		    std::cout << "strang step 3" << std::endl; 


		    // E_x 
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh[k].solu[j] = mesh2Dsolu[k][i][j];
			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep2*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh2Dsolu[k][i][j] = mesh[k].solu[j];		    
		      }


		    std::cout << "strang step 4" << std::endl; 

		    // E_y 
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j+1];
			    else if ( j%DOS == 2 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j-1];
			    else
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j];
			  }
		    

			if ( halfdt = RKUpdate1D (  mesh,  Ncell, (taustep3)*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh2Dsolu[i][k][j+1] = mesh[k].solu[j];
			    else if ( j%DOS == 2 )
			      mesh2Dsolu[i][k][j-1] = mesh[k].solu[j];
			    else
			      mesh2Dsolu[i][k][j] = mesh[k].solu[j];
			  }		    
		      }




		    std::cout << "strang step 5" << std::endl; 
		    // E_x
		    for ( i = 0; i < Ncell; i++ )
		      {
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh[k].solu[j] = mesh2Dsolu[k][i][j];
			if ( halfdt = RKUpdate1D (  mesh,  Ncell, (taustep3+taustep4)*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh2Dsolu[k][i][j] = mesh[k].solu[j];		    
		      }


		    std::cout << "strang step 6" << std::endl; 

		    // E_y
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j+1];
			    else if ( j%DOS == 2 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j-1];
			    else
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j];
			  }

			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep4*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh2Dsolu[i][k][j+1] = mesh[k].solu[j];
			    else if ( j%DOS == 2 )
			      mesh2Dsolu[i][k][j-1] = mesh[k].solu[j];
			    else
			      mesh2Dsolu[i][k][j] = mesh[k].solu[j];
			  }
	    
		      }


		  } else { //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		  std::cout << "strang step 1" << std::endl; 

		    // E_y tau_1
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j+1];
			    else if ( j%DOS == 2 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j-1];
			    else
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j];
			  }


			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep1*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh2Dsolu[i][k][j+1] = mesh[k].solu[j];
			    else if ( j%DOS == 2 )
			      mesh2Dsolu[i][k][j-1] = mesh[k].solu[j];
			    else
			      mesh2Dsolu[i][k][j] = mesh[k].solu[j];
			  }		    
		      }


		    
		    std::cout << "strang step 2" << std::endl; 

		    // E_x tau_1+tau_2
		    for ( i = 0; i < Ncell; i++ )
		      {
	
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh[k].solu[j] = mesh2Dsolu[k][i][j];
			if ( halfdt = RKUpdate1D (  mesh,  Ncell, (taustep1+taustep2)*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh2Dsolu[k][i][j] = mesh[k].solu[j];		    
		      }



		    std::cout << "strang step 3" << std::endl; 

		    // E_y
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j+1];
			    else if ( j%DOS == 2 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j-1];
			    else
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j];
			  }

			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep2*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh2Dsolu[i][k][j+1] = mesh[k].solu[j];
			    else if ( j%DOS == 2 )
			      mesh2Dsolu[i][k][j-1] = mesh[k].solu[j];
			    else
			      mesh2Dsolu[i][k][j] = mesh[k].solu[j];
			  }
		      }


		    std::cout << "strang step 4" << std::endl; 
		    // E_x
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh[k].solu[j] = mesh2Dsolu[k][i][j];
			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep3*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh2Dsolu[k][i][j] = mesh[k].solu[j];		    
		      }




		    std::cout << "strang step 5" << std::endl; 
		    // E_y
		    for ( i = 0; i < Ncell; i++ )
		      {


			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j+1];
			    else if ( j%DOS == 2 )
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j-1];
			    else
			      mesh[k].solu[j] = mesh2Dsolu[i][k][j];
			  }

			if ( halfdt = RKUpdate1D (  mesh,  Ncell, (taustep3+taustep4)*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ ) {
			    if ( j%DOS == 1 )
			      mesh2Dsolu[i][k][j+1] = mesh[k].solu[j];
			    else if ( j%DOS == 2 )
			      mesh2Dsolu[i][k][j-1] = mesh[k].solu[j];
			    else
			      mesh2Dsolu[i][k][j] = mesh[k].solu[j];
			  }		    
		      }



		    std::cout << "strang step 6" << std::endl; 
		    // E_x
		    for ( i = 0; i < Ncell; i++ )
		      {

			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh[k].solu[j] = mesh2Dsolu[k][i][j];
			if ( halfdt = RKUpdate1D (  mesh,  Ncell, taustep4*dt/dx ) )
			  exit(1);
		    
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < N; j++ )
			    mesh2Dsolu[k][i][j] = mesh[k].solu[j];		    
		      }

		}









		/*
		bool halfdt = 1;
		int noofhalf = 0;
		
		do
		  {
		    if ( halfdt = RKUpdate1D (  mesh,  Ncell, lamx ) )
		      {
			dt *= 0.5;
			lamx *= 0.5;
			noofhalf++;
			std::cout << "half time step!\t" << noofhalf << std::endl;
			if ( noofhalf > 6 )
			  exit(1);
			for ( k = 0; k < Ncell; k++ )
			  for ( j = 0; j < Ncell; j++ )
			     mesh[k].solu[i] = mesh[k].oldU[i];
		      }		    
		  }while( halfdt );
		*/



		
		ntime += dt;
		++IT;
		
		std::cout << "IT:"<<IT<<",t = " << ntime << ", dt = " << dt << std::endl;
		
		
 
	}
	
	
	// output
	double Vrho;
	FILE *fp = fopen("dataE.txt","w");
	FILE *fq = fopen("dataV.txt","w");

	for ( k=gpn; k < ncell+gpn; k++ )
	  {
	    for ( i=gpn; i < ncell+gpn; i++ )
	      {
		fprintf(fp,"%16.16f ",mesh2Dsolu[k][i][0]/sqrt(2.));

		Vrho = 0.;
		for ( int kk = 1; kk<=M ; kk++ )
		  Vrho += pow( mesh2Dsolu[k][i][DOS*kk] , 2.0 );
		Vrho *= .5; 
		fprintf(fq,"%16.16f ",Vrho);
	      }
	    fprintf(fp,"\n");
	    fprintf(fq,"\n");
	  }
	
	fclose(fp);	
	fclose(fq);


	//output exact solutions

	for ( k = 0; k < Ncell; k++ )
	{
	  std::cout << k << std::endl;
	  for ( i = 0; i < Ncell; i++)
	    {	
	      for ( j = 0; j < N; j++ )
		mesh2Dsolu[k][i][j] = 0.;

	      for ( kk = 0; kk < 4; kk++ )
		{
		  for ( ii = 0; ii < 4; ii++ )
		    {
		      initial( xmid[k]+tp_g[kk]*dx - tstop,  xmid[i]+tp_g[ii]*dx - tstop, mesh[0].UGL1 );
		      for ( j = 0; j < N; j++ )
			mesh2Dsolu[k][i][j] += tp_w[kk]*tp_w[ii]*mesh[0].UGL1[j];

		    }

		}
		
	    }	  
	}



	fp = fopen("dataE0.txt","w");
	fq = fopen("dataV0.txt","w");

	for ( k=gpn; k < ncell+gpn; k++ )
	  {
	    for ( i=gpn; i < ncell+gpn; i++ )
	      {
		fprintf(fp,"%16.16f ",mesh2Dsolu[k][i][0]/sqrt(2.));

		Vrho = 0.;
		for ( int kk = 1; kk<=M ; kk++ )
		  Vrho += pow( mesh2Dsolu[k][i][DOS*kk] , 2.0 );
		Vrho *= .5; 
		fprintf(fq,"%16.16f ",Vrho);
	      }
	    fprintf(fp,"\n");
	    fprintf(fq,"\n");
	  }
	
	fclose(fp);	
	fclose(fq);





	
	
	delete [] mesh;
	Idel3Dmatrix(Ncell,Ncell,mesh2Dsolu);
	
	return 0;
}



