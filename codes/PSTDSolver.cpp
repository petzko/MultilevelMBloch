#include <iostream>
#include "./LTensor/LTensor.h"
#include <cmath>
#include <complex.h>
#include <fftw3.h>



typedef fftw_complex cdouble; 

#define DIM 1

Marray<double,3> levichevi(3,3,3);

void initialize(int, char ** );

int getLINidx(std::vector<int> idx, int dim);
std::vector<int> getVectorIdx(int dim, int LI);

Marray<cdouble,2> fast_fourier(Marray<cdouble,2>,Marray<cdouble,2> , int);

int main(int argc, char** argv){


	initialize(argc,argv); 

	/** initialize grid points 
	 * & side length for each direction as well as 
	 * calculate the grid spacing for that given geometry...
	 */


	/***
	  if you want to use the convenient tensor notation... 
	  we must have equal nr of grid points for each direction ... 

	 **/

	int grid_size = 64; 
	int grid_size_dim = 1; 
	for(int d =0 ; d< DIM; d++)
		grid_size_dim *= grid_size;
 
	int N[DIM];
	double L[DIM];
	double dx[DIM]; 
	for(int d =0; d< DIM; d++) 
	{
		N[d] = grid_size ;  
		L[d] = 1.0;
		dx[d] = L[d]/(N[d] - 1 ); 
	} 


	// allocate memory for the wavenr array and the coordinates array... 
	Marray<cdouble,2>* k_v = new Marray<cdouble,2>(grid_size_dim,DIM); 
	Marray<cdouble,2>* r_v = new Marray<cdouble,2>(grid_size_dim,DIM); 
	Marray<cdouble,2>* E_v = new Marray<cdouble,2>(grid_size_dim,DIM); 
	Marray<cdouble,2>* H_v = new Marray<cdouble,2>(grid_size_dim,DIM); 


	// initialize the grid points, the wave vectors and the EM field vectors...  
	// for transversal of the whole computational domain we will need dim*N^dim iterations... 

	for(int d= 0 ; d< DIM; d++){	
		for(int n = 0 ; n< N[d]; n++)
		{
			(*k_v)(n,d) = 2.0*M_PI*n/N[d];
			(*r_v)(n,d) = dx[d]*n; 
			(*H_v)(n,d) = 0; 
			(*E_v)(n,d) = 0;

		}
	}

	// time variables  	
	double t = 0 ; 
	// timestep size
	double dt =10E-4;
	double T =1; // end after 1 sec of simulation 
	
	// dummy indices.. for Einstein notation ...
	

	Index<'i'> i_dx; 
	Index<'n'> n_dx;
	double mu = 1.0; double eps = 1.0; 
 	Marray<cdouble,2> DxE(grid_size,DIM);
	Marray<cdouble,2> DxH(grid_size,DIM);  
	
	while(t < T) 
	{
		DxE = fast_fourier((*E_v),(*k_v),grid_size); 
		(*H_v)(n_dx,i_dx) = (*H_v)(n_dx,i_dx) - dt/mu*(DxE(n_dx,i_dx));
			
		DxH = fast_fourier((*H_v),(*k_v),grid_size);
 		(*E_v)(n_dx,i_dx) = (*E_v)(n_dx,i_dx) +dt/eps*(DxH(n_dx,i_dx));  

		t+=dt; 
//		std::cout<<"H @ "<<  t << " = " << (*H_v) <<" \n";
	}



	delete k_v ;
	delete r_v ; 
	delete E_v; 
	delete H_v; 
}

// compute, using the fast fourier transform... 
// the expression NABLA x X_v , approximating the spatial derivative using the fourier differentiation theorem 

Marray<cdouble,1> fast_fourier(Marray<cdouble,2> X_v,Marray<cdouble,2> k_v, int N){
	
// initialize memory for the result... 	
	std::cout<<"FFTT... \n"; 
	Marray<cdouble,2> res(N,DIM); 
	fftw_complex *in_f, *out_f,*in_k;
	fftw_complex *in_b, *out_b; 

	fftw_plan p_f; fftw_plan p_b;
	// prepare the plan for the forward direction fft..
 
	in_f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out_f = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p_f = fftw_plan_dft_1d(N, in_f, out_f, FFTW_FORWARD, FFTW_ESTIMATE);
	
	in_k = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	
	in_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	out_b = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
	p_b = fftw_plan_dft_1d(N, in_b, out_b, FFTW_BACKWARD, FFTW_ESTIMATE);
		
	for(int n = 0; n < N; n ++)
	{
		for(int d = 0; d< DIM; d++)
		{
			in_f[n*DIM + d] = X_v(n,d);
			in_k[n*DIM + d] = k_v(n,d);
		}
	}
	/** 
		1) 
		FFT(in_F)
	*/
	fftw_execute(p_f); /* repeat as needed */


	for(int n = 0; n < N; n ++)
	{
                for(int d = 0; d< DIM; d++)
		{
                        in_b[n*DIM + d] = I*k_v(n,d)*out_f[n*DIM+d]; 
                        
                }
        }

	fftw_execute(p_b);

	for(int n = 0; n < N; n ++)
        {
                for(int d = 0; d< DIM; d++)
                {
			res(n,d) = out_b[n*DIM+d];
		}
	}

	// free allocated memory.... 
	fftw_destroy_plan(p_f);
	fftw_destroy_plan(p_b); 
	fftw_free(in_k); 
	fftw_free(in_f); fftw_free(out_f); fftw_free(in_b); fftw_free(out_b) ; 


	
	/***
		using fftw make fourier transform
	 	1) 
			F(X_v); 
		2)
i			then pairwise multiplicaiton 
			i * k_v * F(X_v);
		3)	
			then ifft
			F^-1(i*k_v*F(X_v) ; 
		4) renormalize the result (fftw does an unnormalized fft) 
			so divide by grid_size... 	
		return the result... 			
	***/

	return res; 

}

void getOffsets(int dim, std::vector<int> N, std::vector<int>& offsets)
{
	offsets.resize(dim,1);
	int tmp = 1; 	
	for(int j = dim-1; j >=0 ; j--)
		offsets[j] = tmp;
		tmp*=N[dim-1-j];

	
} // dim = 3 ->  offset[2] = 1, offset[1] = N[0], offset[0] = N[0]*N[1]


int getLINidx(std::vector<int> idx,std::vector<int> offsets, int dim){
	int I = 0;
	for (int p = 0; j < dim; p++)
		I+= idx[p]*offsets[p];
	
	return I; 
}


std::vector<int> getVectorIdx(int dim, int LI,std::vector<int> offsets){

	std::vector<int> idx (dim,0); 
	for(int p = 0; p<dim, p++)
	{
		idx[p] = LI / offsets[p];
		LI = LI % offsets[p];
	}	

}


void initialize(int argc, char** argv){


        for(int i = 0 ; i < 3; i++)
        {
                for(int j = 0; j < 3; j++)
                {
                        for(int k = 0 ; k < 3; k++)
                        {
                                int value = 1;
                                if( i == j || j == k || i == k)
                                        value = 0 ;
                                else
				{
                                        int nr_inversions = 0;
                                        nr_inversions += i>j ? 1:0;
                                        nr_inversions += i>k ? 1:0;
                                        nr_inversions += j>k ? 1:0;
                                        for(int ctr =0; ctr<nr_inversions; ctr ++ )
                                                value *= -1;
                                }

                                levichevi(i,j,k) = value ;
                        }

                }
	}

}
