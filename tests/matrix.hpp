#ifndef _MATRIX_
#define _MATRIX_

#include <complex>
#include <iostream>
#include <stdexcept>

#define MBLOCH_OUT_WRN(str,FILE,LINE){\
	std::cout<<"\nWARNING: (in FILE:" << FILE << "; LINE:" << LINE << ") :\n" << str << "\n";}

#define MBLOCH_OUT_ERR(str,FILE,LINE){\
	std::cout<<"FILE:" << FILE << " LINE:" << LINE << ": \n ERROR:" << str << "\n";}


namespace mbloch{


	template<unsigned int M, unsigned int N,typename _Tp>
		class Matrix{

			//define the data container... -> column wise ordered 
			// M*N matrix of complex numbers  
			_Tp data[M*N];
			public:

			// initialize an empty complex Nr matrix
			Matrix()
			{
				for(int i = 0; i<M; i++)
					for(int j = 0; j< N; j++)	
						data[i*N + j] = 0.0;
			}

			// overload some operators for convenience 

			// make the asignment operator copy the data... 
			// coppies the data from lhs to rhs...	
			Matrix<M,N,_Tp>& operator= (Matrix<M,N,_Tp> arg)
			{
				if(this != &arg)
					for(int i = 0; i<M; i++)
						for(int j = 0; j< N; j++)
							data[i*N + j] = arg(i,j);	
			};  

			// make the function call operator retrive the i,j th data element
			_Tp& operator()(int i,int j)
			{ 

				if(i >= M || i < 0 ) 
				{	
					MBLOCH_OUT_ERR("Array index out of bounds. Cannot retrieve Matrix element. ",__FILE__,__LINE__)
						throw std::out_of_range("row index out of bounds.");

				}				

				if(j >= N || j < 0 )
				{
					MBLOCH_OUT_ERR("Array index out of bounds. Cannot retrieve Matrix element. ",__FILE__,__LINE__)
						throw std::out_of_range("column index out of bounds.");
				}

				return data[i*N+j];

			}

			// overload the multiplication by matrix operator... 
			template<unsigned int L>
				Matrix<M,L,_Tp> operator*(Matrix<N,L,_Tp>& arg)
				{

					Matrix<M,L,_Tp> res;
					for (int i =0; i< M;i++)
						for(int j = 0; j< L; j++)
							for(int k = 0 ; k< N; k++)
								res(i,j) += data[i*N + k]*arg(k,j);		

					return res;
				}


			//overload the addition operator
			Matrix<M,N,_Tp> operator+(Matrix<M,N,_Tp> arg)
			{

				Matrix<M,N,_Tp> res;
				for (int i =0; i< M;i++)
					for(int j = 0; j< N; j++)
						res(i,j) = data[i*N+j] + arg(i,j);                 	
				return res;
			}

			//overload the subtraction operator
			Matrix<M,N,_Tp> operator-(Matrix<M,N,_Tp> arg)
			{
				Matrix<M,N,_Tp> res;
				for (int i =0; i< M;i++)
					for(int j = 0; j< N; j++)
						res(i,j) = data[i*N+j] - arg(i,j);
				return res;
			}


			//overload the print operator

			friend std::ostream& operator<<(std::ostream& os,const mbloch::Matrix<M,N,_Tp>& obj)
			{

				for(int i = 0; i< M;i++)
				{       
					for(int j = 0 ; j < N;j++)
						os <<obj.data[i*N+j] << " ";
					os<<"\n";
				}

				return os;
			}




		}; 

}


template<unsigned int M, typename _Tp> 
mbloch::Matrix<M,M,_Tp> eye(){
	mbloch::Matrix<M,M,_Tp> res; 
	 for(int i = 0;i<M;i++)
                res(i,i) = 1.0;
	return res; 
}


template<unsigned int M, unsigned int N,typename _Tp>
mbloch::Matrix<M,N,_Tp> operator*(const _Tp & lhs,mbloch::Matrix<M,N,_Tp>& rhs ){

	mbloch::Matrix<M,N,_Tp> res;
	for (int i =0; i< M;i++)
		for(int j = 0; j< N; j++)
			res(i,j) = lhs*rhs(i,j);                      

	return res;


}


template<unsigned int M,unsigned int N,typename _Tp>
mbloch::Matrix<M,N,_Tp> operator*(mbloch::Matrix<M,N,_Tp>& lhs, const _Tp& rhs){
	return rhs*lhs;

}

template<unsigned int M,unsigned int N,unsigned int K,unsigned int L,typename _Tp>
mbloch::Matrix<M*K,N*L,_Tp> tensor_product(mbloch::Matrix<M,N,_Tp> A,mbloch::Matrix<K,L,_Tp> B){

	mbloch::Matrix<M*K,N*L,_Tp> res;

	for(int i = 0 ; i< M; i++)
		for(int j = 0; j< N; j++)
		{	
			for(int k = 0; k< K; k++)
			{	for(int l = 0; l< L; l++)	
				{	
					res(i*K+k,j*L+l) = A(i,j)*B(k,l);
				}
			}

		}

	return res; 
}


template< unsigned int M, unsigned int K,typename _Tp>
mbloch::Matrix<M*K,M*K,_Tp> kronecker_sum(mbloch::Matrix<M,M,_Tp> A,mbloch::Matrix<K,K,_Tp> B){



	// initialize the identity matrices...
	mbloch::Matrix<M,M,_Tp> EyeA = eye<M,_Tp>();
	mbloch::Matrix<K,K,_Tp> EyeB = eye<K,_Tp>();


	return tensor_product(A,EyeB) + tensor_product(EyeA,B);
}

// here the dereferencing operator shall mean conjugate transpose... 
template< unsigned int M,typename _Tp>
mbloch::Matrix<M,M,_Tp> operator*(mbloch::Matrix<M,M,_Tp>& arg){

	mbloch::Matrix<M,M,_Tp> res;
	for(int i = 0; i< M; i++)
		for(int j = 0; j<M;j++)
			res(j,i) = std::conj(arg(i,j));
	return res;

}


#endif
