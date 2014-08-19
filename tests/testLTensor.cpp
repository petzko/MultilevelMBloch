#include <iostream>
#include "./LTensor/LTensor.h"



int main(){

	// create the rank 3 tensor
	Marray<double,3 > eps(3,3,3); 

	for(int i = 0 ; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			for(int k = 0 ; k < 3; k++)
			{

				int value = 1; 
				if( i == j || j == k || i == k)
					value = 0 ; 
				else{
					int nr_inversions = 0; 
					nr_inversions += i>j ? 1:0;
					nr_inversions += i>k ? 1:0;
					nr_inversions += j>k ? 1:0;		
					for(int ctr =0; ctr<nr_inversions; ctr ++ )
						value *= -1;
				}


				eps(i,j,k) = value ;

			}

		}


	}

	std::cout<< eps ; 

	Marray<double,1> A(3); 
	A(0) = 1;
	Marray<double,1> B(3); 
	B(1) = 1; 
	Marray<double,1> C(3); 
	Index<'i'> i; Index<'j'> j; Index<'k'> k; 
	//vector product of A and B 
	C(i) = eps(i,j,k)*A(j)*B(k);
	std::cout << "e1 x e2 = e3 vector product...\n ";
	std::cout << C; 
	
	C(i) = A(i) + 2*B(i) ;
	std::cout<<"e1 + 2* e2 = e3 addmultiply operation test \n"; 
	std::cout<<C; 
	
	
	
	return 0; 

}
