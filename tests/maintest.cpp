#include "matrix.hpp"

int main(){

	mbloch::Matrix<2,2,double> H; 
	mbloch::Matrix<2,2,double> V; 
	mbloch::Matrix<2,4,double> B;
	mbloch::Matrix<2,1,double> vect;
	mbloch::Matrix<2,2,double> Eye; 
	Eye(0,0) = 1.0; 
	Eye(1,1) = 1.0; 

	std::complex<double> I(0.0,1.0);
	H(0,1) = 1.0; 
	H(1,0) = 1.0; 
	H(0,0) = 1.0; 
	H(1,1) = 2.0; 
	std::cout<< "\nPrint H: \n" << H;	
	
	V(0,0) = -1.0;// + I;
	V(1,1) = 1.0;//- I;  
	std::cout<<"\nPrint V:\n" << V;	
	
	mbloch::Matrix<2,2,double> tmp = H;
	std::cout<<"\nPrint H copy:\n" << tmp;
	tmp = H*V;
	std::cout<<"\nPrint H*V:\n" << tmp;
	std::cout<<"\nPrint H*B:\n" << H*B;
	std::cout<<"\nPrint H*vect:\n" << (H*vect);

	std::cout<<"\nPrint H+V:\n" << (H+V);
	std::cout<<"\nPrint H-V:\n" << (H-V);	

//	std::complex<double> nr = 2.0; 
	std::cout<<"\nPrint H*2:\n"<<H*2.0;
	std::cout<<"\nPrint 2*H:\n"<<2.0*H;

	mbloch::Matrix<4,4,double> kronProd= tensor_product(Eye,H);  
	std::cout<<"\nPrint tensor product: \n"<< kronProd;

	std::cout<<"\nPrint Kronecker Sum: \n"<<kronecker_sum(H,V);


	mbloch::Matrix< 2 , 2, std::complex<double> > complexH ; 
	complexH(0,1) = -1.0+I;
	complexH(1,0) = 1.0-I;
	std::cout<<"\nPrint complex Matrix cH:\n"<<complexH;
	std::cout<<"\bPrint complex conjugate matrix cH*:\n"<<*complexH;
	

	return 0; 
}
