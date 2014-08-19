#include <cstdlib>
#include <iostream>
#include <fstream>
#include<iomanip>
#include <string>
#include "ltensor/LTensor.h"
#include <complex>


using namespace std;
#define size 6




int main(void){

	//Simple contraction
	Marray<double,3> A3(3,3,3);
	Marray<double,2> B2(3,3);
	Marray<double,1> C1Test(3);

	A3.sucesion(0,1);
	B2.sucesion(0,1);
	//declare the indexes to be used
	//the char parameter represents the contraction index
	Index<'i'> iG;
	Index<'j'> jG;
	Index<'k'> kG;
	C1Test(iG)=A3(iG,jG,kG)*B2(jG,kG);


	//Create the arrays
	Marray<double,4> A4(3,3,3,3);
	Marray<double,1> B1(3);
	Marray<double,3> C3Test(3,3,3);
 
	//fill 'em
	A4.sucesion(0,1);
	B1.sucesion(0,1);
	Index<'l'> lG;
	
	//perform the operation
	C3Test(iG,jG,kG)=A4(iG,jG,lG,kG)*B1(lG);


}
