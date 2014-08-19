#include <cstdlib>
#include <iostream>
#include <fstream>
#include<iomanip>
#include <string>
#include "ltensor/LTensor.h"
#include <complex>


#define size 6
using namespace std;

int main()
{

    cout<<"+++++++++LTENSOR TEST+++++++++"<<endl;
    cout<<"Assignation routines test"<<endl;

    cout<<"----------Rank 1-----------"<<endl;

    Marray<double,1> a1(size),b1(size+1);

    const Marray<double,1> c1(size,7);

    Marray<int,1> d1(size+2,9);

	Marray<double,1> ds1(size+3,3.3);


	Marray<int,1> ds12(size+3,5);


    b1.sucesion(0);



	cout<<"///////////////////"<<endl;
    cout<<"Operation a1=b1"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;

    a1=b1;


    cout<<"a1:"<<a1<<endl;

    assert(a1==b1);



    cout<<"///////////////////"<<endl;
    cout<<"Operation a1=const c1"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"c1:"<<c1;
    cout<<"Performing op"<<endl;
    a1=c1;
    cout<<"a1:"<<a1<<endl;

    assert(a1==c1);

    cout<<"///////////////////"<<endl;
    cout<<"Operation a1=10"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"Performing op"<<endl;
    a1=10;
    cout<<"a1:"<<a1<<endl;

    for(int i=0;i<a1.get_dim1();i++)
        assert(a1(i)==10);


    cout<<"///////////////////"<<endl;
    cout<<"Operation a1(double)=d1(int)"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"d1:"<<d1;
    cout<<"Performing op"<<endl;
    a1=d1;
    cout<<"a1:"<<a1<<endl;

    for(int i=0;i<a1.get_dim1();i++)
        assert(a1(i)==d1(i));

	cout<<"///////////////////"<<endl;
    cout<<"Operation d1(int)=b1(double)"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"d1:"<<d1;
    cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;
    d1=b1;
    cout<<"d1:"<<d1<<endl;

    for(int i=0;i<d1.get_dim1();i++)
        assert(d1(i)==b1(i));


	cout<<"///////////////////"<<endl;
    cout<<"Operation d1(int)=b1(double)"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"d1:"<<d1;
    cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;
    d1=b1;
    cout<<"d1:"<<d1<<endl;

    for(int i=0;i<d1.get_dim1();i++)
        assert(d1(i)==b1(i));


	cout<<"///////////////////"<<endl;
    cout<<"Operation a1*=2"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"Performing op"<<endl;
    a1*=2;
    cout<<"a1:"<<a1<<endl;


	cout<<"///////////////////"<<endl;
    cout<<"Operation a1+=2"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"Performing op"<<endl;
    a1+=2;
    cout<<"a1:"<<a1<<endl;

	cout<<"///////////////////"<<endl;
    cout<<"Operation a1-=2"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"Performing op"<<endl;
    a1-=2;
    cout<<"a1:"<<a1<<endl;


	cout<<"///////////////////"<<endl;
    cout<<"Operation a1/=2"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
    cout<<"Performing op"<<endl;
    a1/=2;
    cout<<"a1:"<<a1<<endl;


	cout<<"///////////////////"<<endl;
    cout<<"Operation a1(i)=b1(i)"<<endl;
	cout<<"///////////////////"<<endl;
	Index <'i'> i;
	b1.resize(a1.get_dim1());
	b1.sucesion(0,1);
	cout<<"a1:"<<a1;
	cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;


    a1(i)=b1(i);
    cout<<"a1:"<<a1<<endl;

    assert(a1==b1);

	cout<<"///////////////////"<<endl;
    cout<<"Operation a1(i)=b1(i)*2"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
	cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;
    a1(i)=b1(i)*2;
    cout<<"a1:"<<a1<<endl;

    for(int i=0;i<a1.get_dim1();i++)
        assert(a1(i)==b1(i)*2);

	cout<<"///////////////////"<<endl;
    cout<<"Operation a1(i)=b1(i)/2"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
	cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;
    a1(i)=b1(i)/2;
    cout<<"a1:"<<a1<<endl;

    for(int i=0;i<a1.get_dim1();i++)
        assert(a1(i)==b1(i)/2);

    cout<<"///////////////////"<<endl;
    cout<<"Operation a1+=b1"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
	cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;
    a1+=b1;
    cout<<"a1:"<<a1<<endl;


    cout<<"///////////////////"<<endl;
    cout<<"Operation a1-=b1"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
	cout<<"b1:"<<b1;
    cout<<"Performing op"<<endl;
    a1-=b1;
    cout<<"a1:"<<a1<<endl;



	cout<<"///////////////////"<<endl;
    cout<<"Operation a1(double)=ds1(double) (different sizes)"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
	cout<<"ds1:"<<ds1;
    cout<<"Performing op"<<endl;
    a1=ds1;
    cout<<"a1:"<<a1<<endl;

    assert(a1==ds1);


	cout<<"///////////////////"<<endl;
    cout<<"Operation a1(double)=ds12(int) (different sizes/types)"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"a1:"<<a1;
	cout<<"ds12:"<<ds12;
    cout<<"Performing op"<<endl;
    a1=ds12;
    cout<<"a1:"<<a1<<endl;

    for(int i=0;i<a1.get_dim1();i++)
        assert(a1(i)==ds12(i));

	cout<<"///////////////////"<<endl;
    cout<<"Operation ds12(int)=b1(double) (different sizes/types)"<<endl;
	cout<<"///////////////////"<<endl;
    cout<<"b1:"<<b1;
	cout<<"ds12:"<<ds12;
    cout<<"Performing op"<<endl;
    ds12=b1;
    cout<<"ds12:"<<ds12<<endl;

	for(int i=0;i<ds12.get_dim1();i++)
        assert(b1(i)==ds12(i));


	cout<<endl;

//cin.get();


}
