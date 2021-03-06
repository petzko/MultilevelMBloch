#include <iostream>
#define CHECK_OOB
#include<fstream>

using namespace std;

//macros definitions

//order order 3 vs 2

int main() {

    ifstream in;
    ofstream out;

    char a,b,c,d,e;

// type A3(i,j,k)*B(k)
    in.open("3vs2.input");
    out.open("3vs2.cpp");
out<<" #include \"TExprs.h\" \n \
#include <iostream> \n \
#define CHECK_OOB \n \
#include<fstream> \n \
int main(void){   \n \
Marray<double,3> A3(3,3,3);\n \
Marray<double,2> B2(3,3); \n \
Marray<double,1> C1Test(3);\n \
Marray<double,1> C1(3); \n \
A3.sucesion(0,1); \n \
B2.sucesion(0,1); \n \
IndexG<'i'> iG; \n \
IndexG<'j'> jG; \n \
IndexG<'k'> kG; \n \
";




    if(!in.is_open()){
        cout<<"Error";
        return 0;
        }

    in>>a;  in>>b;  in>>c;  in>>d;  in>>d;  in>>e;

    while(!in.eof()){
        out<<"for (int i=0;i<3;i++) { \n \
    for (int j=0;j<3;j++){ \n \
        for(int k=0;k<3;k++){ \n ";
        out<<"              C1(i)+=A3("<<a<<","<<b<<","<<c<<")*B2("<<d<<","<<e<<"); \n";
        out<<"          }\n     }\n } \n\n";

        out<<"C1Test(iG)=A3("<<a<<"G,"<<b<<"G,"<<c<<"G)*B2("<<d<<"G,"<<e<<"G); \n";
        out<<"assert(C1==C1Test); \n\n\n";

        out<<"C1=0;\n\n\n";
        in>>a;  in>>b;  in>>c;  in>>d;  in>>d;  in>>e;

        }

    out<<"}\n";


    in.close();
    out.close();

    }
