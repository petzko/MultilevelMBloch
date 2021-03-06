#include <iostream>
#define CHECK_OOB
#include<fstream>

using namespace std;

//macros definitions

//order order 4 vs 3

int main() {

    ifstream in;
    ofstream out;

    char a,b,c,d,e,f,g;

// type A3(i,j,k)*B(k)
    in.open("4vs3.input");
    out.open("4vs3.cpp");
out<<" #include \"TExprs.h\" \n \
#include <iostream> \n \
#define CHECK_OOB \n \
#include<fstream> \n \
int main(void){   \n \
Marray<double,4> A4(3,3,3,3);\n \
Marray<double,3> B3(3,3,3); \n \
Marray<double,3> C3Test(3,3,3);\n \
Marray<double,3> C3(3,3,3); \n \
A4.sucesion(0,1); \n \
B3.sucesion(0,1); \n \
IndexG<'i'> iG; \n \
IndexG<'j'> jG; \n \
IndexG<'k'> kG; \n \
IndexG<'l'> lG; \n \
IndexG<'m'> mG; \n \
";




    if(!in.is_open()){
        cout<<"Error";
        return 0;
        }

    in>>a;  in>>b;  in>>c;  in>>d;  in>>e;  in>>e; in>>f; in>>g;

    while(!in.eof()){
        out<<"for (int i=0;i<3;i++) { \n \
    for (int j=0;j<3;j++){ \n \
        for (int l=0;l<3;l++){ \n\
            for (int m=0;m<3;m++){ \n\
                for(int k=0;k<3;k++){ \n ";

        out<<"                C3(i,j,m)+=A4("<<a<<","<<b<<","<<c<<","<<d<<")*B3("<<e<<","<<f<<","<<g<<"); \n";
        out<<"           }\n        }\n     }\n        }\n} \n\n";

        out<<"C3Test(iG,jG,mG)=A4("<<a<<"G,"<<b<<"G,"<<c<<"G,"<<d<<"G)*B3("<<e<<"G,"<<f<<"G,"<<g<<"G); \n";
        out<<"assert(C3==C3Test); \n\n\n";

        out<<"C3=0;\n\n\n";
        in>>a;  in>>b;  in>>c;  in>>d;  in>>e;  in>>e; in>>f;   in>>g;

        }

    out<<"}\n";


    in.close();
    out.close();

    }
