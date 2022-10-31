#include<iostream>

#include"Algorithm.h"

using namespace std;

// muliplication number on matrix
Matrix& operator*(const double a, const Matrix& R)
{
    int rows = R.getRows();
    int cols = R.getCols();
    Matrix* res = new Matrix(R);

    for(int i=0;i<rows;i++)
        for(int j=0;j<cols;j++)
            (*res)(i,j)*=a;
    return *res;
}

int main(int argc, char const *argv[])
{
    string fileName = "RM_32_16";
    // Matrix A("Data/"+fileName+".gen");

    Matrix A(3,3,true,10);

    cout<<A<<endl;    

    Matrix B = AlgorithmViterbi::binMatrix(3);
    B.eraseRows(3,B.getRows());

    cout<<B<<endl;

    cout<<A.appendCols(B)<<endl;

    // A.toFile("Data/"+fileName+"_spen.gen");
    return 0;
}

