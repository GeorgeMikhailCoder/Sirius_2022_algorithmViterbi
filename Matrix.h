#include<iostream>
#include<fstream>
#include <string> 
#include<iomanip>
#include <cmath>
using namespace std;

const bool defaultRandom = true;
const int defaultRandomSqueeze = 2;

class Matrix
{
protected:
    int rows;
    int cols;
    double** mass;
public:

// basic class functions
    Matrix(int r, int c, bool isRand=defaultRandom, int randomSqueeze = defaultRandomSqueeze)
    {
        srand(time(NULL));
        rows = r;
        cols = c;
        mass = new double*[rows];
        for(int i=0;i<rows;i++)
        {
            mass[i] = new double[cols];
            for(int j=0;j<cols;j++)
                if(isRand)
                    mass[i][j]=rand()%randomSqueeze;
                else
                    mass[i][j]=0;
        }
    }

    Matrix(const Matrix&R)
    {
        rows = R.rows;
        cols = R.cols;
        mass = new double*[rows];
        for(int i=0;i<rows;i++)
        {
            mass[i] = new double[cols];
            for(int j=0;j<cols;j++)
                mass[i][j]=R(i,j);
        }
    }

    Matrix(string fileName)
    {
        ifstream f(fileName);
        if(!f.is_open())
            throw("Matrix: Wrong input file");
        else
        {
            int r,c;
            f>>c>>r;

            rows = r;
            cols = c;
            mass = new double*[rows];
            for(int i=0;i<rows;i++)
            {
                mass[i] = new double[cols];
                for(int j=0;j<cols;j++)
                {
                    int smb = f.get()-48;
                    if(smb==0 || smb==1)
                        mass[i][j]=smb;
                    else
                        j--;
                }
            }
        }
        f.close();

    }

    ~Matrix()
    {
         for(int i=0;i<rows;i++)
            delete []mass[i];
        delete []mass;
    }
    
    string getSizeS()
    {
        return "("+ to_string(rows) + ", " +  to_string(cols) + ")";
    }

    int getRows() const 
    {
        return rows;
    }

    int getCols() const
    {
        return cols;
    }

    Matrix& slice(int rs, int rf, int cs, int cf)
    {// returns a copy of a part of existing matrix
        
        if(rs>=rows || rf > rows || cs >= cols || cf > cols || rs > rf || cs > cf) // = для последних?
            throw("Matrix slice: incorrect boundaries of slice");
        
        Matrix* res = new Matrix(rf-rs, cf-cs);
        for(int i=0;i<res->rows;i++)
            delete [] res->mass[i];
        delete [] res->mass;

        res->mass=this->mass;

        res->mass += rs;
        res->mass[0] += cs;
        
        return *res;
    }

    Matrix& getStb(int c)
    {
        return this->slice(0,rows,c,c+1);
    }

    Matrix& copy()
    {
        Matrix* res = new Matrix(rows,cols);

        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                (*res)(i,j)=this->mass[i][j];
        return *res;
    }


    Matrix& erase(int rs=-1, int rf=-1, int cs=-1, int cf=-1)
    {
        if(rs>=rows || rf > rows || cs >= cols || cf > cols || rs > rf || cs > cf) // = для последних?
            throw("Matrix erase: incorrect boundaries of erase");
        double **prevMass = mass;
        int prevRows = rows;
        int prevCols = cols;
        
        auto beforeColSpace = [cs=cs, cf=cf](int j)->int {
             return int((cs <= j) && (cs!=-1) && (cf!=-1));
             };
        auto beforeRowSpace = [rs=rs,rf=rf](int i)->int { 
            return int((rs <= i)  && (rs!=-1) && (rf!=-1));
            };

        rows = rows - (rf - rs) * int((rs!=-1) && (rf!=-1));
        cols = cols - (cf - cs) * int((cs!=-1) && (cf!=-1));


        int a(0),b(0);
        mass = new double*[rows];
        for(int i=0;i<rows;i++)
        {
            mass[i] = new double[cols];
            for(int j=0;j<cols;j++)
                {
                    a = i + (rf - rs)*beforeRowSpace(i);
                    b = j + (cf - cs)*beforeColSpace(j);
                    mass[i][j] = prevMass[a][b];
                }
        }
        for(int i=0;i<prevRows;i++)
            delete []prevMass[i];
        delete []prevMass;
        return *this;
    }

    Matrix& eraseRows(int r, int rf = -1)
    {
        if(rf==-1)
            rf = r + 1;
        return erase(r,rf,-1,-1);
    }
    
    Matrix& eraseCols(int c, int cf = -1)
    {
        if(cf == -1)
            cf = c + 1;
        return erase(-1,-1,c,cf);
    }


    Matrix& appendRows(Matrix& rowsMatrix, int r = -1)
    {
        if(r>rows || r < -1)
            throw("Matrix appendRows: incorrect row number");
        if(cols != rowsMatrix.getCols())
            throw("Matrix appendRows: cols number must be equal");

        if(r==-1) 
            r = rows;

        double **prevMass = mass;
        int prevRows = rows;
        
        rows += rowsMatrix.getRows();        

        int i_prev(0), i_new(0);
        mass = new double*[rows];
        for(int i=0;i<rows;i++)
        {
            if(i<r || i>=r+rowsMatrix.getRows())
            {
                mass[i] = prevMass[i_prev];
                i_prev+=1;
            }
            else
            {
                mass[i] = new double[cols];
                for(int j=0;j<cols;j++)
                    mass[i][j] = rowsMatrix(i_new,j);
                i_new+=1;
            }
        }
        return *this;
    }

    Matrix& appendCols(Matrix& colsMatrix, int c = -1)
    {
        if(c>cols || c < -1)
            throw("Matrix appendCols: incorrect col number");
        if(rows != colsMatrix.getRows())
            throw("Matrix appendCols: rows number must be equal");

        if(c==-1) 
            c = cols;


        double **prevMass = mass;
        int prevRows = rows;
        int prevCols = cols;


        cols += colsMatrix.getCols();

        int j_prev(0), j_new(0);
        mass = new double*[rows];
        for(int i=0;i<rows;i++)
            {
                mass[i] = new double[cols];
                for(int j=0;j<cols;j++)
                    if(j<c || j>=c+colsMatrix.getCols())
                    {
                        mass[i][j] = prevMass[i][j_prev];
                        j_prev++;
                    }
                    else
                    {
                        mass[i][j] = colsMatrix(i,j_new);
                        j_new++;
                    }
            }

        for(int i=0;i<rows;i++)
            delete []prevMass[i];
        delete []prevMass;
        return *this;
    }


    void toFile(string fileName)
    {
        ofstream f(fileName);
        if(!f.is_open())
            throw("Matrix: Wrong ounput file");
        else
        {
            f<<cols<<" "<<rows<<endl;
            
            for(int i=0;i<rows;i++)
            {
                for(int j=0;j<cols;j++)
                    f<<mass[i][j];
                f<<endl;
            }
            f.close();
        }
    }

// basic class operators
    double operator ()(int r, int c) const
    {
        return this->mass[r][c];
    }

    double& operator ()(int r, int c)
    {
        return this->mass[r][c];
    }

    Matrix& operator ()(int r)
    {
        return this->slice(r,r+1,0,cols);
    }

    Matrix& operator ()(int rs, int rf, int cs, int cf)
    {
        return this->slice(rs,rf,cs,cf);
    }

    friend ostream& operator <<(ostream& os, const Matrix& R)
    {
        os<<"Matrix with size ("<<R.rows<<", "<<R.cols<<")";
        if(R.rows*R.cols<100)
        {
            os<<endl;
            for(int i=0;i<R.rows;i++)
            {
                for(int j=0;j<R.cols;j++)
                    os<<setw(7)<<setprecision(3)<<R(i,j)<<"\t";
                os<<endl;
            }
        }
        else
            os<<", big matrix"<<endl;
        return os;
    }

    Matrix& operator =(const Matrix& R)
    {
        if(R.rows==rows && R.cols==cols)
        {
            for(int i=0;i<rows;i++)
            {
                for(int j=0;j<cols;j++)
                    mass[i][j]=R(i,j);
            }
        }
        else
        {
            rows = R.rows;
            cols = R.cols;
            mass = new double*[rows];
            for(int i=0;i<rows;i++)
            {
                mass[i] = new double[cols];
                for(int j=0;j<cols;j++)
                    mass[i][j]=R(i,j);
            }
        }

        return *this;
    }

// static matrix functions
    static Matrix& eye(int n)
    {
        Matrix* res = new Matrix(n,n);
        for(int i=0;i<n;i++)
            (*res)(i,i)=1;
        return *res;
    }


// arifmetic operators
    Matrix& operator + (const Matrix& R)
    {
        if(!(cols == R.cols | rows==R.rows))
        {
            cout<<"Wrong addition"<<endl;
            throw("Wrong matrix addition");
            return *this;
        }
        Matrix* res = new Matrix(rows,cols);

        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                (*res)(i,j)=this->mass[i][j] + R(i,j);
        return *res;
    }

    Matrix& operator +=(const Matrix& R)
    {
        if(!(cols == R.cols | rows==R.rows))
            throw("Wrong matrix addition");

        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                this->mass[i][j]=this->mass[i][j] + R(i,j);
        return *this;
    }

    Matrix&  operator * (const Matrix& R)
    {
        if(cols != R.rows)
        {
            cout<<"Wrong multiplication"<<endl;
            throw("Wrong matrix multiplication");
            return *this;
        }
        
        int N = cols;
        int newRows = rows;
        int newCols = R.cols;
        Matrix* res = new Matrix(newRows,newCols);

        for(int i=0;i<newRows;i++)
            for(int j=0;j<newCols;j++)
            {
                double sum = 0;
                for(int k=0;k<N;k++)
                    sum+=this->mass[i][k]*R(k,j);
                (*res)(i,j)=sum;
            }
        return *res;
    }

    Matrix&  operator * (const double a)
    {
        Matrix* res = new Matrix(*this);

        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                (*res)(i,j)*=a;
        return *res;
    }

    Matrix& operator % (const int number)
    {
        Matrix* res = new Matrix(*this);

        for(int i=0;i<rows;i++)
            for(int j=0;j<cols;j++)
                (*res)(i,j)= int((*res)(i,j))%number;
        return *res;
    }
    

};

class VectorS: public Matrix
{
public:
    VectorS(int N,bool isRand=defaultRandom,int randomSqueeze=defaultRandomSqueeze):Matrix(1,N, isRand, randomSqueeze){};
    VectorS(const Matrix& R):Matrix(R){};

    double& operator()(int k)
    {
        return Matrix::operator()(0,k);
    }
};
