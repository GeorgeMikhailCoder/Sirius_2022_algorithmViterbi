#include"Matrix.h"


int pow(const int a, const unsigned int n)
{
    int res = 1;
    for(int i=0;i<n;i++)
        res*=a;
    return res;
}

class AlgorithmViterbi
{

public:
    static::Matrix& getRowsStartEnd(const Matrix& R)
    {
        Matrix mass(R.getRows(),2);

        for(int r=0;r<R.getRows();r++)
        {
            int c_start = 0;
            for(;c_start<R.getCols() && R(r,c_start)==0;c_start++);
            if(c_start==R.getCols())
            {
                mass(r,0)=-1;
                mass(r,1)=-1;  
                continue;
            }
            
            int c_end = R.getCols()-1;
            for(;c_end>=0 && R(r,c_end)==0; c_end--);

            mass(r,0)=c_start;
            mass(r,1)=c_end;            
        }
        return *(new Matrix(mass));
    }

    static bool checkSpen(const Matrix& R)
    {
        Matrix mass(R.getRows(),2);

        for(int r=0;r<R.getRows();r++)
        {
            int c_start = 0;
            for(;c_start<R.getCols() && R(r,c_start)==0;c_start++);
            if(c_start==R.getCols())
                return false;
            
            int c_end = R.getCols()-1;
            for(;c_end>=0 && R(r,c_end)==0; c_end--);
            if(c_end==-1)
                return false;

            mass(r,0)=c_start;
            mass(r,1)=c_end;            
        }

        for(int r = 0; r<mass.getRows(); r++)
            for(int i=r+1;i<mass.getRows(); i++)
                if(mass(r,0) == mass(i,0) || 
                   mass(r,0) == mass(i,1) || 
                   mass(r,1) == mass(i,0) || 
                   mass(r,1) == mass(i,1))
                   return false;
        return true;
    };

    static void spen(Matrix& R)
    {
        if(R.getRows()>R.getCols() && 2*R.getRows()<R.getCols())
            throw("Algorithm: cols must be greater than rows and lower than 2*rows");

        VectorS str0 = R(0);

        // forward gauss
        int k=0;
        for(int c=0; c<R.getRows()+k; c++) // getRows - by the columns while col==row
        {
            int r=c-k;
            int r0 = r;
            for(;r<R.getRows() && R(r,c)==0;r++);
            // r - index of non zero str
            if(r==R.getRows())
            {
                k+=1;
                continue;
            }

            VectorS str0 = R(r0);
            R(r0) = R(r);
            R(r) = str0;
            r=r0;
            for(int i=0;i<R.getRows();i++)
                if(R(i,c)!=0 && i!=r)
                    R(i) = (R(i) + R(r))%2;           
        }

        //cout<<"Before back"<<endl<<R;
        R.toFile("tmp.gen");

        // backward gauss
        VectorS UsedRows = VectorS(R.getRows(),false); // 0 - unused, 1 - used

        for(int c=R.getCols()-1; c>=R.getRows()-k; c--) // cols .. rows
        {
            int r = R.getRows() - 1; // rows .. 1
            int r0 = r;
            
            R.toFile("tmp.gen");
            double a =  UsedRows(r);
            bool f = UsedRows(r) == 1;
            for(;r>=0 && (R(r,c)==0 || UsedRows(r) == 1);r--);        
                    
            // r - first non zero str
            if(r==-1)
                continue;

            UsedRows(r) = 1;
            for(int i=r-1; i>=0; i--)
                if(R(i,c) == 1)
                    R(i) = (R(i) + R(r))%2;
        }

        
        // sorting
        for(int c=0;c<R.getRows();c++)
        {
            int r=c;
            for(;r<R.getRows() && R(r,c)==0;r++);
            if(r==R.getRows())
                continue;
            VectorS str0 = R(c);
            R(c) = R(r);
            R(r) = str0;
        }

    }


    static Matrix& binMatrix(int n)
    {
        Matrix res(pow(2,n), n);
        int half = pow(2,n-1);
        bool smb = false;
        for(int c=0;c<res.getCols();c++)
        {
            int r=0;
            for(;r<res.getRows();)
            {
                for(int i=0;i<half;i++,r++)
                {
                    res(r,c) = int(smb);
                }
                smb=!smb;

                for(int i=0;i<half;i++,r++)
                {
                    res(r,c) = int(smb);
                }
                smb=!smb;
            }
            half /= 2;
        }

        Matrix* resPtr = new Matrix(res);
        return *resPtr;
    }

    static void smezh_i(const Matrix& smezh_prev, const Matrix& xactmass, bool kst,int xst, bool kend,int xend)
    {
        
    }

};