#include"Matrix.h"
#include"ViterbiVerticle.h"

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

    static Matrix& getGrid(Matrix& spenG)
    {
        Matrix mass = getRowsStartEnd(spenG);
        int maxIndex = spenG.getCols();

        // define count verticles
        // cout<<mass;
        // int countVerticles = 1;
        // int cVlayer = 1;
        // for(int i=0;i<maxIndex;i++)
        // {
        //     bool mul2 = false;
        //     for(int r=0;r<mass.getRows() && !mul2;r++)
        //         mul2 |= mass(r,0)==i;
        //     if(mul2)
        //         cVlayer*=2;
        //     bool div2 = false;
        //     for(int r=0;r<mass.getRows() && !div2;r++)
        //         div2 |= mass(r,1)==i;
        //     if(div2)
        //         cVlayer/=2;
        //     countVerticles+=cVlayer;      
        // }

        ViterbiGrid grid(maxIndex + 1);

        int Nact = 0;
        grid.layers[0] = new ViterbiLayer();
        grid.layers[0].Nact = 0;
        
        grid.layers[0].verticles = new ViterbiVerticle();

        
        int r=0;
        for(int i = 1;i<grid.numLayers; i++)
        {

            grid.layers[i] = new ViterbiLayer();


            bool startFlag = false;
            int startXnum=0;
            for(;startXnum<mass.getRows() && !startFlag;startXnum++)
                startFlag |= mass(startXnum,0)==i;
            
            

            bool endFlag = false;
            int endXnum=0;
            for(;endXnum<mass.getRows() && !endFlag;r++)
                endFlag |= mass(endXnum,1)==i;
            
            
            grid.layers[i].Nact = grid.layers[i-1].Nact + int(startFlag) - int(endFlag);

            int Nact = grid.layers[i].Nact;
            if(Nact!=0)
            {
                grid.layers[i].activeX = new int[Nact];

                for(int j=0;j<grid.layers[i-1].Nact;j++)
                    if(grid.layers[i].activeX[j] != endXnum)
                        grid.layers[i].activeX[j] = grid.layers[i].activeX[j];
                
                if(startFlag)
                    grid.layers[i].activeX[Nact-1] = startXnum;


                grid.layers[i].verticles = new ViterbiVerticle[pow(2,Nact)];
                if(!startFlag && !endFlag)
                {
                    for(int j=0;j<Nact;j++)
                        grid.layers[i-1].verticles[j].connectNext(  grid.layers[i].verticles[j]  );
                    
                }
                if(startFlag && !endFlag)
                {
                    for(int j=0;j<Nact;j++)
                        grid.layers[i-1].verticles[j/2].connectNext(   grid.layers[i].verticles[j]  );
                }
                if(!startFlag && endFlag)
                {
                    // do this
                }
                if(startFlag && endFlag)
                {
                    // do this
                }
            }
            

        }

        return spenG;
    }

};

