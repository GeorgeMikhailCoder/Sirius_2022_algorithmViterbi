class ViterbiVerticle
{
public:
    ViterbiVerticle* next[2]; // указатели на следующие вершины в решётке

    ViterbiVerticle()
    {
        next[0] = nullptr;
        next[1] = nullptr;
    };

    void connectNext(ViterbiVerticle& v)
    {
        if(next[0]==nullptr)
            next[0]=&v;
        else if(next[1]==nullptr)
            next[1]=&v;
        else
            throw("ViterbiVerticle connectNext: next array already full");
    }

};

class ViterbiLayer
{
public:
    int* activeX;             // массив индексов активных х
    int Nact;                 // количество активных х
    ViterbiVerticle *verticles; 

    ViterbiLayer()
    {
        verticles = nullptr;
        activeX = nullptr;
        Nact=0;
    }
    ~ViterbiLayer()
    {
        delete[] verticles;
    }

};

class ViterbiGrid
{
public:
    ViterbiLayer* layers;
    int numLayers;

    ViterbiGrid(int N)
    {
        layers = new ViterbiLayer[N];
        numLayers=N;
    }
    ~ViterbiGrid()
    {
        delete []layers;
    }

    void toFile(string fileName)
    {
        ofstream f(fileName);
        if(!f.is_open())
            throw("Matrix: Wrong ounput file");
        else
        {
            f<<"grid start"<<endl;
            f<<"layers count "<<numLayers<<endl;
            
            for(int i=0; i<numLayers; i++)
            {
                f<<endl<<"layer "<<i<<": ";
                f<<"Nact="<<layers[i].Nact<<", active X = [";
                for(int j=0;j<layers[i].Nact-1;j++)
                    f<<layers[i].activeX[j]<<", ";
                if(layers[i].Nact != 0)
                    f<<layers[i].activeX[layers[i].Nact-1]<<"]"<<endl;
                
                for(int j=0;j<pow(2,layers[i].Nact);j++)
                    f<<"Verticle "<<j<<": my adress = "<<layers[i].verticles+j<<", next_1 = "<<layers[i].verticles[j].next[0]<<", next_2 = "<<layers[i].verticles[j].next[1]<<endl;

            }
            
            f.close();
        }
    }
};
