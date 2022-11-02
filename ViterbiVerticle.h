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

};
