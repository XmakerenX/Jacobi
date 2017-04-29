#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include<sys/time.h>
#include <pthread.h>

#define N 500

struct matrix
{
    double* arr;

    matrix()
    {
        // allocate the matrix as one array
        arr = new double[N*N];
    }

    ~matrix()
    {
        delete[] arr;
    }

    // make accessing arr like accessing to a matrix
    double* operator[](int i)
    {
        return arr + N*i;
    }

};

class Thread
{
    pthread_t thread_;
};

//-----------------------------------------------------------------------------
// Name : getTime ()
//-----------------------------------------------------------------------------
double getTime ()
{
    struct timeval x;
    gettimeofday(&x, NULL);
    return (x.tv_sec * 1000000u + x.tv_usec ) / 1.e6;
}

//-----------------------------------------------------------------------------
// Name : maxind ()
// Desc : index of largest off-diagonal elemnt in row k
//-----------------------------------------------------------------------------
int maxind(matrix& S,int k)
{
    int m = k +1;
    for (int i = k + 2; i < N; i++)
    {
        if (std::abs(S[k][i]) > std::abs(S[k][m]) )
        {
            m = i;
        }
    }

    return m;
}

//-----------------------------------------------------------------------------
// Name : update ()
// Desc : update e[k] to its status
//-----------------------------------------------------------------------------
void update(int k, double t, int& state, double e[N], bool changed[N])
{
    double y = e[k];
    e[k] = y + t;

    if (changed[k] && y == e[k])
    {
        changed[k] = false;
        state--;
    }
    else if (!changed[k] && y != e[k])
    {
        changed[k] = true;
        state++;
    }
}

//-----------------------------------------------------------------------------
// Name : rotate ()
// Desc : perform rotation of s[i][j] s[k][l]
//-----------------------------------------------------------------------------
inline void rotate(matrix &S,int k, int l, int i, int j, double& c, double& s)
{
    // | S(kl) |   | c -s | |S(kl) |
    // |       | = |      | |      |
    // | S(ij) |   | s  c | |S(ij) |
    double skl,sij;

    skl = c * S[k][l] -s * S[i][j];
    sij = s * S[k][l] +c * S[i][j];

    S[k][l] = skl;
    S[i][j] = sij;
}

//-----------------------------------------------------------------------------
// Name : jacobi ()
//-----------------------------------------------------------------------------
void jacobi(matrix& S, double e[N], matrix& E)
{
    int i,k,l,m;
    double s,c,t,p,y,d,r;

    int state = N;

    int* ind = new int[N];
    bool* changed = new bool[N];

    int count = 0;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            E[i][j] = 0;

    // E = I
    for (int i = 0; i < N; i++)
        E[i][i] = 1;

    //state = N;

    for (int k = 0; k < N; k++)
    {
        ind[k] = maxind(S,k);
        e[k] = S[k][k];
        changed[k] = true;
    }

    while (state !=  0)
    {
        m = 0;
        for (k = 1; k < N-1; k++)
        {
            if (std::abs(S[k][ind[k]]) > std::abs(S[m][ind[m]]) )
            {
                m = k;
            }
        }

        k = m;
        l = ind[m];
        p = S[k][l];

        // calculate c = cos o, s = sin o
        double pp = p*p;
        y =  ( e[l] - e[k] ) / 2;
        d = std::abs(y) + std::sqrt(pp + y*y);
        r = std::sqrt(pp + d*d);
        c = d / r;
        s = p / r;
        t = pp / d;

        if ( y < 0)
        {
            s = -s;
            t = -t;
        }

        S[k][l] = 0.0f;
        update(k, -t, state, e, changed);
        update(l, t, state, e, changed);
        // rotate rows and columns k and l
        for (int i = 0; i < k; i++)
        {
            rotate(S,i,k,i,l, c, s);
        }

        for (int i = k +1; i < l; i++)
        {
            rotate(S,k,i,i,l, c, s);
        }

        for (int i = l + 1; i < N; i++)
        {
            rotate(S,k,i,l,i, c, s);
        }

        // rotate eigenvectors
        for (int i = 0; i  < N; i++)
        {
            double eik, eil;
            eik = c * E[i][k] - s * E[i][l];
            eil = s * E[i][k] + c * E[i][l];

            E[i][k] = eik;
            E[i][l] = eil;
        }

        // rows k, l have changed, update rows indK, indL
        ind[k] =  maxind(S,k);
        ind[l] =  maxind(S,l);

        count++;
    }

    delete[] ind;
    delete[] changed;

   std::cout << "count = "<< count << "\n\n";
}

//-----------------------------------------------------------------------------
// Name : readMatrixFromFile ()
//-----------------------------------------------------------------------------
bool readMatrixFromFile(matrix& S, std::string filePath)
{
    std::ifstream file(filePath);

    if (!file.is_open())
    {
        std::cout << "Failed to open " << filePath << "\n";
        return false;
    }

    long matIndex = 0;
    double x;

    while (file >> x)
    {
        S.arr[matIndex] = x;
        matIndex++;
        if (matIndex > N*N)
        {
            std::cout <<"File holds too many values for matirx\n";
            std::cout <<"count = " << matIndex << "\n";
            file.close();
            return false;
        }

    }

    if(matIndex != N*N)
    {
        std::cout << "File holds too little values for matrix\n";
        file.close();
        return false;
    }

    file.close();
    return true;
}

//-----------------------------------------------------------------------------
// Name : main ()
//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    std::cout << "the big jacobi test!" << std::endl;
    std::cout.precision(17);

    double* e = new double[N];
    matrix mat;
    matrix E;

//    mat[0][0] = 4;
//    mat[0][1] = -30;
//    mat[0][2] = 60;
//    mat[0][3] = -35;

//    mat[1][0] = -30;
//    mat[1][1] = 300;
//    mat[1][2] = -675;
//    mat[1][3] = 420;

//    mat[2][0] = 60;
//    mat[2][1] = -675;
//    mat[2][2] = 1620;
//    mat[2][3] = -1050;

//    mat[3][0] = -35;
//    mat[3][1] = 420;
//    mat[3][2] = -1050;
//    mat[3][3] = 700;

//    mat[0][0] = 4;
//    mat[0][1] = -30;
//    mat[0][2] = 60;
//    mat[0][3] = -35;

//    mat[1][0] = -30;
//    mat[1][1] = 300;
//    mat[1][2] = -675;
//    mat[1][3] = 420;

//    mat[2][0] = 60;
//    mat[2][1] = -675;
//    mat[2][2] = 1620;
//    mat[2][3] = -1050;

//    mat[3][0] = -35;
//    mat[3][1] = 420;
//    mat[3][2] = -1050;
//    mat[3][3] = 700;

//    mat[0][0] = 3;
//    mat[0][1] = 4;
//    mat[0][2] = 3;
//    mat[0][3] = 2;

//    mat[1][0] = 4;
//    mat[1][1] = 3;
//    mat[1][2] = 2;
//    mat[1][3] = 3;

//    mat[2][0] = 3;
//    mat[2][1] = 2;
//    mat[2][2] = 3;
//    mat[2][3] = 1;

//    mat[3][0] = 2;
//    mat[3][1] = 3;
//    mat[3][2] = 1;
//    mat[3][3] = 3;

    if (readMatrixFromFile(mat, "myFile4.txt") == false)
    {
        return 1;
    }

    std::cout <<"Finished loading matrix\n";

    double t = getTime();

    jacobi(mat, e, E);

    std::cout << "Time: " << getTime () - t << "\n\n";

//    for (int i = 0; i < N; i++)
//    {
//        std::cout << "e" << i+1 << " = " << e[i] << "\n\n";

//        for (int j = 0; j < N; j++)
//            std::cout << E[j][i] <<" ";

//        std::cout << "\n\n";
//    }

    delete[] e;
    return 0;
}
