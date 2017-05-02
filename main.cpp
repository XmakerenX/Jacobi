#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include<sys/time.h>
#include <pthread.h>
#include <unistd.h>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#define N 500

typedef void* (*THREADFUNC )(void* args );

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

struct Thread
{
    Thread (matrix& S, matrix& E, int* ind, int& k, int& l, double& c, double& s, THREADFUNC func)
        :S_(S), E_(E),k_(k), l_(l), c_(c), s_(s)
    {
        threadFunc_ = func;

        ind_ = ind;

        threadReady = false;
        threadDone = false;

        startIndex = 0;
        endIndex = 0;
    }

    ~Thread()
    {
        pthread_cancel(thread_);
    }

    bool start()
    {
        return (pthread_create(&thread_, NULL, threadFunc_, this) == 0);
    }

    void wait()
    {
        pthread_join(thread_, NULL);
    }

    THREADFUNC threadFunc_;

    pthread_t thread_;

    matrix& S_;
    matrix& E_;
    volatile int& k_;
    volatile int& l_;
    volatile double& c_;
    volatile double& s_;

    volatile int startIndex;
    volatile int endIndex;
    volatile bool threadReady;
    volatile bool threadDone;

    int* ind_;

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
// Name : toIterate ()
//-----------------------------------------------------------------------------
bool toIterate(matrix& S)
{
    bool con = false;

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            if (std::abs(S[i][j]) > 0.0001 )
            {
                con = true;
                return true;
            }
        }
    }

    return con;
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
//inline void rotate(matrix &S,int k, int l, int i, int j, double& c, double& s)
inline void rotate(matrix &S,int k, int l, int i, int j, volatile double& c,volatile double& s)
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
// Name : forRotate1 ()
//-----------------------------------------------------------------------------
void* forRotate1(void* arg)
{
    Thread* T = (Thread*)arg;

    while (1)
    {
        if (T->threadReady)
        {
            //std::cout <<"Thread 1 running\n";
            for (int i = T->startIndex; i < T->endIndex; i++)
            {
                rotate(T->S_, i, T->k_, i, T->l_, T->c_, T->s_);
            }
            T->threadReady = false;
            T->threadDone  = true;
            //std::cout << "threadDone1 = " << T->threadDone << "\n";
            //std::cout <<"Thread 1 Done\n";
        }
    }
}

//-----------------------------------------------------------------------------
// Name : forRotate2 ()
//-----------------------------------------------------------------------------
void* forRotate2(void* arg)
{
    Thread* T = (Thread*)arg;

    while (1)
    {
        if (T->threadReady)
        {
            //std::cout <<"Thread 2 running\n";
            for (int i = T->startIndex; i < T->endIndex; i++)
            {
                rotate(T->S_, T->k_, i, i, T->l_, T->c_, T->s_);
            }
            T->threadReady = false;
            T->threadDone=  true;
            //std::cout << "threadDone2 = " << T->threadDone << "\n";
            //std::cout <<"Thread 2 Done\n";
        }
    }
}

//-----------------------------------------------------------------------------
// Name : forRotate3 ()
//-----------------------------------------------------------------------------
void* forRotate3(void* arg)
{
    Thread* T = (Thread*)arg;

    while (1)
    {
        if (T->threadReady)
        {
            //std::cout <<"Thread 2 running\n";
            for (int i = T->startIndex; i < T->endIndex; i++)
            {
                rotate(T->S_, T->k_, i, T->l_, i, T->c_, T->s_);
            }
            T->threadReady = false;
            T->threadDone  =  true;
            //std::cout << "threadDone2 = " << T->threadDone << "\n";
            //std::cout <<"Thread 2 Done\n";
        }
    }
}

//-----------------------------------------------------------------------------
// Name : rotateVector ()
//-----------------------------------------------------------------------------
void* rotateVector(void* arg)
{
    Thread* T = (Thread*)arg;

    while (1)
    {
        if (T->threadReady)
        {
            for (int i = 0; i  < N; i++)
            {
                double newEik, newEil;

                newEik = T->c_ * T->E_[i][T->k_] - T->s_ * T->E_[i][T->l_];
                newEil = T->s_ * T->E_[i][T->k_] + T->c_ * T->E_[i][T->l_];

                T->E_[i][T->k_] = newEik;
                T->E_[i][T->l_] = newEil;
            }
            T->threadReady = false;
            T->threadDone  = true;

        }
    }
}

//-----------------------------------------------------------------------------
// Name : threadMaxind ()
//-----------------------------------------------------------------------------
void* threadMaxind(void* arg)
{
    Thread* T = (Thread*)arg;

    while (1)
    {
        if (T->threadReady)
        {
            T->ind_[T->k_] =  maxind(T->S_,T->k_);
            T->threadReady = false;
            T->threadDone  = true;

        }
    }
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

    double time1;
    double time2 = 0;
    double time3;
    Thread* threads[1];

    threads[0] = new Thread(S,E, ind,k, l ,c, s, rotateVector);
    //threads[1] = new Thread(S,E, ind,k, l ,c, s, threadMaxind);
    //threads[0] = new Thread(S, k, l ,c, s, forRotate2);
    //threads[1] = new Thread(S, k, l ,c, s, forRotate2);
    //threads[2] = new Thread(S, k, l ,c, s, forRotate3);

    threads[0]->start();
    //threads[1]->start();
//    threads[2]->start();

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            E[i][j] = 0;

    // E = I
    for (int i = 0; i < N; i++)
        E[i][i] = 1;

    for (int k = 0; k < N; k++)
    {
        ind[k] = maxind(S,k);
        e[k] = S[k][k];
        changed[k] = true;
    }

    while (toIterate(S))
    //while (state !=  0)
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

        threads[0]->threadReady = true;

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

        //t *= 0.1;

        S[k][l] = 0.0f;

        update(k, -t, state, e, changed);
        update(l, t, state, e, changed);

        // rotate rows and columns k and l
//        threads[0]->startIndex = 0;
//        threads[0]->endIndex = k;
//        threads[0]->threadReady = true;

        for (int i = 0; i < k; i++)
        {
            rotate(S,i,k,i,l, c, s);
        }
        //std::cout << "time1 = " << getTime() - time1 << "\n";
//        threads[0]->startIndex = k + 1;
//        threads[0]->endIndex = l;
//        threads[0]->threadReady = true;
        //time2 = getTime();
        for (int i = k +1; i < l; i++)
        {
            rotate(S,k,i,i,l, c, s);
        }
        //std::cout << "time2 = " << getTime() - time2 << "\n";

//        threads[2]->startIndex = l + 1;
//        threads[2]->endIndex = N;
//        threads[2]->threadReady = true;
        //time3 = getTime();
        for (int i = l + 1; i < N; i++)
        {
            rotate(S,k,i,l,i, c, s);
        }

        //threads[1]->threadReady = true;

        //std::cout << "time3 = " << getTime() - time3 << "\n";

//        for (int i = 0; i < k; i++)
//        {
//            rotate(S,i,k,i,l, c, s);
//        }

        // rotate eigenvectors
        //time1 = getTime();
//        for (int i = 0; i  < N; i++)
//        {
//            double eik, eil;
//            eik = c * E[i][k] - s * E[i][l];
//            eil = s * E[i][k] + c * E[i][l];

//            E[i][k] = eik;
//            E[i][l] = eil;
//        }
//        //std::cout << "time1 = " << getTime() - time1 << "\n";

        //std::cout << "waiting for threads\n";
       //time1 = getTime();
       //while (!threads[0]->threadDone);
       //while (!threads[0]->threadDone || !threads[1]->threadDone);
       //time2 += getTime() - time1;
       //while (!threads[0]->threadDone || !threads[1]->threadDone || !threads[2]->threadDone);
       //std::cout << "time1 = " << getTime() - time1 << "\n";
//       // std::cout << "Done waiting for threads\n";
       //threads[0]->threadDone = false;
       //threads[1]->threadDone = false;
       //threads[2]->threadDone = false;
       //std::cout << "Finished waiting threads\n";

        // rows k, l have changed, update rows indK, indL
        ind[k] =  maxind(S,k);
        ind[l] =  maxind(S,l);

        time1 = getTime();
        while (!threads[0]->threadDone);
        //while (!threads[0]->threadDone || !threads[1]->threadDone);
        //while (!threads[1]->threadDone);
        time2 += getTime() - time1;
        threads[0]->threadDone = false;

        count++;
    }

    delete[] ind;
    delete[] changed;

   std::cout << "count = "<< count << "\n\n";
   std::cout << "time2 = "<< time2 << "\n\n";
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
// Name : generateMatrix ()
//-----------------------------------------------------------------------------
void generateMatrix(matrix& mat, int n)
{
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            mat[i][j] = ( (rand() % 100) - 50);


    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << mat[i][j];
            if ( (j + 1) != n)
                std::cout << ",";
        }
        std::cout << "\n";
    }
}


//-----------------------------------------------------------------------------
// Name : main ()
//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    std::cout << "the big jacobi test!" << std::endl;
    std::cout.precision(16);

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

//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < N; j++)
//        {
//            std::cout << mat[i][j];
//            if ( (j + 1) != N)
//                std::cout << ",";
//        }
//        std::cout << "\n";
//    }

    //srand(unsigned(time(0)));
    //generateMatrix(mat , N);

    std::cout <<"Finished loading matrix\n";

    double t = getTime();

    jacobi(mat, e, E);

    std::cout << "Time: " << getTime () - t << "\n\n";

//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < N; j++)
//        {
//            std::cout << mat[i][j] << "\t";
//        }
//        std::cout << "\n";
//    }


//    for (int i = 0; i < N; i++)
//    {
//        std::cout << "e" << i+1 << " = " << e[i] << "\n\n";

////        for (int j = 0; j < N; j++)
////            std::cout << E[j][i] <<" ";

//        std::cout << "\n\n";
//    }

    delete[] e;
    return 0;
}
