#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include<sys/time.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>

#include <iomanip>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

//#include <math.h>

#include<mkl.h>

#define N 6

typedef void* (*THREADFUNC )(void* args );

struct matrix
{
    double* arr;

    matrix()
    {
        // allocate the matrix as one array
        //arr = new double[N*N];
        arr = (double *)mkl_malloc( N*N*sizeof( double ), 64 );
    }

    ~matrix()
    {
        mkl_free(arr);
        //delete[] arr;
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

bool print = false;

//-----------------------------------------------------------------------------
// Name : sigStop ()
// Desc : Handles the ctrl-z signal
//-----------------------------------------------------------------------------
void sigStop(int p)
{
    print = true;
}

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
// Name : printMatrix ()
//-----------------------------------------------------------------------------
void printMatrix(matrix& S)
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            std::cout << std::left << std::setw(17) << S[i][j] << "\t";

        std::cout << "\n";
    }

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

//    for (int i = k - 1; i >= 0; i--)
//    {
//        if (std::abs(S[k][i]) > std::abs(S[k][m]) )
//        {
//            m = i;
//        }
//    }

    return m;
}

//-----------------------------------------------------------------------------
// Name : calcSum ()
//-----------------------------------------------------------------------------
double calcSum(matrix& S)
{
    double sum = 0;
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
            sum += S[i][j];
    }

    return std::abs(sum);
}


//-----------------------------------------------------------------------------
// Name : toIterate ()
//-----------------------------------------------------------------------------
bool toIterate(matrix& S, double& exp)
{
    bool con = false;

    double sum = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
              sum += S[i][j] * S[i][j];
//            if (std::abs(S[i][j]) > 0.01 )
//            {
//                con = true;
//                return true;
//            }
        }

//        for (int j = i - 1; j >= 0; j--)
//        {
//              sum += S[i][j] * S[i][j];
////            if (std::abs(S[i][j]) > 0.01 )
////            {
////                con = true;
////                return true;
////            }
//        }
    }
    std::cout << "sum = " << std::sqrt(sum) << "\n";

//    for (int i = 0; i < N; i++)
//        std::cout << S[i][i] << " ";
//    std::cout <<"\n";

    if (std::sqrt(sum) < exp)
    //if (std::sqrt(sum) < 1)
        return false;
    else
        return true;

    //return con;
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
// Name : jacobiNeoMkl ()
//-----------------------------------------------------------------------------
void jacobiNeoMkl(matrix& S, double e[N], matrix& E)
{
    int i,k,l,m;
    double s,c,t,p,y,d,r,g;

    double esp;
    int sweepNum = (N * (N - 1)) / 2;

    int count = 0;

    double time1;
    double time2 = 0;
    double time3;

    matrix Qij;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            Qij[i][j] = 0;

    // E = I
    for (int i = 0; i < N; i++)
        Qij[i][i] = 1;

    esp = calcSum(S) * 0.0001;

    k = 0;
    l = 1;

    while (toIterate(S, esp))
    {
        //std::cout << "count = " << count << "\n";

        //std::cout << "k = "<< k << " l = "<< l <<  "\n";

        if (count == 45*45)
            break;

        printMatrix(S);
        std::cout <<"\n";

        std::cout << "S[k][k] = " << S[k][k] << "\n";
        std::cout << "S[l][l] = " << S[l][l] << "\n";
        std::cout << "S[k][l] = " << S[k][l] << "\n";

        g = 100.0*std::abs(S[k][l]);
        //y = (S[l][l] - S[k][k]) / (2 * S[k][l]);
        //y = 0.5 * (S[l][l] - S[k][k]) / (S[k][l]);
        //else
        {

//        y = 0.5 * (S[k][k] - S[l][l]) / (S[k][l]);

//        double h = S[k][k] - S[l][l];

//        t = 1 / (std::abs(y) + std::sqrt(1 + y*y));
//        if (y < 0)
//            t = -t;
////        double t1 = -y + std::sqrt(1 + y*y);
////        double t2 = -y - std::sqrt(1 + y*y);

////        if (std::abs(t1) < std::abs(t2))
////        {
////            std::cout << "t1\n";
////            t = t1;
////        }
////        else
////        {
////            std::cout << "t2\n";
////            t = t2;
////        }


//        c = 1 / (std::sqrt(1 + t*t) );

//        s = c * t;

        p = S[k][l];

        // calculate c = cos o, s = sin o<
        double pp = p*p;
        //std::cout << "(e[l] - e[k]) / 2 =" << (e[l] - e[k]) / 2 << "\n";
        y =  ( S[l][l] - S[k][k] ) / 2;
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

        std::cout << "g = " << g << "\n";
        std::cout << "y = " << y << "\n";
        std::cout << "t = " << t << "\n";
        std::cout << "c = " << c << "\n";
        std::cout << "s = " << s << "\n";

        Qij[k][k] = c;
        Qij[l][l] = c;
        Qij[l][k] = -s;
        Qij[k][l] = s;

//        printMatrix(Qij);
//        std::cout << "\n";

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                     N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, S.arr, N);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                     N, N, N, 1.0, S.arr, N, Qij.arr, N, 0.0, S.arr, N);


        S[k][l] = 0.0;
        S[l][k] = 0.0;

        //Restore Qij to being I matrixs
        Qij[k][k] = 1;
        Qij[l][l] = 1;
        Qij[l][k] = 0;
        Qij[k][l] = 0;

        }

        std::cout << "k = "<< k << " l = "<< l <<  "\n";

        if ( (k < N - 2) && l < N - 1)
        {
            k = k;
            l = l +1;
        }
        else
        {
            if (k < N - 2 && l == N - 1)
            {
                int temp = k;
                k = temp + 1;
                l = temp + 2;
            }
            else
            {
                if ( (k == N - 2) && (l == N - 1))
                {
                    k = 0;
                    l = 1;
                }
                else
                    std::cout <<"Error no condition was entered\n";
            }
        }

        count++;

    }

    for (int i = 0; i < N; i++)
    {
        e[i] = S[i][i];
    }

   std::cout << "count = "<< count << "\n\n";
   std::cout << "time2 = "<< time2 << "\n\n";
}

//-----------------------------------------------------------------------------
// Name : jacobiMkl ()
//-----------------------------------------------------------------------------
void jacobiMkl(matrix& S, double e[N], matrix& E)
{
    int i,k,l,m;
    double s,c,t,p,y,d,r;

    double esp;
    int sweepNum = (N * (N - 1)) / 2;

    matrix Qij;

    int state = N;

    int* ind = new int[N];
    bool* changed = new bool[N];

    int count = 0;

    double time1;
    double time2 = 0;
    double time3;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            E[i][j] = 0;

    // E = I
    for (int i = 0; i < N; i++)
        E[i][i] = 1;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            Qij[i][j] = E[i][j];
        }

//    for (int k = 0; k < N; k++)
//    {
//        e[k] = S[k][k];
//        changed[k] = true;
//    }

    esp = calcSum(S) * 0.0001;

    std::cout << "esp = " << esp << "\n";

    k = 0;
    l = 1;

    while (toIterate(S, esp))
    {
//        m = 0;
//        for (k = 1; k < N-1; k++)
//        {
//            if (std::abs(S[k][ind[k]]) > std::abs(S[m][ind[m]]) )
//            {
//                m = k;
//            }
//        }

//        k = m;
//        l = ind[m];


        for (int j = 0; j < sweepNum; j++)
        {
            if (count == 1300)
                break;

            if (print)
            {

                for (int i = 0; i < N; i++)
                {
                    e[i] = S[i][i];
                }

                bool swapped = false;
                do
                {
                    swapped = false;
                    for (int i = 1; i < N; i++)
                    {
                        if (e[i-1] > e[i])
                        {
                            double temp = e[i -1];
                            e[i - 1] = e[i];
                            e[i] = temp;
                            swapped = true;
                        }
                    }
                } while (swapped);

                for (int i = 0; i < N; i++)
                {
                    std::cout << "e" << i+1 << " = " << e[i] << "\n\n";

                    std::cout << "\n\n";
                }

                print = false;
            }

            double delta = esp / sweepNum;
            if (delta <= std::sqrt ((S[k][l] * S[k][l]) + (S[l][k] * S[l][k])))
            {
        //std::cout << "k = " << k << "\n";
        //std::cout << "l = " << l << "\n";
        //double beta = (S[k][k] - S[l][l]) / (2 * S[k][l]);
        //double beta = (S[k][k] - S[l][l]) / (2 * S[k][l]);
       // double beta = (S[l][l] - S[k][k]) / (2 * S[k][l]);

//        double cot;
//        double cot1 = beta - std::sqrt((beta* beta) + 1);
//        double cot2 = beta + std::sqrt((beta* beta) + 1);

//        if (std::abs(cot1) > std::abs(cot2))
//            cot = cot1;
//        else
//            cot = cot2;

//        std::cout << "beta = " << beta << "\n";
//        std::cout << "t = " << t << "\n";

//        std::cout << "beta = " << beta << "\n";
//        std::cout << "cot = " << cot << "\n";

//        s = 1 / (std::sqrt(cot * cot + 1));
//        c = s * cot;

//        y = (S[l][l] - S[k][k]) / (2 * S[k][l]);
//        if (y > 0)
//        {
//            t = -y + std::sqrt(y*y + 1);
//        }
//        else
//        {
//            t = -y - std::sqrt(y*y + 1);
//        }

//        c = 1 / std::sqrt(1 + t*t);
//        s = c * t;

          //y = (-2 * S[k][l] ) / (S[l][l]- S[k][k]);
          y = (S[l][l]- S[k][k]) / (-2 * S[k][l] );
          if (y > 0)
              t = 1 / (std::abs(y) + std::sqrt(1 + y*y));
          else
              t = -1 / (std::abs(y) + std::sqrt(1 + y*y));

          c = 1 / std::sqrt(1 + t*t);
          s = c * t;


//        p = S[k][l];

//        // calculate c = cos o, s = sin o
//        double pp = p*p;
//        y =  ( S[l][l] - S[k][k] ) / 2;
//        d = std::abs(y) + std::sqrt(pp + y*y);
//        r = std::sqrt(pp + d*d);
//        c = d / r;
//        s = p / r;
//        t = pp / d;

//        if ( y < 0)
//        {
//            s = -s;
//            t = -t;
//        }

//          std::cout << "c = " << c << "\n";
//          std::cout << "s = " << s << "\n";
//          std::cout << "l = " << l << "\n";
//          std::cout << "k = " << k << "\n";

//        if (l > k)
//        {
//            s = -s;
//        }

        Qij[k][k] = c;
        Qij[l][l] = c;
        Qij[l][k] = -s;
        Qij[k][l] = s;


        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                     N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, S.arr, N);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                     N, N, N, 1.0, S.arr, N, Qij.arr, N, 0.0, S.arr, N);

        //time1 = getTime();
        // S = Qij * S
//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                     N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, S.arr, N);

//        // S = S * (Qij)^T
//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//                     N, N, N, 1.0, S.arr, N, Qij.arr, N, 0.0, S.arr, N);


        //std::cout << "time1 = " << getTime() - time1 << "\n";
//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans,
//                     N, N, N, 1.0, S.arr, N, Qij.arr, N, 0.0, S.arr, N);

//        for (int i = 0; i < N; i++)
//        {
//            for (int j = 0; j < N; j++)
//                std::cout << Qij[i][j] << " ";

//            std::cout << "\n";
//        }

//        std::cout <<"###################################\n";
//        for (int i = 0; i < N; i++)
//        {
//            for (int j = 0; j < N; j++)
//                std::cout << Qij[i][j] << " ";

//            std::cout << "\n";
//        }
//        std::cout <<"###################################\n";

        //Restore Qij to being I matrixs
        Qij[k][k] = 1;
        Qij[l][l] = 1;
        Qij[l][k] = 0;
        Qij[k][l] = 0;

        //S[k][l] = 0;
//        std::cout <<"-----------------------------------\n";
//        for (int i = 0; i < N; i++)
//        {
//            for (int j = 0; j < N; j++)
//                std::cout << S[i][j] << "\t";

//            std::cout << "\n";
//        }

//        std::cout <<"-----------------------------------\n";

        //ind[k] =  maxind(S,k);
        //ind[l] =  maxind(S,l);
        }
        if ( (k < N - 2) && l < N - 1)
        {
            k = k;
            l = l +1;
        }
        else
        {
            if (k < N - 2 && l == N - 1)
            {
                int temp = k;
                k = temp + 1;
                l = temp + 2;
            }
            else
            {
                if ( (k == N - 2) && (l == N - 1))
                {
                    k = 0;
                    l = 1;
                }
                else
                    std::cout <<"Error no condition was entered\n";
            }
        }


        count++;
    }
  }

    for (int i = 0; i < N; i++)
    {
        e[i] = S[i][i];
    }

    delete[] ind;
    delete[] changed;

   std::cout << "count = "<< count << "\n\n";
   std::cout << "time2 = "<< time2 << "\n\n";
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

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            E[i][j] = 0;

    // E = I
    for (int i = 0; i < N; i++)
        E[i][i] = 1;

    for (int k = 0; k < N; k++)
    {
        //ind[k] = maxind(S,k);
        e[k] = S[k][k];
        changed[k] = true;
    }

    k = 0;
    l = 1;

    while (state !=  0)
    {
//        if (count == 45*4)
//            break;

          printMatrix(S);
          std::cout << "\n";
          toIterate(S,s);
//        m = 0;
//        for (k = 1; k < N-1; k++)
//        {
//            if (std::abs(S[k][ind[k]]) > std::abs(S[m][ind[m]]) )
//            {
//                m = k;
//            }
//        }

//        k = m;
//        l = ind[m];



//        std::cout << "k = " << k << " l = " << l << "\n";
//        for (int i = 0; i < N; i++)
//            std::cout << "ind " << i << " = "<< ind[i] << " \n";
        //std::cout <<"\n";

//        y = (S[l][l] - S[k][k]) / (2 * S[k][l]);
//        if (y > 0)
//        {
//            t = -y + std::sqrt(y*y + 1);
//        }
//        else
//        {
//            t = -y - std::sqrt(y*y + 1);
//        }

//        c = 1 / std::sqrt(1 + t*t);
//        s = c * t;

        p = S[k][l];

        // calculate c = cos o, s = sin o<
        double pp = p*p;
        //std::cout << "(e[l] - e[k]) / 2 =" << (e[l] - e[k]) / 2 << "\n";
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

        //std::cout << "t = " << t << "\n";
        std::cout << "c = " << c << "\n";
        std::cout << "s = " << s << "\n";

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

//        std::cout <<"---------------------------------------------\n";
//        printMatrix(S);
//        std::cout <<"---------------------------------------------\n";

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

        // rows k, l have changed, update rows indK, indL
        //ind[k] =  maxind(S,k);
        //ind[l] =  maxind(S,l);

        std::cout << "k = "<< k << " l = " << l << "\n";

        if ( (k < N - 2) && l < N - 1)
        {
            k = k;
            l = l +1;
        }
        else
        {
            if (k < N - 2 && l == N - 1)
            {
                int temp = k;
                k = temp + 1;
                l = temp + 2;
            }
            else
            {
                if ( (k == N - 2) && (l == N - 1))
                {
                    k = 0;
                    l = 1;
                }
                else
                    std::cout <<"Error no condition was entered\n";
            }
        }

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
    signal(SIGTSTP, &sigStop);

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

    //if (readMatrixFromFile(mat, "5002.txt") == false)
    //if (readMatrixFromFile(mat, "myMat2.txt") == false)
    //if (readMatrixFromFile(mat, "424.txt") == false)
    //if (readMatrixFromFile(mat, "500b.txt") == false)
    //if (readMatrixFromFile(mat, "102b.txt") == false)
    //if (readMatrixFromFile(mat, "222.txt") == false)
    if (readMatrixFromFile(mat, "6b.txt") == false)
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

    jacobiNeoMkl(mat, e, E);
    //jacobiMkl(mat, e, E);
    //jacobi(mat, e, E);

    std::cout << "Time: " << getTime () - t << "\n\n";

    bool swapped = false;
    do
    {
        swapped = false;
        for (int i = 1; i < N; i++)
        {
            if (e[i-1] > e[i])
            {
                double temp = e[i -1];
                e[i - 1] = e[i];
                e[i] = temp;
                swapped = true;
            }
        }
    } while (swapped);

//    for (int i = 0; i < N; i++)
//    {
//        for (int j = 0; j < N; j++)
//        {
//            std::cout << mat[i][j] << "\t";
//        }
//        std::cout << "\n";
//    }


    for (int i = 0; i < N; i++)
    {
        std::cout << "e" << i+1 << " = " << e[i] << "\n\n";

//        for (int j = 0; j < N; j++)
//            std::cout << E[j][i] <<" ";

        std::cout << "\n\n";
    }

    delete[] e;
    return 0;
}
