#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sys/time.h>
#include <pthread.h>
#include <unistd.h>
#include <signal.h>

#include <iomanip>

#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

#include<mkl.h>

#define N 8

typedef void* (*THREADFUNC )(void* args );
// if to print current results
bool print = false;

//-----------------------------------------------------------------------------
// matrix Class
//-----------------------------------------------------------------------------
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

//-----------------------------------------------------------------------------
// Thread Class
//-----------------------------------------------------------------------------
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

struct matrixPoint
{
    matrixPoint(int newK, int newL)
    {
        k = newK;
        l = newL;
    }

    int k;
    int l;
};

struct Kpairs
{
    Kpairs(std::vector<matrixPoint> newMatrixPoints)
        :matrixPoints(newMatrixPoints)
    {
        ;
    }

    std::vector<matrixPoint> matrixPoints;
};

//-----------------------------------------------------------------------------
// Name : genrateMartixPoints ()
// Desc : generate Martrix pairs to iterate on
//-----------------------------------------------------------------------------
std::vector<Kpairs> genrateMartixPoints(int n)
{
    int m = (n + 1)/2;

    int q,p;
    int c = 1;
    std::vector<matrixPoint> vecPoints;
    std::vector<Kpairs> vecKPairs;

    for (int k = 1; k <= m - 1; k++)
    {
        for (int c = 1; (c + m) <= n; c++)
        {
            q = m - k + c;

            if ( ((m - k + 1) <= q) &&  (q <= (2*m - 2*k)) )
            {
                p = (2*m - 2*k + 1) - q;
            }
            else
                if ( (2*m - 2*k < q) && (q <= (2*m - k - 1)) )
                {
                    p = (4*m - 2*k ) - q;
                }
                else
                    if ( (2*m - k - 1) < q)
                        p = n;


            if (p > q)
            {
                int temp = q;
                q  = p;
                p = temp;
            }

            vecPoints.push_back( matrixPoint(p - 1 , q - 1) );
            //std::cout << "(" << p << ", " << q <<  ")\n";
        }
        vecKPairs.push_back(Kpairs(vecPoints) );
        vecPoints.clear();
    }

    for (int k = m; k <= (2*m - 1); k++)
    {

        for (int c = 0; c < (n -m); c++)
        {
            q = 4*m - n -k + c;

            if (( q < (2*m - k + 1)) )
                p = n;
            else
                if ( ((2*m - k + 1) <= q) &&  (q <=(4*m - 2*k - 1)) )
                    p = ((4*m - 2*k) - q);
                else
                    if ((4*m - 2*k - 1 )< q)
                        p = (6*m - 2*k -1) - q;

            if (p > q)
            {
                int temp = q;
                q  = p;
                p = temp;
            }

            vecPoints.push_back( matrixPoint(p - 1 , q - 1) );
            //std::cout << "(" << p << ", " << q <<  ")\n";
        }
        vecKPairs.push_back(Kpairs(vecPoints) );
        vecPoints.clear();
    }

    return vecKPairs;

}

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
void jacobiNeoMkl(matrix& S, double e[N], matrix& E, std::vector<Kpairs>& kPairs)
{
    int i,k,l,m;
    double s,c,t,p,y,d,r,g;

    double esp;
    int sweepNum = (N * (N - 1)) / 2;

    int count = 0;

    double time1;
    double time2 = 0;
    double time3 = 0;

    matrix Qij;
    matrix QijT;
    matrix temp;

    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
        {
            Qij[i][j] = 0;   
            QijT[i][j] = 0;
        }

    // E = I
    for (int i = 0; i < N; i++)
    {
        Qij[i][i] = 1;
        QijT[i][i] = 1;
    }

    esp = calcSum(S) * 0.001;

    k = 0;
    l = 1;

    bool end = true;

    //while (toIterate(S, esp))
    while(end)
    {

        time1 = getTime();

        if (!toIterate(S, esp))
            break;

        for (int i = 0; i <kPairs.size(); i++)
        {
            Kpairs pair = kPairs[i];

            if (print)
            {
                toIterate(S, s);
                print = false;
            }

            //std::cout << "count = " << count << "\n";
//            if (count == 3992)
//            {
//                std::cout << "ending by count\n";
//                end = false;
//                break;
//            }

            for (int j = 0; j < pair.matrixPoints.size(); j++)
            {
                k = pair.matrixPoints[j].k;
                l = pair.matrixPoints[j].l;

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

                Qij[k][k] = c;
                Qij[l][l] = c;
                Qij[l][k] = -s;
                Qij[k][l] = s;

                QijT[k][k] = c;
                QijT[l][l] = c;
                QijT[k][l] = -s;
                QijT[l][k] = s;
            }

//            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
//                         N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, S.arr, N);

//            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                         N, N, N, 1.0, S.arr, N, Qij.arr, N, 0.0, S.arr, N);

            //time1 = getTime();
            time3 += getTime() - time1;
//            cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
//                         N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, temp.arr, N);
            cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, N, N, 1.0, S.arr, N, QijT.arr, N, 0.0, temp.arr, N );

            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                         N, N, N, 1.0, temp.arr, N, Qij.arr, N, 0.0, S.arr, N);

            time1 = getTime();
            //std::cout << "11111111\n";
            //cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, temp.arr, N );
            //std::cout << "22222222\n";
            //cblas_dsymm(CblasRowMajor, CblasLeft, CblasUpper, N, N, 1.0, temp.arr, N, QijT.arr, N, 0.0, S.arr, N );
            //std::cout << "33333333\n";

            //std::cout <<"time1 = " << getTime() - time1 << "\n";
            //time2 += getTime() - time1;

            for (int j = 0; j < pair.matrixPoints.size(); j++)
            {
                k = pair.matrixPoints[j].k;
                l = pair.matrixPoints[j].l;

                Qij[k][k] = 1;
                Qij[l][l] = 1;
                Qij[l][k] = 0;
                Qij[k][l] = 0;

                QijT[k][k] = 1;
                QijT[l][l] = 1;
                QijT[l][k] = 0;
                QijT[k][l] = 0;
            }

            count++;

            time2 += getTime() - time1;
        }
    }

    for (int i = 0; i < N; i++)
    {
        e[i] = S[i][i];
    }

   std::cout << "count = "<< count << "\n\n";
   std::cout << "time2 = "<< time2 << "\n\n";
   std::cout << "time3 = "<< time3 << "\n\n";
}

//-----------------------------------------------------------------------------
// Name : jacobiMkl ()
//-----------------------------------------------------------------------------
void jacobiMkl(matrix& S, double e[N], matrix& E)
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
    matrix temp;

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

//        if (count == 45*45)
//            break;

//        printMatrix(S);
//        std::cout <<"\n";

//        std::cout << "S[k][k] = " << S[k][k] << "\n";
//        std::cout << "S[l][l] = " << S[l][l] << "\n";
//        std::cout << "S[k][l] = " << S[k][l] << "\n";

        //g = 100.0*std::abs(S[k][l]);

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

//        std::cout << "g = " << g << "\n";
//        std::cout << "y = " << y << "\n";
//        std::cout << "t = " << t << "\n";
//        std::cout << "c = " << c << "\n";
//        std::cout << "s = " << s << "\n";

        Qij[k][k] = c;
        Qij[l][l] = c;
        Qij[l][k] = -s;
        Qij[k][l] = s;

//        printMatrix(Qij);
//        std::cout << "\n";

//        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
//                     N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, S.arr, N);

//        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
//                     N, N, N, 1.0, S.arr, N, Qij.arr, N, 0.0, S.arr, N);

        time1 = getTime();
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans,
                     N, N, N, 1.0, Qij.arr, N, S.arr, N, 0.0, temp.arr, N);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,
                     N, N, N, 1.0, temp.arr, N, Qij.arr, N, 0.0, S.arr, N);

        time2 += getTime() - time1;

        //Restore Qij to being I matrixs
        Qij[k][k] = 1;
        Qij[l][l] = 1;
        Qij[l][k] = 0;
        Qij[k][l] = 0;



        //std::cout << "k = "<< k << " l = "<< l <<  "\n";

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
// Name : jacobi ()
//-----------------------------------------------------------------------------
void jacobi(matrix& S, double e[N], matrix& E)
{
    int i,k,l,m;
    double s,c,t,p,y,d,r;

    int state = N;

    int* ind = new int[N];
    bool* changed = new bool[N];
    double* e2 = new double[N];

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
        if (print)
        {
            toIterate(S, s);
            print = false;
        }
//        if (count == 45*4)
//            break;

/*         printMatrix(S);
         std::cout << "\n"*/;

        p = S[k][l];

        // calculate c = cos o, s = sin o<
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

        //std::cout << "t = " << t << "\n";
        //std::cout << "c = " << c << "\n";
        //std::cout << "s = " << s << "\n";

        S[k][l] = 0.0f;

        update(k, -t, state, e, changed);
        update(l, t, state, e, changed);

        std::cout << "k = " << k << " l = " << l << "\n";

//        std::cout << "Start: for (int i = 0; i < k; i++) \n";
        // rotate rows and columns k and l
        for (int i = 0; i < k; i++)
        {
            std::cout << "("<< i << ", " << k << ")\n";
            std::cout << "("<< i << ", " << l << ")\n";
            rotate(S,i,k,i,l, c, s);
        }

        //std::cout << "End: for (int i = 0; i < k; i++) \n";

        //std::cout << "Start: for (int i = k +1; i < l; i++)\n";

        for (int i = k +1; i < l; i++)
        {
            std::cout << "("<< k << ", " << i << ")\n";
            std::cout << "("<< i << ", " << l << ")\n";
            rotate(S,k,i,i,l, c, s);
        }

        //std::cout << "End: for (int i = k +1; i < l; i++)\n";

        //std::cout << "Start: for (int i = l + 1; i < N; i++)\n";

        for (int i = l + 1; i < N; i++)
        {
            std::cout << "("<< k << ", " << i << ")\n";
            std::cout << "("<< l << ", " << i << ")\n";
            rotate(S,k,i,l,i, c, s);
        }


        //std::cout << "End: for (int i = l + 1; i < N; i++)\n";

//        for (int i = 0; i < N; i++)
//        {
//            e2[i] = e[i];
//        }

//        bool swapped = false;
//        do
//        {
//            swapped = false;
//            for (int i = 1; i < N; i++)
//            {
//                if (e2[i-1] > e2[i])
//                {
//                    double temp = e2[i -1];
//                    e2[i - 1] = e2[i];
//                    e2[i] = temp;
//                    swapped = true;
//                }
//            }
//        } while (swapped);

//        for (int i = 0; i < N; i++)
//        {
//            std::cout << "e" << i+1 << " = " << e[i] << "\n";

//            std::cout << "\n\n";
//        }

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

//        std::cout << "k = "<< k << " l = " << l << "\n";

        time1 = getTime();
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
        time2 += getTime() - time1;

        count++;
    }

    delete[] ind;
    delete[] changed;

   std::cout << "count = "<< count << "\n\n";
   std::cout << "time2 = "<< time2 << "\n\n";
}

//-----------------------------------------------------------------------------
// Name : neoJacobi ()
//-----------------------------------------------------------------------------
void neoJacobi(matrix& S, double e[N], matrix& E, std::vector<Kpairs>& kPairs)
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
        e[k] = S[k][k];
        changed[k] = true;
    }

    while (state !=  0)
    {
//        if (count == 45*4)
//            break;

        for (int i = 0; i < kPairs.size(); i++)
        {
            Kpairs pair = kPairs[i];
            for (int j = 0; j < pair.matrixPoints.size(); j++)
            {
                k = pair.matrixPoints[j].k;
                l = pair.matrixPoints[j].l;

                p = S[k][l];

                // calculate c = cos o, s = sin o<
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

    //if (readMatrixFromFile(mat, "myMat2.txt") == false)
    //if (readMatrixFromFile(mat, "424.txt") == false)
    //if (readMatrixFromFile(mat, "500b.txt") == false)
    //if (readMatrixFromFile(mat, "102b.txt") == false)
    //if (readMatrixFromFile(mat, "222.txt") == false)
    if (readMatrixFromFile(mat, "08b.txt") == false)
    //if (readMatrixFromFile(mat, "6b.txt") == false)
    //if (readMatrixFromFile(mat, "5002.txt") == false)
    //if (readMatrixFromFile(mat, "100b.txt") == false)
    //if (readMatrixFromFile(mat, "5000d.txt") == false)
    //if (readMatrixFromFile(mat, "3000b.txt") == false)
    //if (readMatrixFromFile(mat, "2500b.txt") == false)
    //if (readMatrixFromFile(mat, "2000b.txt") == false)
    //if (readMatrixFromFile(mat, "1500b.txt") == false)
    //if (readMatrixFromFile(mat, "1000b.txt") == false)
    //if (readMatrixFromFile(mat, "900b.txt") == false)
    //if (readMatrixFromFile(mat, "800b.txt") == false)
    //if (readMatrixFromFile(mat, "ert4.txt") == false)
    //if (readMatrixFromFile(mat, "ert2.txt") == false)
    //if (readMatrixFromFile(mat, "500d.txt") == false)
    //if (readMatrixFromFile(mat, "450b.txt") == false)
    //if (readMatrixFromFile(mat, "400b.txt") == false)
    //if (readMatrixFromFile(mat, "395b.txt") == false)
    //if (readMatrixFromFile(mat, "390b.txt") == false)
    //if (readMatrixFromFile(mat, "375b.txt") == false)
    //if (readMatrixFromFile(mat, "350b.txt") == false)
    //if (readMatrixFromFile(mat, "300b.txt") == false)
    //if (readMatrixFromFile(mat, "200b.txt") == false)
    //if (readMatrixFromFile(mat, "150b.txt") == false)
    //if (readMatrixFromFile(mat, "125b.txt") == false)
    //if (readMatrixFromFile(mat, "105b.txt") == false)
    //if (readMatrixFromFile(mat, "101b.txt") == false)
    {
        return 1;
    }

    std::vector<Kpairs> kPairs = genrateMartixPoints(N);

    for (int i = 0; i < kPairs.size(); i++)
    {
        std::vector<matrixPoint>& matrixPoints = kPairs[i].matrixPoints;

        for (int j = 0; j < matrixPoints.size(); j++)
        {
            std::cout << matrixPoints[j].k << " " <<  matrixPoints[j].l << "\n";
        }
        std::cout<<"---------------------\n";
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

    //jacobiNeoMkl(mat, e, E,kPairs);
    //jacobiMkl(mat, e, E);
    jacobi(mat, e, E);
    //neoJacobi(mat, e, E, kPairs);

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
