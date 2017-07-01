#ifndef  _MATRIX_H
#define  _MATRIX_H

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <fstream>
#include <iomanip>

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
        :matrixPoints(newMatrixPoints) {}

    std::vector<matrixPoint> matrixPoints;
};

//-----------------------------------------------------------------------------
// matrix Class
//-----------------------------------------------------------------------------
class matrix
{

public:
    matrix  (int size);
    matrix  (std::string filePath);
    ~matrix ();

    // make accessing arr like accessing a matrix
    double* operator[](int i)
    {
        return arr + N*i;
    }

    double *            jacobi             ();
    double *            jacobiMkl          ();

    bool                readMatrixFromFile (std::string filePath);
    void                printMatrix        ();

    double *            getEigenvalue      ();
    matrix *            getEigenvalueVector();
    int                 getSize            ();

    std::vector<Kpairs> genrateMartixPoints();

private:
    // used by jacobi method
    void                update             (int k, double t, int& state, bool changed[]);
    inline void         rotate             (int k, int l, int i, int j,  double c, double s);
    // used by jacobiMkl method
    double              calcSum            ();
    bool                toIterate          (double exp);

    long                readMartixSize     (std::string filePath);

    double* arr;
    double* e;
    //matrix* E;
    std::vector<Kpairs> kPairs;
    int N; // matrix size
};

#endif  //_MATRIX_H
