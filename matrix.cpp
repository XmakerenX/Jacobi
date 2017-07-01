#include"main.h"
#include"matrix.h"
#include<mkl.h>
#include<sstream>

//-----------------------------------------------------------------------------
// Constructor
// Desc : creates emtpy matrix with the given size
//-----------------------------------------------------------------------------
matrix::matrix(int size)
{
    N = size;
    // allocate the matrix as one array
    arr = (double *)mkl_malloc( N*N*sizeof( double ), 64 );
    e = new double[N];
    kPairs = genrateMartixPoints();
    //E = new matrix(N);
}

//-----------------------------------------------------------------------------
// Constructor
// Desc : create new matrix and init it from the given file
//-----------------------------------------------------------------------------
matrix::matrix(std::string filePath)
{
    N = readMartixSize(filePath);
    if (N == 0)
        std::exit(1);
    // allocate the matrix as one array
    arr = (double *)mkl_malloc( N*N*sizeof( double ), 64 );
    e = new double[N];
    if (!readMatrixFromFile(filePath))
        std::exit(1);
    kPairs = genrateMartixPoints();
    //E = new matrix(N);
}

//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------
matrix::~matrix()
{
    mkl_free(arr);
    delete[] e;
    //delete E;
}

//-----------------------------------------------------------------------------
// Name : update ()
// Desc : update e[k] to its status
//-----------------------------------------------------------------------------
void matrix::update(int k, double t, int& state, bool changed[])
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
inline void matrix::rotate(int k, int l, int i, int j,  double c, double s)
{
    // | S(kl) |   | c -s | |S(kl) |
    // |       | = |      | |      |
    // | S(ij) |   | s  c | |S(ij) |
    double skl,sij;

    matrix& S = *this;

    skl = c * S[k][l] -s * S[i][j];
    sij = s * S[k][l] +c * S[i][j];

    S[k][l] = skl;
    S[i][j] = sij;
}

//-----------------------------------------------------------------------------
// Name : jacobi ()
//-----------------------------------------------------------------------------
double* matrix::jacobi()
{
    int i,k,l,m;
    double s,c,t,p,y,d,r;

    int state = N;

    bool* changed = new bool[N];

    int count = 0;

    matrix& S = *this;

//        for (int i = 0; i < N; i++)
//            for (int j = 0; j < N; j++)
//                E[i][j] = 0;

//        // E = I
//        for (int i = 0; i < N; i++)
//            E[i][i] = 1;

    for (int k = 0; k < N; k++)
    {
        e[k] = S[k][k];
        changed[k] = true;
    }

    // starting k and  l
    k = 0;
    l = 1;

    while (state !=  0)
    {
        // choose the current Skl to be the pivot
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

        // the calc should cause Skl to be zero (or close to it)
        S[k][l] = 0.0f;
        // update e[k] and e[l] and check state
        update(k, -t, state, changed);
        update(l, t, state, changed);

        // rotate rows and columns k and l
        for (int i = 0; i < k; i++)
            rotate(i,k,i,l, c, s);

        for (int i = k +1; i < l; i++)
            rotate(k,i,i,l, c, s);

        for (int i = l + 1; i < N; i++)
            rotate(k,i,l,i, c, s);


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

        // same as doing
        // for (k = 0; k < N;  k++)
        //  for (l = k + 1; k < N; k++)
        if ( (k < N - 2) && l < N - 1)
            l = l +1;
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

    delete[] changed;

    return e;
}

//-----------------------------------------------------------------------------
// Name : jacobiMkl ()
//-----------------------------------------------------------------------------
double* matrix::jacobiMkl()
{
    int i,k,l,m;
    double s,c,t,p,y,d,r,g;

    double esp;
    int sweepNum = (N * (N - 1)) / 2;

    matrix& S = *this;

    matrix Qij(N);
    matrix QijT(N);
    matrix temp(N);

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

    esp = calcSum() * 0.001;

    k = 0;
    l = 1;

    bool end = true;

    while (toIterate(esp))
    {
        for (int i = 0; i <kPairs.size(); i++)
        {
            Kpairs pair = kPairs[i];
            for (int j = 0; j < pair.matrixPoints.size(); j++)
            {
                k = pair.matrixPoints[j].k;
                l = pair.matrixPoints[j].l;

                p = S[k][l];

                // calculate c = cos o, s = sin o
                double pp = p*p;
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

            cblas_dsymm(CblasRowMajor, CblasRight, CblasUpper, N, N, 1.0, S.arr, N, QijT.arr, N, 0.0, temp.arr, N );
            cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, N, N, N, 1.0, temp.arr, N, Qij.arr, N, 0.0, S.arr, N);

            // restore Qij to unit matrix
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
        }
    }

    for (int i = 0; i < N; i++)
    {
        e[i] = S[i][i];
    }

    return e;
}

//-----------------------------------------------------------------------------
// Name : readMartixSize ()
//-----------------------------------------------------------------------------
long matrix::readMartixSize(std::string filePath)
{
    std::ifstream file(filePath);

    std::string line;

    if (!file.is_open())
    {
        std::cout << "Failed to open "<< filePath << "\n";
        return false;
    }

    long rowSize = 1;
    long colSize = 0;
    double x;

    // read first row
    std::getline(file,line);
    // make sure it read something
    if (line == "")
        if (rowSize != colSize)
        {
            std::cout << "Error: Matrix is not squared\n";
            return 0;
        }

    // count how many numbers in the first row
    std::stringstream strStream(line);
    do
    {
        std::getline(strStream, line , ',');
        colSize++;

    }while(!strStream.eof());

    // count how many columns the matrix have
    do
    {
        std::getline(file,line);
        // make sure it read something
        if (line == "")
                break;

        rowSize++;

    }while(!file.eof());

    file.close();

    // return error when given non-squared matrix
    if (rowSize != colSize)
    {
        std::cout << "Error: Matrix is not squared\n";
        return 0;
    }
    else
        return rowSize;
}

//-----------------------------------------------------------------------------
// Name : readMatrixFromFile ()
//-----------------------------------------------------------------------------
bool matrix::readMatrixFromFile(std::string filePath)
{
    matrix& S = *this;

    std::ifstream file(filePath);
    //std::stringstream strStream;
    std::string str;

    if (!file.is_open())
    {
        std::cout << "Failed to open " << filePath << "\n";
        return false;
    }

    long matIndex = 0;
    double x;

    do
    {
        std::getline(file,str);
        // nothing was read than finish reading
        if (str == "")
            break;

        std::stringstream strStream(str);
        do
        {
            std::getline(strStream, str, ',');
            std::stringstream strStreamComa(str);
            strStreamComa >> x;

            S.arr[matIndex] = x;
            matIndex++;
            if (matIndex > N*N)
            {
                std::cout <<"File holds too many values for matirx\n";
                std::cout <<"count = " << matIndex << "\n";
                file.close();
                return false;
            }

        }while(!strStream.eof());

    }while(!file.eof());

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
// Name : printMatrix ()
//-----------------------------------------------------------------------------
void matrix::printMatrix()
{
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
            std::cout << std::left << std::setw(17) << *this[i][j] << "\t";

        std::cout << "\n";
    }

}

//-----------------------------------------------------------------------------
// Name : calcSum ()
//-----------------------------------------------------------------------------
double matrix::calcSum()
{
    matrix& S = *this;

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
bool matrix::toIterate(double exp)
{
    matrix& S = *this;
    bool con = false;

    double sum = 0;

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
              sum += S[i][j] * S[i][j];
        }
    }

    if (std::sqrt(sum) < exp)
    //if (std::sqrt(sum) < 1)
        return false;
    else
        return true;
}

//-----------------------------------------------------------------------------
// Name : genrateMartixPoints ()
// Dsec : genrate a vecotr of n/2 points that jacboi rotation  can be done on
//        them in one go
// **** : only needed for the mkl implementation
//-----------------------------------------------------------------------------
std::vector<Kpairs> matrix::genrateMartixPoints()
{
    int m = (N + 1)/2;

    int q,p;
    int c = 1;
    std::vector<matrixPoint> vecPoints;
    std::vector<Kpairs> vecKPairs;

    for (int k = 1; k <= m - 1; k++)
    {
        for (int c = 1; (c + m) <= N; c++)
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
                        p = N;


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
        for (int c = 0; c < (N -m); c++)
        {
            q = 4*m - N -k + c;

            if (( q < (2*m - k + 1)) )
                p = N;
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
// Name : getEigenvalue ()
//-----------------------------------------------------------------------------
double* matrix::getEigenvalue()
{
    return e;
}

//-----------------------------------------------------------------------------
// Name : getEigenvalueVector ()
//-----------------------------------------------------------------------------
matrix* matrix::getEigenvalueVector()
{
    //return E;
    return NULL;
}

//-----------------------------------------------------------------------------
// Name : getSize ()
//-----------------------------------------------------------------------------
int matrix::getSize()
{
    return N;
}
