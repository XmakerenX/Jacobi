#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <sys/time.h>
#include <pthread.h>
#include <unistd.h>
#include <iomanip>
#include <mkl.h>

#include "main.h"
#include "matrix.h"
#include "thread.h"

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
// Name : main ()
//-----------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    bool mklMode = false;

    if (argc >= 4 || argc == 1)
    {
        std::cout << "wrong number of parameters was given\n";
        std::cout << "launch as:\n   Jacboi <file_path>\n";
        std::cout << "or Jacobi <file_path> (mkl for the mkl implementation)\n";
        return 1;
    }

    if (argc == 3)
        if (std::string(argv[2]) == "mkl")
            mklMode = true;

    std::cout << "the big jacobi test!" << std::endl;
    std::cout.precision(16);

    double* e;
    matrix mat(argv[1]);

    std::cout <<"Finished loading matrix\n";

    double t = getTime();

    if (!mklMode)
        e = mat.jacobi();
    else
        e = mat.jacobiMkl();

    std::cout << "Time: " << getTime () - t << "\n\n";

    // hackish sorting of the eigenValues returned
    bool swapped = false;
    do
    {
        swapped = false;
        for (int i = 1; i < mat.getSize(); i++)
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

    for (int i = 0; i < mat.getSize(); i++)
    {
        std::cout << "e" << i+1 << " = " << e[i] << "\n\n";

//        for (int j = 0; j < N; j++)
//            std::cout << E[j][i] <<" ";
    }

    std::cout << "Time: " << getTime () - t << "\n\n";

    return 0;
}
