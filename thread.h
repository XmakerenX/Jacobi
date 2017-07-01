#ifndef  _THREAD_H
#define  _THREAD_H

#include "matrix.h"

typedef void* (*THREADFUNC )(void* args );

class Thread
{

public:
    Thread (matrix& S, matrix& E, int* ind, int& k, int& l, double& c, double& s, THREADFUNC func);
    ~Thread();
    bool start();
    void wait();
    void* rotateVector(void* arg);

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


#endif  //_THREAD_H
