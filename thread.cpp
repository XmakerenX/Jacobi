#include "thread.h"

//-----------------------------------------------------------------------------
// Constructor
//-----------------------------------------------------------------------------
Thread::Thread (matrix& S, matrix& E, int* ind, int& k, int& l, double& c, double& s, THREADFUNC func)
:S_(S), E_(E),k_(k), l_(l), c_(c), s_(s)
{
    threadFunc_ = func;

    ind_ = ind;

    threadReady = false;
    threadDone = false;

    startIndex = 0;
    endIndex = 0;
}

//-----------------------------------------------------------------------------
// Destructor
//-----------------------------------------------------------------------------
Thread::~Thread()
{
    pthread_cancel(thread_);
}

//-----------------------------------------------------------------------------
// start()
//-----------------------------------------------------------------------------
bool Thread::start()
{
    return (pthread_create(&thread_, NULL, threadFunc_, this) == 0);
}

//-----------------------------------------------------------------------------
// wait()
//-----------------------------------------------------------------------------
void Thread::wait()
{
    pthread_join(thread_, NULL);
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
            for (int i = 0; i  < T->S_.getSize(); i++)
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
