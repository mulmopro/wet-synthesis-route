#ifndef EIGEN_SOLVER_H
#define EIGEN_SOLVER_H

extern "C"
{
    void dsteqr_
    (
        char &, int *, double *, double *, double *, int *, double *, int *
    );
}

#endif
