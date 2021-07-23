#ifndef LIN_SOLVER_H
#define LIN_SOLVER_H

extern "C"
{
    void dgesv_
    (
        int *, int *, double *, int *, int *, double *, int *, int *
    );
}

#endif
