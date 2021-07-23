#ifndef EIGEN_VEC_H
#define EIGEN_VEC_H

extern void dsteqr
(
    char* choice, int* n1, double* a, double* b, double* eigVector, int* n2,
    double* work, int* info
);

#endif

#ifndef LIN_SOLVER_H
#define LIN_SOLVER_H
extern void dgesv
(
    int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb,
    int* info
);

#endif
