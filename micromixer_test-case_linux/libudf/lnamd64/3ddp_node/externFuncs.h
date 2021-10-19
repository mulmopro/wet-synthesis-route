/* ----------------------------------------------------------------------------
Copyright © 2021 Politecnico di Torino

This file is part of WetSynthRoute.

WetSynthRoute is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

WetSynthRoute is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with WetSynthRoute.  If not, see <https://www.gnu.org/licenses/>.
---------------------------------------------------------------------------- */


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
