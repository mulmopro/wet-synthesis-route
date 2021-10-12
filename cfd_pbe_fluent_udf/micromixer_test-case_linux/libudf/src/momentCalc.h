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


int nodeReduction(const double* zeta, int n)
{
    int i;
    int rN = 1;

    if (zeta[0] <= 0)
    {
        /* printf("\nWarning:\nNumber of nodes is reduced to %u \n", rN); */
        return rN;
    }
    
    for (i=2; i<2*n-1; i+=2)
    {
        if (zeta[i] <= 0 || zeta[i-1] <= 0)
        {
            /* printf("\nWarning:\nNumber of nodes is reduced to %u \n", rN); */
            break;
        }
        rN++;
    }

    return rN;
}


void adaptiveWheeler(double* nodes, double* weights, const double* moments)
{
    /* Number of moments */
    const int nMoms = 2*N_NODES;

    int i, j;

    /* Defining sigma matrix as a one-dimensional array */
    double sigma[N_NODES*(N_NODES + 1)];

    /* Assigning the bottom row of the 2D sigma matrix */
    for (i=0; i<nMoms; i++)
    {
        sigma[i] = moments[i];
    }

    /* Defining array for diagonal elements of the Jacobi matrix */
    double a[N_NODES];

    /* Assigning the first element of array "a" */
    double atmp = sigma[1]/sigma[0];
    a[0] = atmp;

    /* Assigning the second row of the sigma matrix */
    for (i=0; i<(nMoms - 2); i++)
    {
        sigma[i + nMoms] = sigma[i + 2] - atmp*sigma[i + 1];
    }

    /* Defining array for off-diagonal elements of the Jacobi matrix */
    double b[N_NODES-1];

    /* Assigning the first element of the array "b" */
    double btmp = sigma[nMoms] / sigma[0];
    b[0] = btmp;

    /* Defining zeta array for realizability check */
    double zeta[2*N_NODES - 1];

    zeta[0] = atmp;

    double zetatmp = btmp / atmp;
    zeta[1] = zetatmp;

    if (zetatmp > 0)
    {
        /* Assigning the second element of the array "a" */
        atmp = sigma[nMoms + 1] / sigma[nMoms] - atmp;
        a[1] = atmp;

        zeta[2] = atmp - zetatmp;
    }

    /* Loop over the next rows of the sigma matrix */
    for (j=2; j<N_NODES; j++)
    {
        if (zetatmp <= 0)
        {
            break;
        }

        int nCol = nMoms - 2*j;

        /* Auxiliary indexes for mapping the sigma matrix from 2D to 1D */
        int fb = nMoms*j + j*(1 - j);
        int fb_1 = fb - nCol - 2;
        int fb_2 = fb_1 - nCol - 4;

        /* Loop over the columns of each row and assign the elements */
        for (i=0; i<nCol; i++)
        {
            sigma[i + fb] = sigma[i + 2 + fb_1] - atmp * sigma[i + 1 + fb_1]
                - btmp * sigma[i + 2 + fb_2];
        }

        atmp = sigma[fb + 1]/sigma[fb] - sigma[fb_1 + 1]/sigma[fb_1];
        btmp = sigma[fb] / sigma[fb_1];
        
        a[j] = atmp;
        b[j-1] = btmp;

        double zetatmp = btmp / zeta[2*j-2];
        zeta[2*j-1] = zetatmp;
        zeta[2*j] = atmp - zeta[2*j-1];
    }

    /* Reduce number of nodes (n) if the moments are non-realizable */
    int rN;
    rN = nodeReduction(zeta, N_NODES);

    /* Reduced "a" and "b" arrays (in case of non-realizable moment set) */
    double *aR = (double *)malloc(rN * sizeof(double));
    double *bR = (double *)malloc((rN - 1) * sizeof(double));
    /* Declaring eigen vector of the Jacobi matrix */
    double *eigVector = (double *)malloc(rN * rN * sizeof(double));
    double *work = (double *)malloc((2*rN - 2) * sizeof(double));

    if (aR == NULL || bR == NULL || eigVector == NULL || work == NULL)
    {
        printf("\nError in PD algorithm:\n");
        printf("malloc of arrays with type \"double\" failed!\n");
        printf("Please execute the UDF again.\n");
        /* exit(1); */
    }

    for (i=0; i<rN-1; i++)
    {
        aR[i] = a[i];
        bR[i] = -1.0 * pow(b[i], 0.5);
    }

    aR[rN-1] = a[rN-1];

    memset(eigVector, 0, rN*rN*sizeof(eigVector[0]));
    
    int info;
    char choice = 'I';
    
    /* Eigen vector calculation of the Jacobi matrix */
    dsteqr(&choice, &rN, aR, bR, eigVector, &rN, work, &info);
    
    /* Calculating nodes and weights using the eigen vectors */
    for (i=0; i<rN; i++)
    {
        if (aR[i] > 0)
        {
            nodes[i] = aR[i];
            weights[i] = moments[0] * pow(eigVector[i*rN], 2);
        }
        else if (a[i] < 0)
        {
            printf("\nError:\nNegative abscissa is detected\n");

            nodes[i] = 1.0e-32;
            weights[i] = moments[0] * pow(eigVector[i*rN], 2);
        }
        else
        {
            printf("\nWarning:\nZero abscissa is detected\n");

            nodes[i] = 1.0e-32;
            weights[i] = moments[0] * pow(eigVector[i*rN], 2);
        }
    }

    free(aR);
    free(bR);
    free(eigVector);
    free(work);

    for (i=rN; i<N_NODES; i++)
    {
        nodes[i] = SMALL_SIZE;
        weights[i] = 0.0;
    }
}


void PD(double* nodes, double* weights, const double* moments)
{
    int i, j;

    /* Defining an auxiliary index */
    int fIndex = 2*N_NODES + 1;

    /* If calculation is valid, find nodes and weights using PD algorithm */
    /* if(validCell)
    { */

    /* Defining P matrix as a one-dimensional array */
    double P[(2*N_NODES + 1)*(N_NODES + 1)];

    /* Assigning the elements corresponding to the first column of P */
    P[0] = 1.0; /* Assigning the first element to 1.0 */
    for (i=1; i<fIndex; i++)
    {
        P[i] = 0.0;
    }
    
    /* Assigning the elements corresponding to the second column of P */
    for (i=0; i<(fIndex - 1); i++)
    {
        P[i + fIndex] = pow(-1.0, i) * moments[i];
    }

    /* Defining zeta array */
    double zeta[2*N_NODES];

    zeta[0] = 0.0;

    /*  Assigning other elements of the P and
        calculating the elements of the zeta at the same time */
    for (j=2; j<(fIndex); j++)
    {
        /* Auxiliary indexes for mapping the 2D P matrix to 1D P array */
        int fb = j*(2*fIndex + 1 - j) / 2;
        int fb_1 = fb - fIndex + j - 1;
        int fb_2 = fb_1 - fIndex + j - 2;

        for (i=0; i<(fIndex - j); i++)
        {
            P[i + fb] = 
                        (
                            P[fb_1] * P[i + 1 + fb_2]
                          - P[fb_2] * P[i + 1 + fb_1]
                        );
        }
        
        double product = P[fb_1] * P[fb_2];

        if (product > 0)
        {
            zeta[j - 1] = P[fb]/product;
        }
        /* else
        {
            zeta[1] = 0.0;
        } */
    }

    /* Declaring diagonal and off-diagonal elements of Jacobi matrix */
    double a[N_NODES];
    double b[N_NODES-1];

    /* Declaring eigen vector of the Jacobi matrix */
    double eigVector[N_NODES*N_NODES];
    memset(eigVector, 0, N_NODES*N_NODES*sizeof(eigVector[0]));
    
    double work[2*N_NODES - 2];
    int info;
    char choice = 'I';

    /* Calculating the elements of the Jacobi matrix */
    for (i=0; i<(N_NODES - 1); i++)
    {
        a[i] = zeta[2*i + 1] + zeta[2*i];

        double tempProd = zeta[2*i + 2] * zeta[2*i + 1];

        if (tempProd < 0)
        {
            printf("\nError:\nMoments are non-realizable\n");
        }

        b[i] = -1.0 * pow(tempProd, 0.5);
    }

    /*  Calculating last element of the diagonal out of the loop
        since it has one element more than off-diagonal */
    a[N_NODES-1] = zeta[fIndex - 2] + zeta[fIndex - 3];
    
    int n = N_NODES;
    /* Eigen vector calculation of the Jacobi matrix */
    dsteqr(&choice, &n, a, b, eigVector, &n, work, &info);
    
    /* Calculating nodes and weights using the eigen vectors */
    for (i=0; i<N_NODES; i++)
    {
        if (a[i] > 0)
        {
            nodes[i] = a[i];
            weights[i] = moments[0] * pow(eigVector[i*N_NODES], 2);
        }
        else if (a[i] < 0)
        {
            printf("\nError:\nNegative abscissa is detected\n");

            nodes[i] = 1.0e-32;
            weights[i] = moments[0] * pow(eigVector[i*N_NODES], 2);
        }
        else
        {
            printf("\nWarning:\nZero abscissa is detected\n");

            nodes[i] = 1.0e-32;
            weights[i] = moments[0] * pow(eigVector[i*N_NODES], 2);
        }
    }
}
