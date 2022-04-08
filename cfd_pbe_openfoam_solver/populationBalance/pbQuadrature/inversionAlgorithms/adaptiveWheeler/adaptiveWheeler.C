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

This offering is not approved or endorsed by OpenCFD Limited (ESI Group)
and the OpenFOAM Foundation.
---------------------------------------------------------------------------- */

#include "adaptiveWheeler.H"
#include "eigenSolver.H"
#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace inversionAlgorithms
{
    defineTypeNameAndDebug(adaptiveWheeler, 0);
    addToRunTimeSelectionTable(inversionAlgorithm, adaptiveWheeler, populationBalance);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inversionAlgorithms::adaptiveWheeler::adaptiveWheeler
(
    const populationBalance& pb,
    PtrList<volScalarField>& nodes,
    PtrList<volScalarField>& weights
)
:
    inversionAlgorithm(pb, nodes, weights)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inversionAlgorithms::adaptiveWheeler::~adaptiveWheeler()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inversionAlgorithms::adaptiveWheeler::inversion
(
    const Field<scalar>& validCells
)
{
    // Number of moments
    const label nMoms = 2*numOfNodes_;

    const volScalarField& M0 = moments_[0];

    // Loop over the cells
    forAll(validCells, celli)
    {
        // Number of nodes
        int n = numOfNodes_;

        // If calculation is valid, find nodes and weights using PD algorithm
        if(validCells[celli]>0)
        {
            // Defining sigma matrix as a one-dimensional array
            // with initial zero elements
            List<scalar> sigma(n*(n + 1), Zero);
            
            word allMoments("");
            // Assigning the bottom row of the 2D sigma matrix
            for (int i=0; i<nMoms; i++)
            {
                sigma[i] = moments_[i][celli];

                allMoments += moments_[i].name() + " = " 
                            + Foam::name(sigma[i]) + " | ";
            }

            // Info<< allMoments << endl;

            // Defining array for diagonal elements of the Jacobi matrix
            // with initial zero elements
            List<scalar> a(n, Zero);

            // Assigning the first element of array "a"
            scalar atmp = sigma[1]/sigma[0];
            a[0] = atmp;

            // Assigning the second row of the sigma matrix
            for (int i=0; i<(nMoms - 2); i++)
            {
                sigma[i + nMoms] = sigma[i + 2] - atmp*sigma[i + 1];
            }

            // Defining array for off-diagonal elements of the Jacobi matrix
            // with initial zero elements
            List<scalar> b(n-1, Zero);

            // Assigning the first element of the array "b"
            scalar btmp = sigma[nMoms] / sigma[0];
            b[0] = btmp;

            // Defining zeta array for realizability check
            List<scalar> zeta(nMoms - 1, Zero);

            zeta[0] = atmp;

            scalar zetatmp = btmp / atmp;
            zeta[1] = zetatmp;

            if (zetatmp > 0)
            {
                // Assigning the second element of the array "a"
                atmp = sigma[nMoms + 1] / sigma[nMoms] - atmp;
                a[1] = atmp;

                zeta[2] = atmp - zetatmp;
            }

            // Loop over the next rows of the sigma matrix
            for (int j=2; j<n; j++)
            {
                if (zetatmp <= 0)
                {
                    break;
                }

                label nCol = nMoms - 2*j;

                // Auxiliary indexes for mapping the sigma matrix from 2D to 1D
                label fb = nMoms*j + j*(1 - j);
                label fb_1 = fb - nCol - 2;
                label fb_2 = fb_1 - nCol - 4;

                // Loop over the columns of each row and assign the elements
                for (int i=0; i<nCol; i++)
                {
                    sigma[i + fb] = 
                                    sigma[i + 2 + fb_1]
                                  - atmp * sigma[i + 1 + fb_1]
                                  - btmp * sigma[i + 2 + fb_2];
                }

                atmp = sigma[fb + 1]/sigma[fb]
                     - sigma[fb_1 + 1]/sigma[fb_1];
                btmp = sigma[fb] / sigma[fb_1];
                
                a[j] = atmp;
                b[j-1] = btmp;

                scalar zetatmp = btmp / zeta[2*j-2];
                zeta[2*j-1] = zetatmp;
                zeta[2*j] = atmp - zeta[2*j-1];
            }

            // Reduce number of nodes (n) if the moments are unrealizable
            n = nodeReduction(zeta, n);

            // Reduced "a" and "b" arrays (in case of unrealizable moment set)
            scalar aR[n];
            scalar bR[n-1];

            for (int i=0; i<n-1; i++)
            {
                aR[i] = a[i];
                bR[i] = -1.0 * Foam::sqrt(b[i]);
            }

            aR[n-1] = a[n-1];

            // Declaring eigen vector of the Jacobi matrix
            scalar eigVector[n][n];
            memset(eigVector, 0, n*n*sizeof eigVector[0][0]);
            
            scalar work[2*n - 2];
            int info;
            char choice='I';
            
            // Eigen vector calculation of the Jacobi matrix
            dsteqr_(choice, &n, aR, bR, &eigVector[0][0], &n, work, &info);
            
            // if (realizabileMoms)
            // {
                // Calculating nodes and weights using the eigen vectors
                for (int i=0; i<n; i++)
                {
                    if (aR[i] > 0)
                    {
                        nodes_[i][celli] = aR[i];
                        weights_[i][celli] = M0[celli]
                                        * Foam::pow(eigVector[i][0], 2);
                    }
                    else
                    {
                        WarningInFunction << "Inversion of moments failed for "
                            << "cell " << celli << endl
                            << "located at " << nodes_[i].mesh().C()[celli]
                            << endl << "with the following moments:" << endl
                            << allMoments << endl;
                    }
                    // Info<<celli<<endl;
                    // Info<< nodes_[i][celli]<< endl;
                    // Info<< weights_[i][celli]<< endl;
                }
            // }
            // else
            // {
                for (int i=n; i<numOfNodes_; i++)
                // for (int i=0; i<numOfNodes_; i++)
                {
                    nodes_[i][celli] = dSmall_[i];
                    weights_[i][celli] = 0.0;
                    // Info<<celli<<endl;
                    // Info<< nodes_[i][celli]<< endl;
                    // Info<< weights_[i][celli]<< endl;
                }
            // }
        }
        /* Assigning nodes and weights to arbitrary values, if the calculation
           is not valid*/ 
        else
        {
            for (int i=0; i<n; i++)
            {
                nodes_[i][celli] = dSmall_[i];
                weights_[i][celli] = 0.0;
            }
        }
    }
    forAll(nodes_, nodei)
    {
        nodes_[nodei].correctBoundaryConditions();
        weights_[nodei].correctBoundaryConditions();
    }
}


Foam::List<Foam::List<Foam::scalar>>
Foam::inversionAlgorithms::adaptiveWheeler::inversion
(
    const List<scalar>& moments
) const
{
    List<List<scalar>> nodesAndWeights(2);

    List<scalar>& nodes = nodesAndWeights[0];
    List<scalar>& weights = nodesAndWeights[1];

    const scalar M0 = moments[0];

    if (M0 > 1e-6 && moments[1] > 0.0)
    {
        // Number of moments
        const label nMoms = 2*numOfNodes_;

        // Number of nodes
        int n = numOfNodes_;

        // Defining sigma matrix as a one-dimensional array
        // with initial zero elements
        List<scalar> sigma(n*(n + 1), Zero);
        
        // Assigning the bottom row of the 2D sigma matrix
        for (int i=0; i<nMoms; i++)
        {
            sigma[i] = moments[i];
        }

        // Defining array for diagonal elements of the Jacobi matrix
        // with initial zero elements
        List<scalar> a(n, Zero);

        // Assigning the first element of array "a"
        scalar atmp = sigma[1]/sigma[0];
        a[0] = atmp;

        // Assigning the second row of the sigma matrix
        for (int i=0; i<(nMoms - 2); i++)
        {
            sigma[i + nMoms] = sigma[i + 2] - atmp*sigma[i + 1];
        }

        // Defining array for off-diagonal elements of the Jacobi matrix
        // with initial zero elements
        List<scalar> b(n-1, Zero);

        // Assigning the first element of the array "b"
        scalar btmp = sigma[nMoms] / sigma[0];
        b[0] = btmp;

        // Defining zeta array for realizability check
        List<scalar> zeta(nMoms - 1, Zero);

        zeta[0] = atmp;

        scalar zetatmp = btmp / atmp;
        zeta[1] = zetatmp;

        if (zetatmp > 0)
        {
            // Assigning the second element of the array "a"
            atmp = sigma[nMoms + 1] / sigma[nMoms] - atmp;
            a[1] = atmp;

            zeta[2] = atmp - zetatmp;
        }

        // Loop over the next rows of the sigma matrix
        for (int j=2; j<n; j++)
        {
            if (zetatmp <= 0)
            {
                break;
            }

            label nCol = nMoms - 2*j;

            // Auxiliary indexes for mapping the sigma matrix from 2D to 1D
            label fb = nMoms*j + j*(1 - j);
            label fb_1 = fb - nCol - 2;
            label fb_2 = fb_1 - nCol - 4;

            // Loop over the columns of each row and assign the elements
            for (int i=0; i<nCol; i++)
            {
                sigma[i + fb] = 
                                sigma[i + 2 + fb_1]
                                - atmp * sigma[i + 1 + fb_1]
                                - btmp * sigma[i + 2 + fb_2];
            }

            atmp = sigma[fb + 1]/sigma[fb]
                    - sigma[fb_1 + 1]/sigma[fb_1];
            btmp = sigma[fb] / sigma[fb_1];
            
            a[j] = atmp;
            b[j-1] = btmp;

            scalar zetatmp = btmp / zeta[2*j-2];
            zeta[2*j-1] = zetatmp;
            zeta[2*j] = atmp - zeta[2*j-1];
        }

        // Reduce number of nodes (n) if the moments are unrealizable
        n = nodeReduction(zeta, n);

        // Reduced "a" and "b" arrays (in case of unrealizable moment set)
        scalar aR[n];
        scalar bR[n-1];

        for (int i=0; i<n-1; i++)
        {
            aR[i] = a[i];
            bR[i] = -1.0 * Foam::sqrt(b[i]);
        }

        aR[n-1] = a[n-1];

        // Declaring eigen vector of the Jacobi matrix
        scalar eigVector[n][n];
        memset(eigVector, 0, n*n*sizeof eigVector[0][0]);
        
        scalar work[2*n - 2];
        int info;
        char choice='I';
        
        // Eigen vector calculation of the Jacobi matrix
        dsteqr_(choice, &n, aR, bR, &eigVector[0][0], &n, work, &info);

        nodes.setSize(n);
        weights.setSize(n);

        // Calculating nodes and weights using the eigen vectors
        for (int i=0; i<n; i++)
        {
            if (aR[i] > 0)
            {
                nodes[i] = aR[i];
                weights[i] = M0 * Foam::pow(eigVector[i][0], 2);
            }
            else
            {
                WarningInFunction << "Inversion of moments failed!" << endl;
            }
        }
    }
    else
    {
        nodes.setSize(numOfNodes_);
        weights.setSize(numOfNodes_);

        for (int i=0; i<numOfNodes_; i++)
        {
            nodes[i] = dSmall_[i];
            weights[i] = 0.0;
        }
    }

    return nodesAndWeights;
}


void Foam::inversionAlgorithms::adaptiveWheeler::inversionBoundary()
{
    // Number of moments
    const label nMoms = 2*numOfNodes_;

    const volScalarField::Boundary& bM0 = moments_[0].boundaryField();

    forAll(bM0, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(bM0[patchi]))
        {
            forAll(bM0[patchi], facei)
            {
                // Number of nodes
                int n = numOfNodes_;

                bool validMoments = true;

                // Defining sigma matrix as a one-dimensional array
                // with initial zero elements
                List<scalar> sigma(n*(n + 1), Zero);
                
                word allMoments("");
                // Assigning the bottom row of the 2D sigma matrix
                for (int i=0; i<nMoms; i++)
                {
                    const fvPatchScalarField& pMi
                        = moments_[i].boundaryField()[patchi];

                    if (pMi[facei] <= 0.0)
                    {
                        validMoments = false;
                    }

                    sigma[i] = pMi[facei];

                    allMoments += moments_[i].name() + " = " 
                                + Foam::name(pMi[facei]) + " | ";
                }

                // Info<< allMoments << endl;

                if (validMoments)
                {
                    // Defining array for diagonal elements of the Jacobi matrix
                    // with initial zero elements
                    List<scalar> a(n, Zero);

                    // Assigning the first element of array "a"
                    scalar atmp = sigma[1]/sigma[0];
                    a[0] = atmp;
                    
                    // Assigning the second row of the sigma matrix
                    for (int i=0; i<(nMoms - 2); i++)
                    {
                        sigma[i + nMoms] = sigma[i + 2] - atmp*sigma[i + 1];
                    }

                    // Defining array for off-diagonal elements of the Jacobi matrix
                    // with initial zero elements
                    List<scalar> b(n-1, Zero);

                    // Assigning the first element of the array "b"
                    scalar btmp = sigma[nMoms] / sigma[0];
                    b[0] = btmp;

                    // Defining zeta array for realizability check
                    List<scalar> zeta(nMoms - 1, Zero);

                    zeta[0] = atmp;

                    scalar zetatmp = btmp / atmp;
                    zeta[1] = zetatmp;

                    if (zetatmp > 0)
                    {
                        // Assigning the second element of the array "a"
                        atmp = sigma[nMoms + 1] / sigma[nMoms] - atmp;
                        a[1] = atmp;

                        zeta[2] = atmp - zetatmp;
                    }

                    // Loop over the next rows of the sigma matrix
                    for (int j=2; j<n; j++)
                    {
                        if (zetatmp <= 0)
                        {
                            break;
                        }

                        label nCol = nMoms - 2*j;

                        // Auxiliary indexes for mapping the sigma matrix from
                        // 2D to 1D
                        label fb = nMoms*j + j*(1 - j);
                        label fb_1 = fb - nCol - 2;
                        label fb_2 = fb_1 - nCol - 4;

                        // Loop over the columns of each row
                        // and assign the elements
                        for (int i=0; i<nCol; i++)
                        {
                            sigma[i + fb] = 
                                            sigma[i + 2 + fb_1]
                                        - atmp * sigma[i + 1 + fb_1]
                                        - btmp * sigma[i + 2 + fb_2];
                        }
                        
                        atmp = sigma[fb + 1]/sigma[fb]
                            - sigma[fb_1 + 1]/sigma[fb_1];
                        btmp = sigma[fb] / sigma[fb_1];
                        
                        a[j] = atmp;
                        b[j-1] = btmp;

                        scalar zetatmp = btmp / zeta[2*j-2];
                        zeta[2*j-1] = zetatmp;
                        zeta[2*j] = atmp - zeta[2*j-1];
                    }

                    // Reduce number of nodes (n) if the moments are unrealizable
                    n = nodeReduction(zeta, n);

                    // Reduced "a" and "b" arrays (in case of unrealizable moment set)
                    scalar aR[n];
                    scalar bR[n-1];

                    for (int i=0; i<n-1; i++)
                    {
                        aR[i] = a[i];
                        bR[i] = -1.0 * Foam::sqrt(b[i]);
                    }

                    aR[n-1] = a[n-1];

                    // Declaring eigen vector of the Jacobi matrix
                    scalar eigVector[n][n];
                    memset(eigVector, 0, n*n*sizeof eigVector[0][0]);
                    
                    scalar work[2*n - 2];
                    int info;
                    char choice='I';
                    
                    // Eigen vector calculation of the Jacobi matrix
                    dsteqr_(choice, &n, aR, bR, &eigVector[0][0], &n, work, &info);
                    
                    // Calculating nodes and weights using the eigen vectors
                    for (int i=0; i<n; i++)
                    {
                        if (aR[i] > 0)
                        {
                            // Better to find a way for referencing out of loop!
                            nodes_[i].boundaryFieldRef()[patchi][facei] = aR[i];
                            weights_[i].boundaryFieldRef()[patchi][facei] =
                                bM0[patchi][facei] * Foam::pow(eigVector[i][0], 2);
                            // Info<< nodes_[i].boundaryField()[patchi][facei]<< endl;
                            // Info<< weights_[i].boundaryField()[patchi][facei]<< endl;
                        }
                        else
                        {
                            FatalErrorInFunction << "negative node " << aR[i]
                                << " calculated for face " << facei
                                << " on the patch " << patchi << endl
                                << "located at "
                                << nodes_[0].mesh().Cf().boundaryField()[patchi][facei]
                                << endl << "with following moments:" << endl
                                << allMoments << endl << exit(FatalError);
                        }
                    }

                    for (int i=n; i<numOfNodes_; i++)
                    {
                        nodes_[i].boundaryFieldRef()[patchi][facei] = dSmall_[i];
                        weights_[i].boundaryFieldRef()[patchi][facei] = 0.0;
                        // Info<<celli<<endl;
                        // Info<< nodes_[i][celli]<< endl;
                        // Info<< weights_[i][celli]<< endl;
                    }
                }
                else
                {
                    for (int i=0; i<n; i++)
                    {
                        nodes_[i].boundaryFieldRef()[patchi][facei] = dSmall_[i];
                        weights_[i].boundaryFieldRef()[patchi][facei] = 0.0;
                    }
                }
                
            }
        }
    }
}


// ************************************************************************* //
