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

#include "PD.H"
#include "eigenSolver.H"
#include "fixedValueFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace inversionAlgorithms
{
    defineTypeNameAndDebug(PD, 0);
    addToRunTimeSelectionTable(inversionAlgorithm, PD, populationBalance);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inversionAlgorithms::PD::PD
(
    const populationBalance& pb,
    PtrList<volScalarField>& nodes,
    PtrList<volScalarField>& weights
)
:
    inversionAlgorithm(pb, nodes, weights)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inversionAlgorithms::PD::~PD()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inversionAlgorithms::PD::inversion
(
    const Field<scalar>& validCells
)
{
    int n = numOfNodes_;

    // Defining an auxiliary index
    label fIndex = 2*n + 1;

    // Loop over the cells
    forAll(validCells, celli)
    {
        // If calculation is valid, find nodes and weights using PD algorithm
        if(validCells[celli]>0)
        {
            // Defining P matrix as a one-dimensional array
            // with initial zero elements
            List<scalar> P(fIndex*(fIndex + 1)/2, Zero);

            // Assigning the first element to 1.0
            P[0] = 1.0;
            
            word allMoments("");
            // Assigning the elements corresponding to the second column of P
            for (int i=0; i<(fIndex - 1); i++)
            {
                P[i + fIndex] = Foam::pow(-1.0, i) * moments_[i][celli];

                allMoments += moments_[i].name() + " = " 
                            + Foam::name(moments_[i][celli]) + " | ";
            }

            // Info<< allMoments << endl;

            // Defining zeta array
            List<scalar> zeta(fIndex - 1, Zero);

            // Assigning other elements of the P and
            // calculating the elements of the zeta at the same time
            for (int j=2; j<(fIndex); j++)
            {
                // Auxiliary indexes for mapping the 2D P matrix to 1D P array
                label fb = j*(2*fIndex + 1 - j) / 2;
                label fb_1 = fb - fIndex + j - 1;
                label fb_2 = fb_1 - fIndex + j - 2;
        
                for (int i=0; i<(fIndex - j); i++)
                {
                    P[i + fb] = 
                              (
                                  P[fb_1] * P[i + 1 + fb_2]
                                - P[fb_2] * P[i + 1 + fb_1]
                              );
                }
                
                scalar product = P[fb_1] * P[fb_2];
 
                if (product > 0)
                {
                    zeta[j - 1] = P[fb]/product;
                }
                // else
                // {
                //     zeta[1] = 0.0;
                // }
            }

            // Declaring diagonal and off-diagonal elements of Jacobi matrix
            scalar a[n];
            scalar b[n-1];

            // Declaring eigen vector of the Jacobi matrix
            scalar eigVector[n][n];
            memset(eigVector, 0, n*n*sizeof eigVector[0][0]);
            
            scalar work[fIndex - 3];
            int info;
            char choice='I';

            // Calculating the elements of the Jacobi matrix
            for (int i=0; i<(n - 1); i++)
            {
                a[i] = zeta[2*i + 1] + zeta[2*i];
                b[i] = -1.0 * Foam::sqrt(zeta[2*i + 2] * zeta[2*i + 1]);
            }

            // Calculating last element of the diagonal out of the loop
            // since it has one element more than off-diagonal 
            a[n-1] = zeta[fIndex - 2] + zeta[fIndex - 3];
            
            // Eigen vector calculation of the Jacobi matrix
            dsteqr_(choice, &n, a, b, &eigVector[0][0], &n, work, &info);
            
            // Calculating nodes and weights using the eigen vectors
            for (int i=0; i<n; i++)
            {
                if (a[i] < 0)
                {
                    FatalErrorInFunction << "negative node " << a[i]
                        << " calculated for cell " << celli << endl
                        << "located at " << nodes_[i].mesh().C()[celli]
                        << endl << "with following moments:" << endl
                        << allMoments << endl << exit(FatalError);
                }
                nodes_[i][celli] = a[i];
                weights_[i][celli] = moments_[0][celli] *
                                        Foam::pow(eigVector[i][0], 2);
                // Info<<celli<<endl;
                // Info<< nodes_[i][celli]<< endl;
                // Info<< weights_[i][celli]<< endl;
            }
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


void Foam::inversionAlgorithms::PD::inversionBoundary()
{
    int n = numOfNodes_;

    // Defining an auxiliary index
    label fIndex = 2*n + 1;

    const volScalarField::Boundary& bM0 = moments_[0].boundaryField();

    forAll(bM0, patchi)
    {
        if (isA<fixedValueFvPatchScalarField>(bM0[patchi]))
        {
            forAll(bM0[patchi], facei)
            {
                // Defining P matrix as a one-dimensional array
                // with initial zero elements
                List<scalar> P(fIndex*(fIndex + 1)/2, Zero);

                // Assigning the first element to 1.0
                P[0] = 1.0;
                
                word allMoments("");
                // Assigning the elements corresponding to the second column of P
                for (int i=0; i<(fIndex - 1); i++)
                {
                    const fvPatchScalarField& pMi
                        = moments_[i].boundaryField()[patchi];
                    
                    P[i + fIndex] = Foam::pow(-1.0, i) * pMi[facei];

                    allMoments += moments_[i].name() + " = " 
                                + Foam::name(pMi[facei])
                                + " | ";
                }

                // Info<< allMoments << endl;

                // Defining zeta array
                List<scalar> zeta(fIndex - 1, Zero);

                // Assigning other elements of the P and
                // calculating the elements of the zeta in the same loop
                for (int j=2; j<(fIndex); j++)
                {
                    // Auxiliary indexes for mapping the 2D P matrix to 1D P array
                    label fb = j*(2*fIndex + 1 - j) / 2;
                    label fb_1 = fb - fIndex + j - 1;
                    label fb_2 = fb_1 - fIndex + j - 2;
            
                    for (int i=0; i<(fIndex - j); i++)
                    {
                        P[i + fb] = 
                                (
                                    P[fb_1] * P[i + 1 + fb_2]
                                    - P[fb_2] * P[i + 1 + fb_1]
                                );
                    }
                    
                    scalar product = P[fb_1] * P[fb_2];
    
                    if (product > 0)
                    {
                        zeta[j - 1] = P[fb]/product;
                    }
                    // else
                    // {
                    //     zeta[1] = 0.0;
                    // }
                }

                // Declaring diagonal and off-diagonal elements of Jacobi matrix
                scalar a[n];
                scalar b[n-1];

                // Declaring eigen vector of the Jacobi matrix
                scalar eigVector[n][n];
                memset(eigVector, 0, n*n*sizeof eigVector[0][0]);
                
                scalar work[fIndex - 3];
                int info;
                char choice='I';

                // Calculating the elements of the Jacobi matrix
                for (int i=0; i<(n - 1); i++)
                {
                    a[i] = zeta[2*i + 1] + zeta[2*i];
                    b[i] = -1.0 * Foam::sqrt(zeta[2*i + 2] * zeta[2*i + 1]);
                }

                // Calculating last element of the diagonal out of the loop
                // since it has one element more than off-diagonal 
                a[n-1] = zeta[fIndex - 2] + zeta[fIndex - 3];

                // Eigen vector calculation of the Jacobi matrix
                dsteqr_(choice, &n, a, b, &eigVector[0][0], &n, work, &info);
                
                // Calculating nodes and weights using the eigen vectors
                for (int i=0; i<n; i++)
                {
                    if (a[i] < 0)
                    {
                        FatalErrorInFunction << "negative node " << a[i]
                            << " calculated for face " << facei << endl
                            << "located at "
                            << nodes_[0].mesh().Cf().boundaryField()[patchi][facei]
                            << endl << "with following moments:" << endl
                            << allMoments << endl << exit(FatalError);
                    }
                    // Better to find a way for referencing out of loop!
                    nodes_[i].boundaryFieldRef()[patchi][facei] = a[i];
                    weights_[i].boundaryFieldRef()[patchi][facei] =
                        bM0[patchi][facei] * Foam::pow(eigVector[i][0], 2);
                    // Info<< nodes_[i].boundaryField()[patchi][facei]<< endl;
                    // Info<< weights_[i].boundaryField()[patchi][facei]<< endl;
                }
            }
        }
    }
}



// ************************************************************************* //
