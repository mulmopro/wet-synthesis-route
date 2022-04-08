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

#ifndef inversionAlgorithm_H
#define inversionAlgorithm_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalance.H"
#include "quadratureMethod.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class inversionAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class inversionAlgorithm
{
protected:

    // Protected data

        //- Reference to the moments
        const PtrList<volScalarField>& moments_;

        //- Reference to the diameters for invalid cells
        const List<scalar>& dSmall_;

        //- Minimum value for alpha used to check moment no. 3
        const dimensionedScalar& alphaMin_;

        //- Shape factor
        const scalar kv_;

        //- Reference to the number of nodes
        const uint8_t numOfNodes_;

        //- Reference to nodes of the NDF
        PtrList<volScalarField>& nodes_;

        //- Reference to weights of the NDF
        PtrList<volScalarField>& weights_;

        //- Node reduction based on the zeta values
        int nodeReduction(const List<scalar>&, int) const;

public:

    //- Runtime type information
    TypeName("inversionAlgorithm");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            inversionAlgorithm,
            populationBalance,
            (
                const populationBalance& pb,
                PtrList<volScalarField>& nodes,
                PtrList<volScalarField>& weights
            ),
            (pb, nodes, weights)
        );

    // Constructors

        inversionAlgorithm
        (
            const populationBalance& pb,
            PtrList<volScalarField>& nodes,
            PtrList<volScalarField>& weights
        );


    //- Destructor
    virtual ~inversionAlgorithm();


    // Selectors

        static autoPtr<inversionAlgorithm> New
        (
            const populationBalance& pb,
            PtrList<volScalarField>& nodes,
            PtrList<volScalarField>& weights
        );


    // Member Functions

        //- Calculates nodes and weights of the reconstructed NDF
        virtual void inversion(const Field<scalar>&) = 0;

        //- Calculates nodes and weights for a single cell
        virtual List<List<scalar>> inversion(const List<scalar>&) const = 0;

        //- Calculates nodes and weights of the reconstructed NDF on the
        //- fixed value boundaries
        virtual void inversionBoundary() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
