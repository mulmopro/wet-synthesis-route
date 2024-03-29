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

#ifndef noPBM_H
#define noPBM_H

#include "populationBalance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace populationBalances
{

/*---------------------------------------------------------------------------*\
                      Class noPBM Declaration
\*---------------------------------------------------------------------------*/

class noPBM
:
    public populationBalance
{
    // Private data

        //- Minimum value for alpha used to check moment no. 3
        dimensionedScalar alphaMin_;

        //- Small diameters for the invalid cells
        List<scalar> dSmall_;

        //- Number of nodes
        const uint8_t numOfNodes_;

        //- Number of moments
        const uint8_t numOfMoments_;

        //- Number of coordinates
        const uint8_t numOfCoord_;

        // Precipitation rate
        DimensionedField<scalar, volMesh> precRate_;

        //- List of the moments
        PtrList<volScalarField> moments_;

        //- List of the cation concentration ratios
        PtrList<DimensionedField<scalar, volMesh>> cationConcRatios_;

public:

    //- Runtime type information
    TypeName("none");

    // Constructors

        //- Construct from fvMesh
        noPBM
        (
            const fvMesh& mesh,
            const incompressible::momentumTransportModel& turbulence,
            const solutionNMC& solution
        );


    //- Destructor
        virtual ~noPBM();


    // Member Functions

        bool highOrderScheme() const
        {
            return false;
        };

        // Access

            //- Return the growth model
            const growthModel& growth() const;

            //- Return the nucleation rate model
            const nucleationRateModel& nucleationRate() const;

            //- Return the nucleate size model
            const nucleateSizeModel& nucleateSize() const;

            //- Return the list of aggregation models
            const aggregationList& aggregation() const;

            //- Return the list of breakage models
            const breakageList& breakage() const;

            //- Return the quadrature model
            const quadratureMethod& quadrature() const;

            //- Return the phi field
            const surfaceScalarField& phi() const;

            //- Return the moments
            const PtrList<volScalarField>& moments() const;
            
            //- Return non-const access to the moments
            PtrList<volScalarField>& moments();

            //- Calculate the precipitation rate
            const DimensionedField<scalar, volMesh>& precRate() const;

            //- Return the supersaturation
            const DimensionedField<scalar, volMesh>& superSat() const;

            //- Return the cation concentration ratios
            const PtrList<DimensionedField<scalar, volMesh>>&
                cationConcRatios() const;

            //- Tagging cells for which source calculation and updating
            //  diameter are valid
            const DimensionedField<scalar, volMesh>& sourceCorrCoeff() const;

            //- Return the number of moments
            uint8_t numOfMoments() const;

            //- Return the number of nodes
            uint8_t numOfNodes() const;

            //- Return the number of coordinates
            uint8_t numOfCoord() const;

            //- Return Minimum value for alpha used to check moment no. 3
            const dimensionedScalar& alphaMin() const;

            //- Return the small diameters for the invalid cells
            const List<scalar>& dSmall() const;

            //- Return the crystal density
            scalar crystalRho() const;

            //- Return the crystal molecular weight
            scalar crystalMW() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalances
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
