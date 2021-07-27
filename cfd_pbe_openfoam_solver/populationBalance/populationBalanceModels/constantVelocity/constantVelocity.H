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

#ifndef constantVelocity_H
#define constantVelocity_H

#include "populationBalance.H"
#include "multivariateSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionNMC;
class quadratureMethod;
class growthModel;
class nucleationRateModel;
class nucleateSizeModel;
class aggregationList;
class breakageList;

namespace populationBalances
{

/*---------------------------------------------------------------------------*\
                      Class constantVelocity Declaration
\*---------------------------------------------------------------------------*/

class constantVelocity
:
    public populationBalance
{
    // Private data

        //- Reference to the phi
        const surfaceScalarField& phi_;
        
        //- Reference to the supersaturation field
        const DimensionedField<scalar, volMesh>& superSat_;

        //- Reference to the supersaturation field
        const PtrList<DimensionedField<scalar, volMesh>>& cationConcRatios_;

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

        //- Crystal density
        dimensionedScalar crystalRho_;

        //- Crystal molecular weight
        dimensionedScalar crystalMW_;

        // Precipitation rate
        DimensionedField<scalar, volMesh> precRate_;

        // Used to tag cells for which source calculation and updating diameter
        // is valid
        DimensionedField<scalar, volMesh> sourceCorrCoeff_;

        // Needed by the convection scheme in the moment tansport equations
        multivariateSurfaceInterpolationScheme<scalar>::fieldTable mFields_;

        //- List of the moments
        PtrList<volScalarField> moments_;

        //- Growth model
        autoPtr<growthModel> growth_;

        //- Nucleation rate model
        autoPtr<nucleationRateModel> nucleationRate_;

        //- Nucleation size model
        autoPtr<nucleateSizeModel> nucleateSize_;

        //- List of aggregation models
        autoPtr<aggregationList> aggregationList_;
        
        //- List of breakage models
        autoPtr<breakageList> breakageList_;

        //- Quadrature method
        autoPtr<quadratureMethod> quadrature_;

public:

    //- Runtime type information
    TypeName("constantVelocity");

    // Constructors

        //- Construct from fvMesh
        constantVelocity
        (
            const fvMesh& mesh,
            const incompressible::momentumTransportModel& turbulence,
            const solutionNMC& solution
        );


    //- Destructor
        virtual ~constantVelocity();


    // Member Functions

        //- Solve for the PBE
        virtual void transport_moments();

        //- Solve for the PBE
        virtual void correct();

        //- Read pbProperties dictionary
        virtual bool read();

        bool highOrderScheme() const
        {
            return false;
        };

        // Access

            //- Return the growth rate model
            const growthModel& growth() const
            {
                return growth_();
            }

            //- Return the nucleation rate model
            const nucleationRateModel& nucleationRate() const
            {
                return nucleationRate_();
            }

            //- Return the nucleate size model
            const nucleateSizeModel& nucleateSize() const
            {
                return nucleateSize_();
            }

            //- Return the phi field
            const surfaceScalarField& phi() const
            {
                return phi_;
            }

            //- Return the moments
            const PtrList<volScalarField>& moments() const
            {
                return moments_;
            }
            
            //- Return non-const access to the moments
            PtrList<volScalarField>& moments()
            {
                return moments_;
            }

            //- Return the precipitation rate
            const DimensionedField<scalar, volMesh>& precRate() const
            {
                return precRate_;
            };

            //- Return the supersaturation
            const DimensionedField<scalar, volMesh>& superSat() const
            {
                return superSat_;
            }

            //- Return the cation concentration ratios
            const PtrList<DimensionedField<scalar, volMesh>>&
                cationConcRatios() const
            {
                return cationConcRatios_;
            }

            //- Tagging cells for which source calculation and updating diameter
            //  are valid
            const DimensionedField<scalar, volMesh>& sourceCorrCoeff() const
            {
                return sourceCorrCoeff_;
            }

            //- Return the number of moments
            const uint8_t& numOfMoments() const
            {
                return numOfMoments_;
            }

            //- Return the number of nodes
            const uint8_t& numOfNodes() const
            {
                return numOfNodes_;
            }

            //- Return the number of coordinates
            const uint8_t& numOfCoord() const
            {
                return numOfCoord_;
            }

            //- Return the small diameters for the invalid cells
            const List<scalar>& dSmall() const
            {
                return dSmall_;
            }

            //- Return the crystal density
            scalar crystalRho() const
            {
                return crystalRho_.value();
            }

            //- Return the crystal molecular weight
            scalar crystalMW() const
            {
                return crystalMW_.value();
            }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace populationBalances
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
