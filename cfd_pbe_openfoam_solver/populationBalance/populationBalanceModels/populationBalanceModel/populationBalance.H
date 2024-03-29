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

#ifndef populationBalance_H
#define populationBalance_H

#include "IOdictionary.H"
#include "volFields.H"
#include "runTimeSelectionTables.H"

#include "kinematicMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionNMC;
class populationBalance;
class growthModel;
class nucleationRateModel;
class nucleateSizeModel;
class aggregationList;
class breakageList;
class quadratureMethod;

Ostream& operator<<(Ostream&, const populationBalance&);

/*---------------------------------------------------------------------------*\
                      Class populationBalance Declaration
\*---------------------------------------------------------------------------*/

class populationBalance
:
    public IOdictionary
{
    // Private data

protected:

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the turbulence
        const incompressible::momentumTransportModel& turbulence_;

public:
    //- Runtime type information
    TypeName("populationBalance");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            populationBalance,
            fvMesh,
            (
                const fvMesh& mesh,
                const incompressible::momentumTransportModel& turbulence,
                const solutionNMC& solution
            ),
            (mesh, turbulence, solution)
        );


    // Constructors

        //- Construct from fvMesh
        populationBalance
        (
            const fvMesh& mesh,
            const incompressible::momentumTransportModel& turbulence,
            const solutionNMC& solution
        );


    //- Destructor
        virtual ~populationBalance();


    // Selectors

    static autoPtr<populationBalance> New
    (
        const fvMesh& mesh,
        const incompressible::momentumTransportModel& turbulence,
        const solutionNMC& solution
    );

    // Public data

        const scalar kv_ = constant::mathematical::pi/6.0;

    // Member Functions

        //- Solve for the PBE
        virtual void transport_moments();

        //- Correct boundary condition of moments
        virtual void correct();

        //- update nodes and weights of the quadrature
        virtual void update_quadrature();
        
        // Return true if PBM uses quasi-high-order scheme
        virtual bool highOrderScheme() const = 0;

        //- Maximum allowable delta T in case of employing high order schemes
        virtual void momDeltaT(scalar& deltaT) const;

        //- Read pbProperties dictionary
        virtual bool read();

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the turbulence
            inline const incompressible::momentumTransportModel& turbulence()
                const;

            //- Return the growth model
            virtual const growthModel& growth() const = 0;

            //- Return the nucleation rate model
            virtual const nucleationRateModel& nucleationRate() const = 0;

            //- Return the nucleate size model
            virtual const nucleateSizeModel& nucleateSize() const = 0;

             //- Return the list of aggregation models
            virtual const aggregationList& aggregation() const = 0;

            //- Return the list of breakage models
            virtual const breakageList& breakage() const = 0;

            //- Return the quadrature model
            virtual const quadratureMethod& quadrature() const = 0;

            //- Return the phi field
            virtual const surfaceScalarField& phi() const = 0;

            //- Return the moments
            virtual const PtrList<volScalarField>& moments() const = 0;
            
            //- Return non-const access to the moments
            virtual PtrList<volScalarField>& moments() = 0;

            //- Calculate the precipitation rate
            virtual const DimensionedField<scalar, volMesh>&
            precRate() const = 0;

            //- Return the supersaturation
            virtual const DimensionedField<scalar, volMesh>& superSat()
                const = 0;

            //- Return the cation concentration ratios
            virtual const PtrList<DimensionedField<scalar, volMesh>>& 
                cationConcRatios() const = 0;

            //- Tagging cells for which source calculation and updating
            //  diameter are valid
            virtual const DimensionedField<scalar, volMesh>&
            sourceCorrCoeff() const = 0;

            //- Return the number of moments
            virtual uint8_t numOfMoments() const = 0;

            //- Return the number of nodes
            virtual uint8_t numOfNodes() const = 0;

            //- Return the number of coordinates
            virtual uint8_t numOfCoord() const = 0;

            //- Return the small diameters for the invalid cells
            virtual const List<scalar>& dSmall() const = 0;

            //- Return Minimum value for alpha used to check moment no. 3
            virtual const dimensionedScalar& alphaMin() const = 0;

            //- Return the crystal density
            virtual scalar crystalRho() const = 0;

            //- Return the crystal molecular weight
            virtual scalar crystalMW() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalanceI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
