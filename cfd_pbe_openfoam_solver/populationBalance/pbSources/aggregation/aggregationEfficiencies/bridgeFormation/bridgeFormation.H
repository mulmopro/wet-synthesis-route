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

#ifndef bridgeFormation_H
#define bridgeFormation_H

#include "aggregationEfficiency.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace aggregationEfficiencies
{

/*---------------------------------------------------------------------------*\
                         Class bridgeFormation Declaration
\*---------------------------------------------------------------------------*/

class bridgeFormation
:
    public aggregationEfficiency
{
    // Private data

        //- Reference to turbulence model
        const incompressible::momentumTransportModel& turbulence_;

        //- Reference to growth model
        const growthModel& growth_;

        //- Yield strength of crystal
        const dimensionedScalar A_;
        const scalar A_v_;

        //- Density of liquid
        const dimensionedScalar rhoLiq_;
        const scalar rhoLiq_v_;


public:

    // Constructors

        //- Construct from a dictionary and a turbulence model
        bridgeFormation
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~bridgeFormation();


    // Member Functions

        //- Aggregation Efficiency
         virtual tmp<DimensionedField<scalar, volMesh>> efficiency
        (
            const DimensionedField<scalar, volMesh>&,
            const DimensionedField<scalar, volMesh>&
        ) const;

        virtual scalar efficiency
        (
            scalar, scalar, const PhysChemData&
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace aggregationEfficiencies
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "bridgeFormation.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
