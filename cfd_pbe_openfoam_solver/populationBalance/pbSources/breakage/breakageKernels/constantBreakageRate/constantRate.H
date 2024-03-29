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

#ifndef constantBreakageRate_H
#define constantBreakageRate_H

#include "breakageKernel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace breakageKernels
{

/*---------------------------------------------------------------------------*\
                           Class constantRate Declaration
\*---------------------------------------------------------------------------*/

class constantRate
:
    private breakageKernel
{
    // Private data

        //- Constant breakage rate
        const dimensionedScalar C_;
        const scalar C_v_;

public:

    // Constructors

        //- Construct from a dictionary and a turbulence model
        constantRate
        (
            const dictionary& dict,
            const incompressible::momentumTransportModel& turbulence
        );


    //- Destructor
    virtual ~constantRate();


    //- Member Functions

        //- Breakage Frequency
        virtual tmp<DimensionedField<scalar, volMesh>> frequency
        (
            const DimensionedField<scalar, volMesh>&
        ) const;

        //- Breakage Frequency
        virtual scalar frequency(scalar, const PhysChemData&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace breakageKernels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "constantRate.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
