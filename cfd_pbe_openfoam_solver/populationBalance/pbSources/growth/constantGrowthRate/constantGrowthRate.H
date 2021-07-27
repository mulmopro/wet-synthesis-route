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

#ifndef constantGrowthRate_H
#define constantGrowthRate_H

#include "growthModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class populationBalance;

namespace growthModels
{

/*---------------------------------------------------------------------------*\
                           Class constantGrowthRate Declaration
\*---------------------------------------------------------------------------*/

class constantGrowthRate
:
    public growthModel
{
    //- Private members

    const DimensionedField<scalar, volMesh>& superSat_;
    
    dimensionedScalar G0_;

    
public:

    //- Runtime type information
    TypeName("constantRate");


    // Constructors

        //- Construct from a dictionary and a population balance model
        constantGrowthRate
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~constantGrowthRate();


    // Member Functions

        //- Constant growth rate
        tmp<DimensionedField<scalar, volMesh>>
            rate(const DimensionedField<scalar, volMesh>&) const;

        realtype rate(realtype superSat) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace growthModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
