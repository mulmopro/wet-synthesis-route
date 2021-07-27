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

#ifndef classicalHeteroHomoNucleation_H
#define classicalHeteroHomoNucleation_H

#include "nucleationRateModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class populationBalance;

namespace nucleationRateModels
{

/*---------------------------------------------------------------------------*\
                           Class classicalHeteroHomoNucleation Declaration
\*---------------------------------------------------------------------------*/

class classicalHeteroHomoNucleation
:
    public nucleationRateModel
{
    //- private members

    const DimensionedField<scalar, volMesh>& superSat_;

    dimensionedScalar rateDimension_;
    
    scalar            KJ_1_;

    scalar            KJ_2_;
    
    scalar            BJ_1_;

    scalar            BJ_2_;


public:

    //- Runtime type information
    TypeName("classicalHeteroHomo");

    // Constructors

        //- Construct from a dictionary and a population balance model
        classicalHeteroHomoNucleation
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~classicalHeteroHomoNucleation();


    // Member Functions

        //- Constant nucleation rate
        tmp<DimensionedField<scalar, volMesh>> rate() const;

        realtype rate(realtype superSat) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace nucleationRateModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
