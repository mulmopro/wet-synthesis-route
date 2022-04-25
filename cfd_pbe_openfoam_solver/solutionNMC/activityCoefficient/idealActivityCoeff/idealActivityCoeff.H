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

#ifndef idealActivityCoeff_H
#define idealActivityCoeff_H

#include "activityCoeffModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionNMC;

namespace activityCoeffModels
{

/*---------------------------------------------------------------------------*\
                           Class idealActivityCoeff Declaration
\*---------------------------------------------------------------------------*/

class idealActivityCoeff
:
    public activityCoeffModel
{
    // private members


public:

    //- Runtime type information
    TypeName("idealSolution");


    // Constructors

        //- Construct from a dictionary and a solution
        idealActivityCoeff
        (
            const dictionary& dict,
            const solutionNMC& solution
        );


    //- Destructor
    virtual ~idealActivityCoeff();


    // Member Functions

        //- Calculate ionic strength
        virtual scalar ionic_strength(
            const List<scalar>&, const List<scalar>&) const;

        //- Calculate activity for a cation-OH pair
        virtual scalar pair_activity_coeff(
            label, const List<scalar>&, const List<scalar>&, scalar) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activityCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
