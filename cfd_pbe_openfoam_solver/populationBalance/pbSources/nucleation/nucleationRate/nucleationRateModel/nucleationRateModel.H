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

#ifndef nucleationRateModel_H
#define nucleationRateModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "sundials_types.h"
namespace Foam
{

class populationBalance;

/*---------------------------------------------------------------------------*\
                           Class nucleationRateModel Declaration
\*---------------------------------------------------------------------------*/

class nucleationRateModel
{
protected:

    // Protected data

    const fvMesh& mesh_; 

    const populationBalance& pb_;


public:

    //- Runtime type information
    TypeName("nucleationRateModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        nucleationRateModel,
        dictionary,
        (
            const dictionary& dict,
            const populationBalance& pb
        ),
        (dict, pb)
    );


    // Constructors

        //- Construct from a dictionary and a population balance
        nucleationRateModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~nucleationRateModel();


    // Selectors

        static autoPtr<nucleationRateModel> New
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    // Member Functions

        //- Nucleation rate
        virtual tmp<DimensionedField<scalar, volMesh>> rate() const = 0;

        virtual realtype rate(realtype superSat) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
