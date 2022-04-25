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

#ifndef activityCoeffModel_H
#define activityCoeffModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "volFields.H"
#include "dictionary.H"
#include "runTimeSelectionTables.H"

#include "sundials_types.h"

namespace Foam
{

class solutionNMC;

/*---------------------------------------------------------------------------*\
                           Class activityCoeffModel Declaration
\*---------------------------------------------------------------------------*/

class activityCoeffModel
{
protected:

    // Protected data

    const fvMesh& mesh_; 

    const solutionNMC& solution_;


public:

    //- Runtime type information
    TypeName("activityCoeffModel");


    // Declare runtime construction
    declareRunTimeSelectionTable
    (
        autoPtr,
        activityCoeffModel,
        dictionary,
        (
            const dictionary& dict,
            const solutionNMC& solution
        ),
        (dict, solution)
    );


    // Constructors

        //- Construct from a dictionary and a solution
        activityCoeffModel
        (
            const dictionary& dict,
            const solutionNMC& solution
        );


    //- Destructor
    virtual ~activityCoeffModel();


    // Selectors

        static autoPtr<activityCoeffModel> New
        (
            const dictionary& dict,
            const solutionNMC& solution
        );


    // Member Functions

        //- Calculate ionic strength (if needed by the selected model)
        virtual scalar ionic_strength(
            const List<scalar>&, const List<scalar>&) const = 0;

        //- Calculate activity for a cation-OH pair
        virtual scalar pair_activity_coeff(
            label, const List<scalar>&, const List<scalar>&, scalar) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
