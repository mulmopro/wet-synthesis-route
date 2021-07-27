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

#ifndef aggregation_H
#define aggregation_H

#include "dictionary.H"
#include "populationBalance.H"
#include "runTimeSelectionTables.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aggregation Declaration
\*---------------------------------------------------------------------------*/

class aggregation
{
protected:
    
    // Protected data

        //- Reference to the population balance
        const populationBalance& pb_;


public:

    //- Runtime type information
    TypeName("aggregation");
    
    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            aggregation,
            dictionary,
            (
                const dictionary& dict,
                const populationBalance& pb
            ),
            (dict, pb)
        );

    // Constructors

        //- Construct from the model table, dictionary and pairs
        aggregation
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
        virtual ~aggregation();

    // Selectors

        static autoPtr<aggregation> New
        (
            const dictionary& dict,
            const populationBalance& pb
        );

    // Member Functions

        //- Return aggregation source for moment k
        virtual tmp<DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const = 0;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
