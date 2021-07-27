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

#ifndef breakage_H
#define breakage_H

#include "populationBalance.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class breakage Declaration
\*---------------------------------------------------------------------------*/

class breakage
{
protected:
    // Protected data

        //- Reference to the population balance
        const populationBalance& pb_;


public:

    //- Runtime type information
    TypeName("breakage");
    
    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            breakage,
            dictionary,
            (
                const dictionary& dict,
                const populationBalance& pb
            ),
            (dict, pb)
        );

    // Constructors

        //- Construct from the model table, dictionary and pairs
        breakage
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
        virtual ~breakage();

    // Selectors

        static autoPtr<breakage> New
        (
            const dictionary& dict,
            const populationBalance& pb
        );

    // Member Functions

        //- Return breakage source for moment k
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
