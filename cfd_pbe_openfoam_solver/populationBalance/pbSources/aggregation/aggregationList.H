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

#ifndef aggregationList_H
#define aggregationList_H

#include "aggregation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class aggregationList Declaration
\*---------------------------------------------------------------------------*/

class aggregationList
:
    PtrList<aggregation>
{
private:

    // Private Member Functions

        //- Disallow default bitwise copy construct
        aggregationList(const aggregationList&);

        //- Disallow default bitwise assignment
        void operator=(const aggregationList&);


protected:

    // Protected data

        //- Reference to the population balance
        const populationBalance& pb_;


public:

    //- Constructor
    aggregationList
    (
        const dictionary& dict,
        const populationBalance& pb
    );

    //- Destructor
    ~aggregationList();


    // Member Functions

        //- Calculate the aggregation models
        tmp<DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const;

        //- Calculate the aggregation models
        scalar source
        (
            const List<scalar>&,
            const List<scalar>&,
            int,
            const PhysChemData&
        ) const;

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
