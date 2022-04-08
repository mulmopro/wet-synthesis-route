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

#ifndef aggregationModel_H
#define aggregationModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "aggregation.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class aggregationModel Declaration
\*---------------------------------------------------------------------------*/

template<class CKernel, class CEfficiency>
class aggregationModel
:
    public aggregation,
    public CKernel, 
    public CEfficiency
{

public:

    //- Runtime type information
    TypeName("aggregationModel");


    // Constructors

        //- Construct from components
        aggregationModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
        virtual ~aggregationModel()
        {}


    // Member Functions

        //- Return aggregation source for moment k
        virtual tmp<Foam::DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const;

        //- Return aggregation source for moment k
        virtual scalar source
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

#ifdef NoRepository
    #include "aggregationModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
