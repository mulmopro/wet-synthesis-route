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

#ifndef breakageModel_H
#define breakageModel_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "breakage.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class breakageModel Declaration
\*---------------------------------------------------------------------------*/

template<class BKernel, class BDaughterDist>
class breakageModel
:
    public breakage,
    public BKernel, 
    public BDaughterDist
{

    // Private Member Functions

        //- Disallow default bitwise copy construct
        breakageModel(const breakageModel&);

        //- Disallow default bitwise assignment
        void operator=(const breakageModel&);


public:

    //- Runtime type information
    TypeName("breakageModel");


    // Constructors

        //- Construct from components
        breakageModel
        (
            const dictionary& dict,
            const populationBalance& pb
        );


    //- Destructor
    virtual ~breakageModel()
    {}


    // Member Functions

        //- Return breakage source for moment k
        virtual tmp<Foam::DimensionedField<scalar, volMesh>> source
        (
            const PtrList<volScalarField>&,
            const PtrList<volScalarField>&,
            int
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "breakageModel.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
