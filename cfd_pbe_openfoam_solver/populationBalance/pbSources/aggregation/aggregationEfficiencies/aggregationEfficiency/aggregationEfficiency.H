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

#ifndef aggregationEfficiency_H
#define aggregationEfficiency_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


namespace Foam
{

struct PhysChemData;

/*---------------------------------------------------------------------------*\
                           Class swarmCorrection Declaration
\*---------------------------------------------------------------------------*/

class aggregationEfficiency
{
protected:

    // Protected data

        //- Reference to the population balance
        const populationBalance& pb_;

public:

    // Constructors

        //- Construct from a phase pair
        aggregationEfficiency
        (
            const populationBalance& pb
        );


    //- Destructor
    virtual ~aggregationEfficiency();


    // Member Functions

        virtual tmp<DimensionedField<scalar, volMesh>> efficiency
        (
            const DimensionedField<scalar, volMesh>&,
            const DimensionedField<scalar, volMesh>&
        ) const = 0;

        virtual scalar efficiency
        (
            scalar, scalar, const PhysChemData&
        ) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "aggregationEfficiency.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
