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

#ifndef bromleyActivityCoeff_H
#define bromleyActivityCoeff_H

#include "activityCoeffModel.H"
#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class solutionNMC;

namespace activityCoeffModels
{

/*---------------------------------------------------------------------------*\
                           Class bromleyActivityCoeff Declaration
\*---------------------------------------------------------------------------*/

class bromleyActivityCoeff
:
    public activityCoeffModel
{
    //- Private members

    const scalar A_gamma_;

    label nCation_;
    label nAnion_;

    scalarRectangularMatrix B_;
    scalarRectangularMatrix E_;

    List<scalar> z_c_;
    List<scalar> z_a_;

    DimensionedField<scalar, volMesh> I_;
    
    PtrList<DimensionedField<scalar, volMesh>> pairActivityCoeff_;

    //- Calculate activity for a cation-OH pair
    scalar pure_pair_activity_coeff(
        label id_cation, label id_anion, scalar z_c, scalar z_a, scalar I_s);

    
public:

    //- Runtime type information
    TypeName("Bromley");


    // Constructors

        //- Construct from a dictionary and a solution
        bromleyActivityCoeff
        (
            const dictionary& dict,
            const solutionNMC& solution
        );


    //- Destructor
    virtual ~bromleyActivityCoeff();


    // Member Functions

        //- Calculate ionic strength
        virtual void update_ionic_strength(
            const List<scalar>&, const List<scalar>&, label);

        //- Calculate activity for a cation-OH pair
        virtual scalar pair_activity_coeff(
            label, const List<scalar>&, const List<scalar>&, label);

        //- Get pre-calculated activity coefficient
        virtual scalar pair_activity_coeff(label, label);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activityCoeffModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
