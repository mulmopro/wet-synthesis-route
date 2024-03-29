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

#ifndef solutionNMC_H
#define solutionNMC_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "geometricZeroField.H"

#include "kinematicMomentumTransportModel.H"

#include "sundials_types.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct UserData;
class activityCoeffModel;
class envMixing;

/*---------------------------------------------------------------------------*\
                      Class solutionNMC Declaration
\*---------------------------------------------------------------------------*/

class solutionNMC
:
    public IOdictionary
{

    // Private data

        //- Reference to the mesh
        const fvMesh& mesh_;

        //- Reference to the flux at faces
        const surfaceScalarField& phi_;

        //- Reference to the turbulence
        const incompressible::momentumTransportModel& turbulence_;

        //- Reference to the micromixing model
        const envMixing& micromixing_;

        //- List of the metals
        List<word> metalNames_;

        //- Number of metal cations
        label nMetals_;

        //- Number of total concentrations
        label nTotalConc_;

        //- Number of equilibrium unknowns
        label nComps_;

        //- index of NH3
        label indexNH3_;

        //- index of OH
        label indexOH_;

        //- index of Na
        label indexNa_;

        //- index of SO4
        label indexSO4_;

        //- Maximum number of iterations for Newton-Raphson
        label maxIter_;

        //- Maximum number of iterations for Newton-Raphson
        scalar tol_;

        //- Minimum value for concentrations
        scalar smallConc_;

        //- Constant value for the temperature
        scalar T_;

        //- Minimum value for effective concentrations
        scalar effectiveConc_;

        //- Interval for writing the summary of concentrations in the log
        label writeSummaryInterval_;

        //- Turbulent Schmidt number
        dimensionedScalar turbSc_;

        // Used to tag chemically active cells
        DimensionedField<scalar, volMesh> activeCells_;
        const Field<scalar>& activeCellsField_;

        //- Equilibrium concentration of NH3
        DimensionedField<scalar, volMesh> NH3_;

        //- Equilibrium concentration of OH
        DimensionedField<scalar, volMesh> OH_;

        //- superSaturation
        DimensionedField<scalar, volMesh> superSat_;

        //- pH
        DimensionedField<scalar, volMesh> pH_;

        //- Total concentration of NH3
        volScalarField totalNH3_;

        //- Total concentration of Na
        volScalarField totalNa_;

        //- Total concentration of SO4
        volScalarField totalSO4_;

        //- Pointer list of total concentration of metals
        PtrList<volScalarField> totalMetalConcs_;

        //- Pointer list of equilibrium concentration of metal ions
        PtrList<DimensionedField<scalar, volMesh>> metalCationConcs_;

        //- Pointer list of ratio of metals to total metal concentration
        PtrList<DimensionedField<scalar, volMesh>> cationConcRatios_;

        //- List of molecular diffusivities
        PtrList<dimensionedScalar> D_;

        //- activity coefficient model
        autoPtr<activityCoeffModel> activityModel_;

        //- Dissociation constant of water
        const scalar pKw_ = -0.000113*Foam::pow((T_-273.15), 2) + 0.0371*(T_-273.15) - 14.8;
        const scalar kw_ = Foam::pow(10, pKw_);

        //- Dissociation constant of ammonia
        const scalar pKb_NH3_ = -0.0000422*Foam::pow((T_-273.15), 2) + 0.0038*(T_-273.15) - 4.82; 
        const scalar Kb_NH3_ = Foam::pow(10, pKb_NH3_);

        const scalar ln10_ = Foam::log(10.0);

        //- ksp of metal hydroxides
        const FixedList<scalar, 3> k_sp_ =
            FixedList<scalar, 3>
            (
                Foam::pow
                (
                    10,
                    List<scalar>
                    (
                        {-15.22, -12.70, -14.89}
                    )
                )
            );

        //- Number of ammonia complexes of metals
        const FixedList<label, 3> numOfComplexes_ = {6, 4, 6};

        //- ksp of metal hydroxides
        const FixedList<scalar, 16> kn_NMC_ =
            FixedList<scalar, 16>
            (
                Foam::pow
                (
                    10,
                    List<scalar>
                    (
                        {
                            2.81, 5.08, 6.85, 8.12, 8.93, 9.08,
                            1.00, 1.54, 1.70, 1.30,
                            2.10, 3.67, 4.78, 5.53, 5.75, 5.14
                        }
                    )
                )
            );

        //- Dimension of concentration
        const dimensionSet dimConc = dimMoles/dimVol;

    // Private member functions

        //- Calculate the matrices of the chemical equilibria equations
        void equilibriaEqs(const List<scalar>&, const List<scalar>&,
            scalar, double*, double*) const;

public:

    // Constructors

        //- Construct from fvMesh
        solutionNMC
        (
            const fvMesh&,
            const surfaceScalarField&,
            const incompressible::momentumTransportModel&,
            const envMixing&
        );


    //- Destructor
    virtual ~solutionNMC();


    // Member Functions

        //- Find equilibrium concentrations to calculate the supersaturation 
        bool update(label);

        //- solve quilibrium for a cell
        void solve
        (
            const realtype *y, UserData *aux_data, realtype *cationConcRatios
        ) const;

        //- Correct boundary condition of total concentrations
        void correct();

        //- Solve for the total concentrations
        void transport_species();

        //- Read solutionProperties dictionary
        bool read();

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the metal names
            inline const List<word>& metalNames() const;

            //- Return the number of metals
            inline label nMetals() const;

            //- Return the effective concentration
            inline scalar effectiveConc() const;

            //- Return the total concentration of metals
            inline const PtrList<volScalarField>& totalMetalConcs() const;

            //- Return non-const access to the total concentration of metals
            inline PtrList<volScalarField>& totalMetalConcs();

            //- Return the total concentration of NH3
            inline const volScalarField& totalNH3() const;

            //- Return the total concentration of inertCharges
            inline const volScalarField& totalNa() const;

            //- Return the total concentration of inertCharges
            inline const volScalarField& totalSO4() const;

            //- Return the equilibrium concentration of metal ions
            inline const PtrList<DimensionedField<scalar, volMesh>>&
            metalCationConcs() const;

            //- Return the chemically active cells
            inline const Field<scalar>& activeCellsField() const;

            //- Return the equilibrium concentration of NH3
            inline const DimensionedField<scalar, volMesh>& NH3() const;

            //- Return the equilibrium concentration of OH
            inline const DimensionedField<scalar, volMesh>& OH() const;

            //- Return the supersaturation
            inline const DimensionedField<scalar, volMesh>& superSat() const;

            //- Return the supersaturation
            inline const PtrList<DimensionedField<scalar, volMesh>>&
                cationConcRatios() const;

            //- Return the pH
            inline const DimensionedField<scalar, volMesh>& pH() const;

            //- Return the flux at faces
            inline const surfaceScalarField& phi() const;

            //- Return the turbulence
            inline const incompressible::momentumTransportModel& turbulence()
                const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solutionNMCI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
