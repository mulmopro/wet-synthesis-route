/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::envMixing

Description

SourceFiles
    envMixing.C

\*---------------------------------------------------------------------------*/

#ifndef envMixing_H
#define envMixing_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "geometricZeroField.H"

#include "kinematicMomentumTransportModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class envMixing Declaration
\*---------------------------------------------------------------------------*/

class envMixing
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

        //- Number of pure environments
        label nPureEnvs_;

        //- Correction coefficient for the mixing rate
        scalar correctionCoeff_;

        //- Turbulent Schmidt number
        dimensionedScalar turbSc_;

        //- Probability of reacting environment
        volScalarField reactingEnvP_;

        //- Pointer list of pure environment probabilities
        PtrList<volScalarField> pureEnvP_;

        //- Pointer list of probability fluxes
        PtrList<DimensionedField<scalar, volMesh>> pFluxes_;

        //- correlation coefficients
        const FixedList<scalar, 7> a_cphi = {
            0.4093, 0.6015, 0.5851, 0.09472, -0.3903, 0.1461, -0.01604};


    // Private member functions

        //- Calculate gamma
        const tmp<DimensionedField<scalar, volMesh>> gamma_calc() const;

        //- Calculate the mixing rates for each environment
        void p_flux();

public:

    // Constructors

        //- Construct from fvMesh
        envMixing
        (
            const fvMesh&,
            const surfaceScalarField&,
            const incompressible::momentumTransportModel&
        );


    //- Destructor
    virtual ~envMixing();


    // Member Functions

        //- Solve for the environments
        void solve();

        //- Read micromixing dictionary
        bool read();

        // Access

            //- Return the mesh
            inline const fvMesh& mesh() const;

            //- Return the number of environments
            inline label nPureEnvs() const;

            //- Return the probability of reacting environment
            inline const volScalarField& reactingEnvP() const;

            //- Return the probability of pure environments
            inline const PtrList<volScalarField>& pureEnvP() const;

            //- Return the flux at faces
            inline const surfaceScalarField& phi() const;

            //- Return the turbulence
            inline const incompressible::momentumTransportModel& turbulence()
                const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "envMixingI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
