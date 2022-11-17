/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2015 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::envMixing::mesh() const
{
    return mesh_;
}


inline Foam::label Foam::envMixing::nPureEnvs() const
{
    return nPureEnvs_;
}


inline const Foam::volScalarField& Foam::envMixing::reactingEnvP() const
{
    return reactingEnvP_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::envMixing::pureEnvP() const
{
    return pureEnvP_;
}


inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::envMixing::pFluxes() const
{
    return pFluxes_;
}


inline const Foam::PtrList<Foam::dimensionedScalar>&
Foam::envMixing::envConcs() const
{
    return envConcs_;
}


inline const Foam::surfaceScalarField& Foam::envMixing::phi() const
{
    return phi_;
}


inline const Foam::incompressible::momentumTransportModel&
Foam::envMixing::turbulence() const
{
    return turbulence_;
}

// ************************************************************************* //
