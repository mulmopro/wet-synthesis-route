/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "odeSolverList.H"

#include "solutionNMC.H"
#include "populationBalance.H"
#include "envMixing.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::odeSolverList::odeSolverList
(
    const fvMesh& mesh,
    const incompressible::momentumTransportModel& turbulence,
    const populationBalance& pb,
    const solutionNMC& solution,
    const envMixing& micromixing,
    label num_threads
)
:
    IOdictionary
    (
        IOobject
        (
            "odeSolver",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    PtrList<odeSolver>(num_threads),

    mesh_(mesh)
{
    forAll(*this, index)
    {
        set(index, turbulence, pb, solution, micromixing);
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::odeSolverList::~odeSolverList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::odeSolverList::set
(
    label index,
    const incompressible::momentumTransportModel& turbulence,
    const populationBalance& pb,
    const solutionNMC& solution,
    const envMixing& micromixing
)
{
    PtrList<odeSolver>::set
    (
        index,
        new odeSolver(mesh_, turbulence, pb, solution, micromixing, *this)
    );
}


bool Foam::odeSolverList::read()
{
    if (regIOobject::read())
    {
        bool allOk = true;

        forAll(*this, i)
        {
            odeSolver& ode = operator[](i);

            bool ok = ode.read(*this);

            allOk = (allOk && ok);
        }

        return allOk;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
