/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2019 OpenFOAM Foundation
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
    Foam::odeSolverList

Description
    List container for ode solvers

SourceFiles
    odeSolverList.C

\*---------------------------------------------------------------------------*/

#ifndef odeSolverList_H
#define odeSolverList_H

#include "IOdictionary.H"
#include "odeSolver.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class solutionNMC;
class populationBalance;
class envMixing;

/*---------------------------------------------------------------------------*\
                         Class odeSolverList Declaration
\*---------------------------------------------------------------------------*/

class odeSolverList
:
    public IOdictionary,
    public PtrList<odeSolver>
{
    // Private data

        const fvMesh& mesh_;


    // Private member functions

        // - Remove ambiguity by using size function of PtrList explicitly
        inline label size() const;

        // - Remove ambiguity by using set function of PtrList explicitly
        void set
        (
            label,
            const incompressible::momentumTransportModel&,
            const populationBalance&,
            const solutionNMC&,
            const envMixing&
        );

public:

    // Constructors

        odeSolverList
        (
            const fvMesh& mesh,
            const incompressible::momentumTransportModel& turbulence,
            const populationBalance& pb,
            const solutionNMC& solution,
            const envMixing& micromixing,
            label num_threads
        );

        //- Disallow default bitwise copy construction
        odeSolverList(const odeSolverList&) = delete;


    //- Destructor
    ~odeSolverList();


    // Member Functions

        // I-O

            //- Read dictionary
            bool read();


    // Member Operators
 
         //- Return element const reference
         inline const odeSolver& operator[](const label) const;
 
         //- Return element reference
         inline odeSolver& operator[](const label);

        //- Disallow default bitwise assignment
        void operator=(const odeSolverList&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "odeSolverListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
