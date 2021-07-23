/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "quadratureMethod.H"
#include "inversionAlgorithm.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(quadratureMethod, 0);
    defineRunTimeSelectionTable(quadratureMethod, populationBalance);
}

// * * * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::quadratureMethod::quadratureMethod(const populationBalance& pb)
:
    moments_(pb.moments()),
    numOfNodes_(pb.numOfNodes())
{
    wordList momPatchTypes = moments_[0].boundaryField().types();

    for (label count=0; count < numOfNodes_; count++)
    {
        word nodeName("X" + Foam::name(count));
        nodes_.append
        (
            new volScalarField
            (
                IOobject
                (
                    nodeName,
                    pb.mesh().time().timeName(),
                    pb.mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                pb.mesh(),
                dimensionedScalar("0", dimLength, 0),
                momPatchTypes
            )
        );
        nodes_[count].correctBoundaryConditions();

        word weightName("W" + Foam::name(count));
        weights_.append
        (
            new volScalarField
            (
                IOobject
                (
                    weightName,
                    pb.mesh().time().timeName(),
                    pb.mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                pb.mesh(),
                dimensionedScalar("0", dimless/dimVolume, 0),
                momPatchTypes
            )
        );
        weights_[count].correctBoundaryConditions();
    }

    invAlgm_ = inversionAlgorithm::New(pb, nodes_, weights_);

    invAlgm_->inversionBoundary();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::quadratureMethod::~quadratureMethod()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
