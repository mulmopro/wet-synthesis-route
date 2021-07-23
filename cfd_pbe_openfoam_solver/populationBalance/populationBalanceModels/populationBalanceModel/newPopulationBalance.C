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

#include "populationBalance.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::populationBalance> Foam::populationBalance::New
(
    const fvMesh& mesh,
    const incompressible::momentumTransportModel& turbulence,
    const solutionNMC& solution
)
{
    // Access th pbProperties, but do not register the dictionary
    // otherwise it is registered in the database twice
    IOdictionary pbDict
    (
        IOobject
        (
            "pbProperties",
            mesh.time().constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE,
            false
        )
    );

    word PBMType;

    if (pbDict.lookupOrDefault<Switch>("populationBalance", false))
    {
        PBMType = "constantVelocity";
    }
    else
    {
        PBMType = "none";
    }

    Switch HOScheme(pbDict.lookupOrDefault<Switch>("highOrderScheme", false));

    word schemeType = "";

    if (HOScheme)
    {
        schemeType = "HO";

        Info<< "Selecting Population Balance Model: " << PBMType
            << " (Quasi-high-order schemes enabled)" << endl;
    }
    else
    {
        Info<< "Selecting Population Balance Model: " << PBMType << endl;
    }

    fvMeshConstructorTable::iterator cstrIter =
        fvMeshConstructorTablePtr_->find(PBMType + schemeType);

    if (cstrIter == fvMeshConstructorTablePtr_->end())
    {
        wordList validModelTypes
        (
            fvMeshConstructorTablePtr_->sortedToc()
        );

        cstrIter = fvMeshConstructorTablePtr_->find(PBMType);

        if (cstrIter == fvMeshConstructorTablePtr_->end())
        {
            FatalErrorInFunction
                << "Unknown Population Balance Model "
                << PBMType << endl << endl
                << "Valid Population Balance Models are : " << endl;
            
            forAll(validModelTypes, i)
            {
                word validModelType = validModelTypes[i];

                if
                (
                    validModelType.length() > 2
                 && validModelType.substr(validModelType.length() - 2)=="HO"
                )
                {
                    label elementIndex(-1);

                    forAll(validModelTypes, j)
                    {
                        if
                        (
                            validModelTypes[j]
                         == validModelType.substr(0, validModelType.length() - 2)
                        )
                        {
                            elementIndex = j;
                            break;
                        }
                    }

                    if (elementIndex < 0)
                    {
                        FatalErrorInFunction
                            << validModelType.substr(0, validModelType.length() - 2)
                            << endl;
                    }
                }
                else
                {
                    FatalErrorInFunction
                        << validModelType << endl;
                }
            }

            FatalError << exit(FatalError);
        }
        else
        {
            FatalErrorInFunction
                << "Quasi-high-order schemes are not implemented "
                << "in the Population Balance Model "
                << PBMType << endl << endl
                << "Quasi-high-order schemes can be used with the following "
                << "Population Balance Models: " << endl;

            forAll(validModelTypes, i)
            {
                word validModelType = validModelTypes[i];

                if
                (
                    validModelType.length() > 2
                 && validModelType.substr(validModelType.length() - 2)=="HO"
                )
                {
                    FatalErrorInFunction
                        << validModelType.substr(0, validModelType.length() - 2)
                        << endl;
                }
            }

            FatalError << exit(FatalError);
        }
    }

    return cstrIter()(mesh, turbulence, solution);
}


// ************************************************************************* //
