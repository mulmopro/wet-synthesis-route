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

#include "wordIOList.H"
#include "aggregation.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::aggregation> Foam::aggregation::New
(
    const dictionary& dict,
    const populationBalance& pb
)
{

    word modelType;

    const int nCmpt = 2;
    const char* cmptNames[nCmpt] = {"kernel", "efficiency"};

    // Construct the name of the aggregation model from the components
    modelType =
        "aggregationModel<"
      + word(dict.lookup("type")) + ','
      + word(dict.subDict("efficiency").lookup("type")) + '>';

    Info<< "    " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
            << "Unknown aggregationModel type" << nl
            << modelType << nl << nl
            << "Valid aggregationModel types are:" << nl << nl;

            // Get the list of all the suitable thermo packages available
            wordList validModelTypes
            (
                dictionaryConstructorTablePtr_->sortedToc()
            );

            // Build a table of the aggregation model constituent parts
            // Note: row-0 contains the names of constituent parts
            List<wordList> validModelTypesCmpts
            (
                validModelTypes.size() + 1
            );

            validModelTypesCmpts[0].setSize(nCmpt);
            forAll(validModelTypesCmpts[0], j)
            {
                validModelTypesCmpts[0][j] = cmptNames[j];
            }

            // Split the aggregation model names into their constituent parts
            forAll(validModelTypes, i)
            {
                wordList cmpts(nCmpt);
            
                cmpts[0] =
                    validModelTypes[i].substr
                    (
                        validModelTypes[i].find('<', 0) + 1,
                        validModelTypes[i].find(',', 0) -
                        validModelTypes[i].find('<', 0) - 1
                    );
                   
                cmpts[1] =
                    validModelTypes[i].substr
                    (
                        validModelTypes[i].find(',', 0) + 1,
                        validModelTypes[i].find('>', 0) -
                        validModelTypes[i].find(',', 0) - 1
                    );
                
                    validModelTypesCmpts[i+1] = cmpts;
            }

            // Print the table of available models
            // in terms of their constituent parts
            printTable(validModelTypesCmpts, FatalError);

            FatalError<< exit(FatalError);
        }

    return autoPtr<aggregation>(cstrIter()(dict, pb));
}
