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

#include "wordIOList.H"
#include "breakage.H"

// * * * * * * * * * * * * * * * * Selector  * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::breakage> Foam::breakage::New
(
    const dictionary& dict,
    const populationBalance& pb
)
{

    word modelType;

    const int nCmpt = 2;
    const char* cmptNames[nCmpt] = {"kernel", "daughterEfficiency"};

    // Construct the name of the breakage model from the components
    modelType =
        "breakageModel<"
      + word(dict.lookup("type")) + ','
      + word(dict.subDict("daughterDistribution").lookup("type")) + '>';

    Info<< "    " << modelType << endl;

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);
    
    if (cstrIter == dictionaryConstructorTablePtr_->end())
        {
            FatalErrorInFunction
            << "Unknown breakageModel type" << nl
            << modelType << nl << nl
            << "Valid breakageModel types are:" << nl << nl;

            // Get the list of all the suitable thermo packages available
            wordList validModelTypes
            (
                dictionaryConstructorTablePtr_->sortedToc()
            );

            // Build a table of the breakage model constituent parts
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

            // Split the breakage model names into their constituent parts
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
    
    return autoPtr<breakage>(cstrIter()(dict, pb));
}
