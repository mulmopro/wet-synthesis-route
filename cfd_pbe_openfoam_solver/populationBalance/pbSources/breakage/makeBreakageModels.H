/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2015 OpenFOAM Foundation
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

#ifndef makeBreakageModels_H
#define makeBreakageModels_H

#include "breakageModel.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// This function is used to define all possible combinations of kernels and
// daughter distributions for the breakage model
#define makeBreakageModel(Kernel, daughterDist)                               \
                                                                              \
    namespace Foam                                                            \
    {                                                                         \
        typedef                                                               \
        breakageModel                                                         \
        <                                                                     \
            breakageKernels::Kernel,                                          \
            daughterDistributions::daughterDist                               \
        > breakageModel##Kernel##daughterDist;                                \
                                                                              \
        defineTemplateTypeNameAndDebugWithName                                \
        (                                                                     \
            breakageModel##Kernel##daughterDist,                              \
            "breakageModel<"#Kernel","#daughterDist">",                       \
            0                                                                 \
        );                                                                    \
                                                                              \
        addToRunTimeSelectionTable                                            \
        (                                                                     \
            breakage,                                                         \
            breakageModel##Kernel##daughterDist,                              \
            dictionary                                                        \
        );                                                                    \
    }

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    
#endif

// ************************************************************************* //


