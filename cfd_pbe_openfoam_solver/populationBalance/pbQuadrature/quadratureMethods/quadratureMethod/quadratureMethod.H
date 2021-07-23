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

Class
    Foam::quadratureMethod

Description

SourceFiles
    quadratureMethod.C
    newQuadratureMethod.C

\*---------------------------------------------------------------------------*/

#ifndef quadratureMethod_H
#define quadratureMethod_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "populationBalance.H"
#include "runTimeSelectionTables.H"

namespace Foam
{

class inversionAlgorithm;

/*---------------------------------------------------------------------------*\
                           Class quadratureMethod Declaration
\*---------------------------------------------------------------------------*/

class quadratureMethod
{
protected:

    // Protected data

        //- Reference to the moments
        const PtrList<volScalarField>& moments_;

        //- Reference to the number of nodes
        const uint8_t& numOfNodes_;

        //- Nodes of the NDF
        PtrList<volScalarField> nodes_;

        //- Weights of the NDF
        PtrList<volScalarField> weights_;

        //- Inversion Algorithm
        autoPtr<inversionAlgorithm> invAlgm_;

public:

    //- Runtime type information
    TypeName("quadratureMethod");


    // Declare runtime construction

        declareRunTimeSelectionTable
        (
            autoPtr,
            quadratureMethod,
            populationBalance,
            (const populationBalance& pb),
            (pb)
        );

    // Constructors

        quadratureMethod(const populationBalance& pb);


    //- Destructor
    virtual ~quadratureMethod();


    // Selectors

        static autoPtr<quadratureMethod> New(const populationBalance& pb);


    // Member Functions

        //- Calculates nodes and weights of the reconstructed NDF
        virtual void nodesAndWeights(const Field<scalar>&) = 0;

        //- Returns const access to the nodes
        const PtrList<volScalarField>& nodes() const
        {
            return nodes_;
        }

        //- Returns non-const access to the nodes
        PtrList<volScalarField>& nodes()
        {
            return nodes_;
        }

        //- Returns const access to the weights
        const PtrList<volScalarField>& weights() const
        {
            return weights_;
        }

        //- Returns non-const access to the weights
        PtrList<volScalarField>& weights()
        {
            return weights_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
