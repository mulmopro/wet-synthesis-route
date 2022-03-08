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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const Foam::fvMesh& Foam::solutionNMC::mesh() const
{
    return mesh_;
}


inline const Foam::List<Foam::word>& Foam::solutionNMC::metalNames() const
{
    return metalNames_;
}


inline Foam::label Foam::solutionNMC::nMetals() const
{
    return nMetals_;
}


inline Foam::scalar Foam::solutionNMC::effectiveConc() const
{
    return effectiveConc_;
}


inline Foam::PtrList<Foam::volScalarField>&
Foam::solutionNMC::totalMetalConcs()
{
    return totalMetalConcs_;
}


inline const Foam::PtrList<Foam::volScalarField>&
Foam::solutionNMC::totalMetalConcs() const
{
    return totalMetalConcs_;
}


inline const Foam::volScalarField& Foam::solutionNMC::totalNH3() const
{
    return totalNH3_;
}


inline const Foam::volScalarField& Foam::solutionNMC::totalNa() const
{
    return totalNa_;
}


inline const Foam::volScalarField& Foam::solutionNMC::totalSO4() const
{
    return totalSO4_;
}


inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::solutionNMC::metalCationConcs() const
{
    return metalCationConcs_;
}


inline const Foam::Field<Foam::scalar>& Foam::solutionNMC::activeCellsField()
const
{
    return activeCellsField_;
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::solutionNMC::NH3() const
{
    return NH3_;
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::solutionNMC::OH() const
{
    return OH_;
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::solutionNMC::pH() const
{
    return pH_;
}


inline const Foam::DimensionedField<Foam::scalar, Foam::volMesh>&
Foam::solutionNMC::superSat() const
{
    return superSat_;
}


inline const Foam::PtrList<Foam::DimensionedField<Foam::scalar, Foam::volMesh>>&
Foam::solutionNMC::cationConcRatios() const
{
    return cationConcRatios_;
}


inline const Foam::surfaceScalarField& Foam::solutionNMC::phi() const
{
    return phi_;;
}


inline const Foam::incompressible::momentumTransportModel&
Foam::solutionNMC::turbulence() const
{
    return turbulence_;
}

// ************************************************************************* //
