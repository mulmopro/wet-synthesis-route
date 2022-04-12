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

#ifndef UserData_H
#define UserData_H

#include "sundials_types.h"
#include "solutionNMC.H"
#include "populationBalance.H"
#include "quadratureMethod.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

struct PhysChemData
{
    scalar epsilon;
    const scalar nu;
    scalar growthRate;
};

class growthModel;
class nucleationRateModel;
class nucleateSizeModel;
class aggregationList;
class breakageList;
class inversionAlgorithm;


/*---------------------------------------------------------------------------*\
                      Struct userData Declaration
\*---------------------------------------------------------------------------*/

struct UserData{
    UserData
    (
        label celli,
        solutionNMC& solution,
        const populationBalance& pb,
        const incompressible::momentumTransportModel& turbulence
    ): solution_(solution), growth_(pb.growth()), nucRate_(pb.nucleationRate()),
    nucleateSize_(pb.nucleateSize()), aggregation_(pb.aggregation()),
    breakage_(pb.breakage()), invAlgm_(pb.quadrature().invAlgorithm()),
    cell_id(celli), nMetals(solution.nMetals()), nMoments(pb.numOfMoments()),
    totalNH3(solution.totalNH3()[celli]), totalNa(solution.totalNa()[celli]),
    totalSO4(solution.totalSO4()[celli]),
    effectiveConc(solution.effectiveConc()), kv(pb.kv_),
    crystalRho(pb.crystalRho()), crystalMW(pb.crystalMW()), superSat(0.0),
    concNH3(solution.NH3()[celli]), concOH(solution.OH()[celli]),
    physChemData_{
        turbulence.epsilon()().internalField()[celli],
        turbulence.nu()().internalField()[celli],
        0.0
    }
    {
        for(int i=0; i<nMetals; i++)
        {
            metalCationConcs[i] = solution.metalCationConcs()[i][celli];
        }
    }

    solutionNMC& solution_;
    const growthModel& growth_;
    const nucleationRateModel& nucRate_;
    const nucleateSizeModel& nucleateSize_;
    const aggregationList& aggregation_;
    const breakageList& breakage_;
    const inversionAlgorithm& invAlgm_;
    label cell_id;
    const label nMetals, nMoments;
    realtype totalNH3, totalNa, totalSO4;
    const realtype effectiveConc, kv, crystalRho, crystalMW;
    realtype superSat, concNH3, concOH;
    realtype metalCationConcs[3];

    PhysChemData physChemData_;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
