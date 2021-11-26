const label nMetals(solution_nmc.nMetals());
const label nMoments(pb_ref.numOfMoments());
const label nY = nMetals + nMoments;

const scalar kv(pb_ref.kv_);
const scalar crystalRho(pb_ref.crystalRho());
const scalar crystalMW(pb_ref.crystalMW());

const scalar effectiveConc(solution_nmc.effectiveConc());

const realtype odeEndTime(runTime.deltaTValue());

forAll(activeCellsField, celli)
{
    if (activeCellsField[celli] > 0)
    {
        int i;
        scalar value;
        realtype *y = new realtype[nY];
        Switch validConc(true);

        for(i=0; i<nMetals; i++)
        {
            value = totalMetalConcs[i][celli];

            if (value < effectiveConc)
            {
                validConc = false;
                break;
            }

            y[i] = value;
        }

        if (validConc)
        {
            // Info<<"cell id: "<<celli<<endl;
            UserData *aux_data = new UserData
            (
                solution_nmc,
                pb_ref.growth(),
                pb_ref.nucleationRate(),
                pb_ref.nucleateSize()
            );

            for(i=0; i<nMoments; i++)
            {
                y[i + nMetals] = moments[i][celli];
            }

            for(i=0; i<nMetals; i++)
            {
                aux_data->metalCationConcs[i] = metalCationConcs[i][celli];
            }

            aux_data->cell_id = celli;
            aux_data->nMetals = nMetals;
            aux_data->nMoments = nMoments;
            aux_data->kv = kv;
            aux_data->crystalRho = crystalRho;
            aux_data->crystalMW = crystalMW;
            aux_data->totalNH3 = totalNH3[celli];
            aux_data->inertCharges = inertCharges[celli];
            aux_data->concNH3 = NH3[celli];
            aux_data->concOH = OH[celli];
            aux_data->superSat = 0.0;
            aux_data->effectiveConc = effectiveConc;

            try
            {
                ode.solve(y, 0.0, odeEndTime, aux_data);

                for(i = 0; i < nMetals; i++)
                {
                    totalMetalConcs[i][celli] = y[i];
                }

                for(i = 0; i < nMoments; i++)
                {
                    moments[i][celli] = y[i + nMetals];
                }
            }
            catch(...)
            {
                delete aux_data;
                delete[] y; 
                throw;
            }

            delete aux_data;
        }

        delete[] y;
    }
}
