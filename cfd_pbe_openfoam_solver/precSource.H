const label nMetals(solution_nmc.nMetals());
const label nMoments(pb_ref.numOfMoments());
const label nY = nMetals + nMoments;

const scalar effectiveConc(solution_nmc.effectiveConc());

const realtype odeEndTime(runTime.deltaTValue());

forAll(activeCellsField, celli)
{
    if (activeCellsField[celli] > 0)
    {
        int i;
        realtype *y = new realtype[nY];

        scalar totalConc;
        scalar cationTotalConc(0.0);
        for(i=0; i<nMetals; i++)
        {
            totalConc = totalMetalConcs[i][celli];

            cationTotalConc += totalConc;
            y[i] = totalConc;
        }

        if (cationTotalConc > effectiveConc)
        {
            // Info<<"cell id: "<<celli<<endl;
            UserData *aux_data = new UserData
            (
                celli,
                solution_nmc,
                pb_ref
            );

            for(i=0; i<nMoments; i++)
            {
                y[i + nMetals] = moments[i][celli];
            }

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
