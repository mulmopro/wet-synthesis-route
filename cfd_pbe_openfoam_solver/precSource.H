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
            // cout<<"cell id: "<<celli<<'\n';

            ode.updateUserData(celli, solution_nmc, pb_ref, turb_ref);

            for(i=0; i<nMoments; i++)
            {
                y[i + nMetals] = moments[i][celli];
            }

            try
            {
                ode.solve(y, 0.0, odeEndTime);

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
                delete[] y;
                throw;
            }
        }

        delete[] y;
    }
}
