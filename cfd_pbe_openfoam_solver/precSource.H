const realtype odeEndTime(runTime.deltaTValue());

int error_flag = 0;

#ifdef _OPENMP
#pragma omp parallel default(none) \
        shared(nMetals, nMoments, nY, nCells, odeEndTime, activeCellsField, \
            totalMetalConcs, moments, odes, solution_nmc, error_flag)
{
    odeSolver& ode = odes[omp_get_thread_num()];
#else
    odeSolver& ode = odes[0];
#endif

    int i;
    bool validTotalConc;

    List<realtype> y(nY);

    #ifdef _OPENMP
    #pragma omp for schedule(dynamic, 4)
    #endif
    for (int celli=0; celli<nCells; celli++)
    {
        if (activeCellsField[celli] > 0)
        {
            // update the supersaturation, pH and equilibrium concentrations
            // before integration (if total concentrations are valid)
            validTotalConc = solution_nmc.update(celli);

            if (validTotalConc)
            {
                // cout<<"cell id: "<<celli<<'\n';

                ode.updateUserData(celli);

                for(i=0; i<nMetals; i++)
                {
                    y[i] = totalMetalConcs[i][celli];
                }

                for(i=0; i<nMoments; i++)
                {
                    y[i + nMetals] = moments[i][celli];
                }

                try
                {
                    ode.solve(y, 0.0, odeEndTime);

                    // No data race, but risk of false sharing
                    for(i = 0; i < nMetals; i++)
                    {
                        totalMetalConcs[i][celli] = y[i];
                    }

                    // No data race, but risk of false sharing
                    for(i = 0; i < nMoments; i++)
                    {
                        moments[i][celli] = y[i + nMetals];
                    }
                }
                catch(...)
                {
                    #ifdef _OPENMP
                    #pragma omp atomic
                    #endif
                        error_flag++;
                }
            }
        }
    }

#ifdef _OPENMP
}
#endif

if (error_flag > 0)
{
    throw;
}
