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
---------------------------------------------------------------------------- */


#include "udf.h"
#include "defMacros.h"


int indexP4, indexRegister, indexDissRate; /* defined in DEFINE_EXECUTE_ON_LOADING */

double a_cphi[N_CP];

DEFINE_EXECUTE_ON_LOADING(on_loading_precNMC, libname)
{
    indexP4 = N_UDS_E;
    indexRegister = indexP4 + 1;
    indexDissRate = indexRegister + 1;

    a_cphi[0] = 0.4093; a_cphi[1] = 0.6015; a_cphi[2] = 0.5851;
    a_cphi[3] = 0.09472; a_cphi[4] = -0.3903; a_cphi[5] = 0.1461;
    a_cphi[6] = -0.01604;

    int n_UDS_req = N_UDS_E;
    int n_UDM_req = indexDissRate + 1;

    if (N_UDS < n_UDS_req || N_UDM < n_UDM_req)
    {
        Message("\nThe use of the loaded library '%s' requires %d UDS and "
            "%d UDM.\nPlease make sure that enough UDS and UDM are allocated "
            "before running the simulation.\n", libname, n_UDS_req, n_UDM_req);
    }

    char envName[3]; /* Metals index 1; NaOH index 2; NH3 index 3. */
    int i;
    for (i=0; i<N_UDS_E; i++)
    {
        sprintf(envName, "P%d", i+1);
        Set_User_Scalar_Name(i, envName);
    }

    char fluxName[6];
    for (i=0; i<N_UDS_E; i++)
    {
        sprintf(fluxName, "r_p_%d", i + 1);
        Set_User_Memory_Name(i, fluxName); 
    }

    Set_User_Memory_Name(indexP4, "P4");
    Set_User_Memory_Name(indexRegister, "cell_mark");
    Set_User_Memory_Name(indexDissRate, "diss-rate-liq");

    # if !PARALLEL
        Message("\nThe name of %d UDSs and %d UDMs are updated\n",
            n_UDS_req, n_UDM_req);
    # endif

    # if RP_NODE
        if (I_AM_NODE_ZERO_P)
        {
            Message("\nThe name of %d UDSs and %d UDMs are updated\n",
                n_UDS_req, n_UDM_req);
        }
    # endif
}


DEFINE_DIFFUSIVITY(env_diffusivity, c, t, i)
{
    return TURB_VISCOSITY / SC_TURB;
}


DEFINE_ADJUST(adjust, domain)
{
    #if !RP_HOST
        Thread *t;
        cell_t c;

        int i, j, count_p;

        double env_p_array[N_UDS_E];
        
        double env_p, react_env_p_;

        double epsilon, kappa, nu, reynolds_l, log10_re_l, c_phi, gamma;

        /* looping over all cells*/ 
        thread_loop_c(t, domain)
        {
            begin_c_loop_int(c,t)
            {
                if (C_UDMI(c, t, indexRegister) > 0.0)
                {                   
                    /* Evaluation of gamma */
                    epsilon = DISS_RATE(indexDissRate);
                    kappa = TURB_KIN_ENERGY;
                    nu = MU_LIQ / C_R(c, t);

                    if (kappa > 0.0 && epsilon > 0.0)
                    {
                        reynolds_l = kappa / (pow(epsilon * nu, 0.5));
                        log10_re_l = log10(reynolds_l);
                        c_phi = 0.0;

                        if (reynolds_l > 0.2)
                        {
                            if (reynolds_l < 12853)
                            {
                                for (i=0; i<N_CP; i++)
                                {
                                    c_phi += a_cphi[i] * pow(log10_re_l, i);
                                }
                            }
                            else
                            {
                                c_phi = 2.0;
                            }
                        }
                        else
                        {
                            c_phi = 0.0;
                        }

                        gamma = MIX_CORR * c_phi * epsilon / kappa / 2.0;
                    }
                    else
                    {
                        gamma = 0.0;
                    }

                    react_env_p_ = 1.0;
                    count_p = 0;
                    for (i=0; i<N_UDS_E; i++)
                    {
                        env_p = C_UDSI(c, t, i);

                        if (env_p < 1e-20)
                        {
                            env_p_array[i] = 0.0;
                        }
                        else if (env_p > 0.9999)
                        {
                            count_p = 1;
                            react_env_p_ = 0.0;

                            for (j=0; j<N_UDS_E; j++)
                            {
                                env_p_array[j] = 0.0;
                                C_UDMI(c, t, j) = 0.0;
                            }
                            env_p_array[i] = 1.0;

                            break;
                        }
                        else
                        {
                            env_p_array[i] = env_p;
                            react_env_p_ -= env_p;
                            count_p++;
                        }
                    }

                    if (count_p > 1)
                    {
                        for (i=0; i<N_UDS_E; i++)
                        {
                            env_p = env_p_array[i];

                            if (env_p < 1e-20)
                            {
                                C_UDMI(c, t, i) = 0.0;
                            }
                            else /* && env_p < 0.9999 */
                            {
                                C_UDMI(c, t, i) = gamma * env_p * (1 - env_p);
                            }
                        }
                    }
                    else if (count_p == 1 && react_env_p_ > 1e-20)
                    {
                        for (i=0; i<N_UDS_E; i++)
                        {
                            env_p = env_p_array[i];

                            if (env_p < 1e-20)
                            {
                                C_UDMI(c, t, i) = 0.0;
                            }
                            else /* && env_p < 0.9999 */
                            {
                                C_UDMI(c, t, i) = gamma * env_p * (1 - env_p);
                            }
                        }
                    }

                    REACT_ENV_P = react_env_p_;
                }
                else /* indexRegister <= 0.0 (cell is outside the marked zone) */
                {
                    for (i=0; i<indexRegister; i++)
                    {
                        C_UDMI(c, t, i) = 0.0;
                    }
                }
            }
            end_c_loop_int(c,t)
        }
    #endif

}

