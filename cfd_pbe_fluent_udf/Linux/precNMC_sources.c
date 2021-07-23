#include "udf.h"
#include "defMacros.h"
#include "externVars.h"

DEFINE_SOURCE(conc_source, c, t, dS, eqn)
{
    #if !RP_HOST
        int udsIndex = eqn - EQ_UDS;

        double source = 0.0;

        if (C_UDMI(c, t, indexRegister) > 0.0)
        {
            if (udsIndex < N_METALS)
            {
                source = -1.0*PREC_RATE*C_UDMI(c, t, udsIndex + startUDMIccr);
            }

            source += C_UDMI(c, t, startUDMIrp + env_c_rpIndex[udsIndex])
                * env_conc[udsIndex];

        }

        source *= C_R(c, t);

        dS[eqn] = 0.0;

        return source;
    #endif
}


DEFINE_SOURCE(env_source, c, t, dS, eqn)
{
    #if !RP_HOST
        int envIndex = eqn - EQ_UDS - N_UDS_C;

        double source = 0.0;

        /* printf("\n Gamma: %f \n", GAMMA); */
        
        if (C_UDMI(c, t, indexRegister) > 0.0)
        {
            source = -C_UDMI(c, t, envIndex + startUDMIrp)*C_R(c, t);
        }

        dS[eqn] = 0.0;

        return source;
    #endif
}


DEFINE_SOURCE(mom_source, c, t, dS, eqn)
{
    #if !RP_HOST
        int i;
        int momIndex = eqn - EQ_UDS - N_UDS_C - N_UDS_E;
        double superSat = SUPERSATURATION;
        double rhoLiq = C_R(c, t);
        double source = 0.0;
        double source_a = 0.0;
        double source_b = 0.0;

        if (C_UDMI(c, t, indexRegister) > 0.0 && superSat > 1.0)
        {
            if (momIndex == 3)
            {
                for(i=0; i<N_NODES; i++)
                {
                    /* moment order will be multiplied later */
                    source += C_UDMI(c, t, i + startUDMIw)
                        * C_UDMI(c, t, i + startUDMIg)
                        * pow(C_UDMI(c, t, i + startUDMIn), 2);
                }
            }
            else
            {
                int j;
                double epsilon = DISS_RATE(indexDissRate);
                double mu = MU_LIQ;
                double nu = mu / rhoLiq;
                double L_i, L_j, w_i, L_iToPow3, L_iToPowK;

                for(i=0; i<N_NODES; i++)
                {
                    w_i = C_UDMI(c, t, i + startUDMIw);
                    L_i = C_UDMI(c, t, i + startUDMIn);
                    L_iToPow3 = pow(L_i, 3);
                    L_iToPowK = pow(L_i, momIndex);

                    /* moment order will be multiplied later */
                    source += w_i * C_UDMI(c, t, i + startUDMIg)
                        * pow(L_i, momIndex - 1);

                    source_b += w_i * breakage(L_i, epsilon, nu)
                        * (BR_DAUGHTER_DIST - L_iToPowK);
                    
                    for(j=0; j<N_NODES; j++)
                    {
                        L_j = C_UDMI(c, t, j + startUDMIn);

                        source_a += w_i * C_UDMI(c, t, j + startUDMIw)
                            * aggregation(superSat, L_i, L_j, epsilon, rhoLiq, mu, nu)
                            * (
                                0.5*pow(L_iToPow3 + pow(L_j, 3), momIndex / 3.0)
                                - L_iToPowK
                            );
                    }
                }
            }

            source *= momIndex;

            /* Adding the nucleation, aggregation and breakage contributions */
            source += NUC_RATE*pow(NUCLEATE_SIZE, momIndex)
                + source_a + source_b;

            source *= rhoLiq * REACT_ENV_P;
        }

        dS[eqn] = 0.0;

        return source;
    #endif
}
