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
#include "externFuncs.h"

/* Declaration of external variables */
int startUDMIccr; /* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIn; /* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIw; /* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIg; /* defined in DEFINE_EXECUTE_ON_LOADING */
int startUDMIrp; /* defined in DEFINE_EXECUTE_ON_LOADING */

int numOfComplexes[N_METALS] = {
    N_COMPLEXES_NI, N_COMPLEXES_MN, N_COMPLEXES_CO
};

double kn_NMC[N_COMPLEXES]; /* defined in DEFINE_EXECUTE_ON_LOADING */
double k_sp[N_METALS]; /* defined in DEFINE_EXECUTE_ON_LOADING */
double Kb_NH3; /* defined in DEFINE_EXECUTE_ON_LOADING */
double kw = 1e-14;
double B[N_CATIONS][N_ANIONS]; /* defined in DEFINE_EXECUTE_ON_LOADING */
double E[N_CATIONS][N_ANIONS]; /* defined in DEFINE_EXECUTE_ON_LOADING */

double D_molecular[N_UDS_C];

double a_cphi[N_CP];

/* Aggregation model parameters to be set by the scheme variables */
real c_adj_h;
real A_p;

real env_conc[N_UDS_C];
int env_c_rpIndex[N_UDS_C];

double ln10;

int indexNH3 = N_METALS;
int indexOH; /* defined in DEFINE_EXECUTE_ON_LOADING */
int indexNa; /* defined in DEFINE_EXECUTE_ON_LOADING */
int indexSO4; /* defined in DEFINE_EXECUTE_ON_LOADING */

int indexP4, indexRegister, indexDissRate; /* defined in DEFINE_EXECUTE_ON_LOADING */

int nrhs = 1, lda = N_COMPS, ldb = N_COMPS;

#include "momentCalc.h"
#include "chemicalEquilibria.h"
#include "particleProcesses.h"


DEFINE_EXECUTE_ON_LOADING(on_loading_precNMC, libname)
{
    startUDMIccr = 11;
    startUDMIn = startUDMIccr + N_METALS;
    startUDMIw = startUDMIn + N_NODES;
    startUDMIg = startUDMIw + N_NODES;

    indexOH = N_COMPS - 1;
    indexNa = N_UDS_C - 2;
    indexSO4 = N_UDS_C - 1;

    startUDMIrp = startUDMIg + N_NODES;
    indexP4 = startUDMIrp + N_UDS_E;
    indexRegister = indexP4 + 1;
    indexDissRate = indexRegister + 1;

    ln10 = log(10.0);

    kn_NMC[0] = POW10(2.81); kn_NMC[1] = POW10(5.08); kn_NMC[2] = POW10(6.85);
    kn_NMC[3] = POW10(8.12); kn_NMC[4] = POW10(8.93); kn_NMC[5] = POW10(9.08);

    kn_NMC[6] = POW10(1.00); kn_NMC[7] = POW10(1.54); kn_NMC[8] = POW10(1.70);
    kn_NMC[9] = POW10(1.30);

    kn_NMC[10] = POW10(2.10); kn_NMC[11] = POW10(3.67); kn_NMC[12] = POW10(4.78);
    kn_NMC[13] = POW10(5.53); kn_NMC[14] = POW10(5.75); kn_NMC[15] = POW10(5.14);

    k_sp[0] = POW10(-15.22); k_sp[1] = POW10(-12.70); k_sp[2] = POW10(-14.89);

    Kb_NH3 = POW10(-4.8);

    B[0][0] = 0.1056; B[0][1] = -0.080;
    B[1][0] = 0.1226; B[1][1] = -0.097;
    B[2][0] = 0.1244; B[2][1] = -0.085;
    B[3][0] = -0.0204; B[3][1] = 0.0747;
    B[4][0] = -0.0287; B[4][1] = 0.0540;
    B[5][0] = 0.0606; B[5][1] = 0.0605;

    E[0][0] = 0.00524; E[0][1] = 0.0;
    E[1][0] = 0.00599; E[1][1] = 0.0;
    E[2][0] = 0.00498; E[2][1] = 0.0;
    E[3][0] = 0.0; E[3][1] = 0.0;
    E[4][0] = 0.0; E[4][1] = 0.0;
    E[5][0] = 0.0; E[5][1] = 0.0;

    D_molecular[0] = 2e-9; D_molecular[1] = 2e-9; D_molecular[2] = 2e-9;
    D_molecular[3] = 2e-9; D_molecular[4] = 2e-9; D_molecular[5] = 2e-9;

    a_cphi[0] = 0.4093; a_cphi[1] = 0.6015; a_cphi[2] = 0.5851;
    a_cphi[3] = 0.09472; a_cphi[4] = -0.3903; a_cphi[5] = 0.1461;
    a_cphi[6] = -0.01604;

    env_c_rpIndex[0] = 0; env_c_rpIndex[1] = 0; env_c_rpIndex[2] = 0;
    env_c_rpIndex[3] = 1; env_c_rpIndex[4] = 2; env_c_rpIndex[5] = 0;

    int n_UDS_req = N_UDS_C + N_UDS_E + 2*N_NODES;
    int n_UDM_req = indexDissRate + 1;

    if (N_UDS < n_UDS_req || N_UDM < n_UDM_req)
    {
        Message("\nThe use of the loaded library '%s' requires %d UDS and "
            "%d UDM.\nPlease make sure that enough UDS and UDM are allocated "
            "before running the simulation.\n", libname, n_UDS_req, n_UDM_req);
    }

    Set_User_Scalar_Name(0, "totC_Ni");
    Set_User_Scalar_Name(1, "totC_Mn");
    Set_User_Scalar_Name(2, "totC_Co");
    Set_User_Scalar_Name(3, "totC_NH3");
    Set_User_Scalar_Name(4, "totC_Na");
    Set_User_Scalar_Name(5, "totC_SO4");

    char envName[3]; /* Metals index 1; NaOH index 2; NH3 index 3. */
    int i;
    for (i=0; i<N_UDS_E; i++)
    {
        sprintf(envName, "P%d", i+1);
        Set_User_Scalar_Name(N_UDS_C + i, envName);
    }

    char momentName[3];
    for (i=0; i<2*N_NODES; i++)
    {
        sprintf(momentName, "M%d", i);
        Set_User_Scalar_Name(N_UDS_C + N_UDS_E + i, momentName);
    }

    Set_User_Memory_Name(0, "eqC_Ni");
    Set_User_Memory_Name(1, "eqC_Mn");
    Set_User_Memory_Name(2, "eqC_Co");
    Set_User_Memory_Name(3, "eqC_NH3");
    Set_User_Memory_Name(4, "eqC_OH");
    Set_User_Memory_Name(5, "superSat");
    Set_User_Memory_Name(6, "pH");
    Set_User_Memory_Name(7, "nucRate");
    Set_User_Memory_Name(8, "nuclSize");
    Set_User_Memory_Name(9, "SMD");
    Set_User_Memory_Name(10, "precRate");
    Set_User_Memory_Name(11, "cRatio_Ni");
    Set_User_Memory_Name(12, "cRatio_Mn");
    Set_User_Memory_Name(13, "cRatio_Co");

    char nodeName[3];
    char growthName[5];
    for (i=0; i<N_NODES; i++)
    {
        sprintf(nodeName, "L%d", i + 1);
        Set_User_Memory_Name(startUDMIn + i, nodeName);
        sprintf(nodeName, "w%d", i + 1);
        Set_User_Memory_Name(startUDMIw + i, nodeName);
        sprintf(growthName, "G_L%d", i + 1);
        Set_User_Memory_Name(startUDMIg + i, growthName);
    }

    char fluxName[6];
    for (i=0; i<N_UDS_E; i++)
    {
        sprintf(fluxName, "r_p_%d", i + 1);
        Set_User_Memory_Name(startUDMIrp + i, fluxName); 
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


DEFINE_DIFFUSIVITY(conc_diffusivity, c, t, i)
{
    return C_R(c, t)*D_molecular[i] + TURB_VISCOSITY / SC_TURB;
}

DEFINE_DIFFUSIVITY(env_diffusivity, c, t, i)
{
    return TURB_VISCOSITY / SC_TURB;
}


DEFINE_ADJUST(adjust, domain)
{
    int interruptFlag = 0;

    # if !RP_NODE
        if (RP_Variable_Exists_P("aggregation/c-t") &&
            RP_Variable_Exists_P("aggregation/a-p") &&
            RP_Variable_Exists_P("env_conc/ni") &&
            RP_Variable_Exists_P("env_conc/mn") &&
            RP_Variable_Exists_P("env_conc/co") &&
            RP_Variable_Exists_P("env_conc/nh3") &&
            RP_Variable_Exists_P("env_conc/na") &&
            RP_Variable_Exists_P("env_conc/so4"))
        {
            c_adj_h = RP_Get_Real("aggregation/c-t");
            A_p = RP_Get_Real("aggregation/a-p");
            env_conc[0] = RP_Get_Real("env_conc/ni");
            env_conc[1] = RP_Get_Real("env_conc/mn");
            env_conc[2] = RP_Get_Real("env_conc/co");
            env_conc[3] = RP_Get_Real("env_conc/nh3");
            env_conc[4] = RP_Get_Real("env_conc/na");
            env_conc[5] = RP_Get_Real("env_conc/so4");
        }
        else
        {
            Error("\nScheme variables are not defined.\n");
        }
    # endif

    host_to_node_real_2(c_adj_h, A_p);
    host_to_node_real(env_conc, N_UDS_C);

    #if !RP_HOST
        Thread *t;
        cell_t c;

        int i;
        int numOfMoments = 2*N_NODES;

        double totalConcs[N_UDS_C];
        double pConcs[N_COMPS], equilConcs[N_COMPS];
        double cationConcRatios[N_METALS];

        double env_p[N_UDS_E];

        double nodes[N_NODES], weights[N_NODES], growthRates[N_NODES];
        double moments[2*N_NODES];

        double equilConc, totalConc, cationTotalConc, cationConcRatio, conc_Na;
        double pH, superSat, nuclRate, nuclSize, dm3dt;
        double conc_OH, powConcs_NMC, k_sp_NMC;
        double epsilon, kappa, nu, reynolds_l, log10_re_l, c_phi, gamma;

        cxboolean validConc, validMoment;


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
                    nu = MU_LIQ / C_R (c, t);

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
                    
                    REACT_ENV_P = 1.0;
                    for (i=0; i<N_UDS_E; i++)
                    {
                        env_p[i] = C_UDSI(c, t , i + N_UDS_C);

                        REACT_ENV_P -= env_p[i];

                        if (env_p[i] > 0.0 && env_p[i] < 1.0)
                        {
                            C_UDMI(c, t, i + startUDMIrp) = gamma * env_p[i] * (1 - env_p[i]);
                        }
                        else
                        {
                            C_UDMI(c, t, i + startUDMIrp) = 0.0;
                        }
                    }

                    if (interruptFlag == 0 && REACT_ENV_P > 0.0001)
                    {

                        for (i=0; i<N_UDS_C; i++)
                        {
                            totalConcs[i] = C_UDSI(c, t, i) / REACT_ENV_P;
                        }

                        validConc = TRUE;
                        for (i=0; i<N_METALS; i++)
                        {
                            if (totalConcs[i] < EFFECTIVE_CONC)
                            {
                                validConc = FALSE;
                                break;
                            }
                        }
                        
                        if (validConc)
                        {
                            if (totalConcs[indexNH3] > EFFECTIVE_CONC)
                            {
                                cationTotalConc = 0.0;
                                for (i=0; i<N_METALS; i++)
                                {
                                    equilConc = C_UDMI(c, t, i);
                                    totalConc = totalConcs[i];
                                    if (equilConc > 0.0 && equilConc < totalConc)
                                    {
                                        pConcs[i] = -1.0*log10(equilConc);
                                    }
                                    else
                                    {
                                        pConcs[i] = -1.0*log10(totalConc);
                                    }
                                    cationTotalConc += totalConcs[i];
                                }

                                equilConc = C_UDMI(c, t, indexOH);
                                conc_Na = totalConcs[indexNa];
                                if (equilConc > 1e-7)
                                {
                                    pConcs[indexOH] = -1.0*log10(equilConc);
                                }
                                else if (conc_Na > 0.001)
                                {
                                    pConcs[indexOH] = -1.0*log10(conc_Na);
                                }
                                else
                                {
                                    pConcs[indexOH] = -1.0*log10(0.001);
                                }

                                equilConc = C_UDMI(c, t, indexNH3);
                                totalConc = totalConcs[indexNH3];
                                if (equilConc > EFFECTIVE_CONC && equilConc < totalConc)
                                {
                                    pConcs[indexNH3] = -1.0*log10(equilConc);
                                }
                                else
                                {
                                    pConcs[indexNH3] = -1.0*log10(totalConc);
                                }

                                for (i=0; i<N_METALS; i++)
                                {
                                    cationConcRatio = totalConcs[i] / cationTotalConc;
                                    cationConcRatios[i] = cationConcRatio;
                                    C_UDMI(c, t, i + startUDMIccr) = cationConcRatio;
                                }

                                solveEquilibria(totalConcs, pConcs, cationTotalConc,
                                    cationConcRatios, equilConcs, &pH, &superSat, &interruptFlag);

                                if (interruptFlag == 0)
                                {
                                    for (i=0; i<N_COMPS; i++)
                                    {
                                        C_UDMI(c, t, i) = equilConcs[i];
                                    }

                                    SUPERSATURATION = superSat;
                                    PH = pH;

                                    nuclRate = nucleation(superSat);
                                }
                                else
                                {
                                    superSat = 0.0;
                                    SUPERSATURATION = 0.0;
                                    PH = 7.0;

                                    nuclRate = 0.0;
                                }
                            }
                            else /* totalConcs[indexNH3] <= EFFECTIVE_CONC */
                            {
                                cationTotalConc = 0.0;
                                for (i=0; i<N_METALS; i++)
                                {
                                    cationTotalConc += totalConcs[i];
                                }

                                for (i=0; i<N_METALS; i++)
                                {
                                    cationConcRatio = totalConcs[i] / cationTotalConc;
                                    cationConcRatios[i] = cationConcRatio;
                                    C_UDMI(c, t, i + startUDMIccr) = cationConcRatio;
                                }

                                conc_OH = totalConcs[indexNa] - 2*totalConcs[indexSO4] + 2*cationTotalConc;

                                if (conc_OH > 1e-7)
                                {
                                    k_sp_NMC = 1.0;
                                    powConcs_NMC = 1.0;
                                    for (i=0; i<N_METALS; i++)
                                    {
                                        cationConcRatio = cationConcRatios[i];
                                        k_sp_NMC *= pow(k_sp[i], cationConcRatio);
                                        powConcs_NMC *= pow(totalConcs[i], cationConcRatio);
                                    }

                                    superSat = pow(powConcs_NMC*conc_OH*conc_OH / k_sp_NMC, 1.0/3.0);

                                    C_UDMI(c, t, indexOH) = conc_OH;
                                    SUPERSATURATION = superSat;
                                    PH = 14 + log10(conc_OH);

                                    nuclRate = nucleation(superSat);
                                }
                                else
                                {
                                    superSat = 0.0;
                                    C_UDMI(c, t, indexOH) = 1e-7;
                                    SUPERSATURATION = 0.0;
                                    PH = 7.0;

                                    nuclRate = 0.0;
                                }
                            }
                        }
                        else /* validConc = False */
                        {
                            superSat = 0.0;
                            /* SUPERSATURATION = 0.0; */

                            nuclRate = 0.0;

                            for (i=0; i<startUDMIrp; i++)
                            {
                                C_UDMI(c, t, i) = 0.0;
                            }
                        }

                        NUC_RATE = nuclRate;

                        nuclSize = nucleateSize(superSat);
                        NUCLEATE_SIZE = nuclSize;

                        validMoment = TRUE;
                        for (i=0; i<numOfMoments; i++)
                        {
                            moments[i] = C_UDSI(c, t, i + N_UDS_C + N_UDS_E) / REACT_ENV_P;

                            if (moments[i] <= 0.0)
                            {
                                validMoment = FALSE;
                            }
                        }

                        if (validMoment && moments[3] > SMALL_M3)
                        {
                            adaptiveWheeler(nodes, weights, moments);

                            if (moments[2] > 0.0)
                            {
                                SMD = moments[3] / moments[2];
                            }
                            else
                            {
                                SMD = SMALL_SIZE;
                            }

                            for (i=0; i<N_NODES; i++)
                            {
                                C_UDMI(c, t, i + startUDMIn) = nodes[i];
                                C_UDMI(c, t, i + startUDMIw) = weights[i];

                                growthRates[i] = growth(superSat, nodes[i]);
                                C_UDMI(c, t, i + startUDMIg) = growthRates[i];
                            }

                            dm3dt = 0.0;
                            for(i=0; i<N_NODES; i++)
                            {
                                dm3dt += growthRates[i]*weights[i]*pow(nodes[i], 2);
                            }
                            dm3dt *= 3.0;
                            dm3dt += nuclRate*pow(nuclSize, 3);

                            PREC_RATE = (KV*RHO_CRYST / MW_CRYST)*dm3dt*REACT_ENV_P;
                        }
                        else
                        {                       
                            SMD = SMALL_SIZE;

                            for (i=0; i<N_NODES; i++)
                            {
                                C_UDMI(c, t, i + startUDMIn) = SMALL_SIZE;
                                C_UDMI(c, t, i + startUDMIw) = 0.0;

                                C_UDMI(c, t, i + startUDMIg) = 0.0;
                            }

                            dm3dt = nuclRate*pow(nuclSize, 3);

                            PREC_RATE = (KV*RHO_CRYST / MW_CRYST)*dm3dt*REACT_ENV_P;
                        }
                    }
                    else /* interruptFlag != 0 && REACT_ENV_P <= 0.0 */
                    {
                        for (i=0; i<startUDMIrp; i++)
                        {
                            C_UDMI(c, t, i) = 0.0;
                        }
                    }
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

    # if RP_NODE /* Does nothing in Serial */
        interruptFlag = PRF_GISUM1(interruptFlag);
    # endif

    node_to_host_int_1(interruptFlag);

    # if !RP_NODE
        if (interruptFlag > 0) {
            RP_Set_Integer("interrupt/flag", 1);
        }
    # endif
}
