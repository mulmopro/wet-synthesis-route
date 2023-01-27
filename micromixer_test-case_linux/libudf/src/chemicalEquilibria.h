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


void equilibriaEqs(double* pConcs, const double* totalConcs, double cationTotalConc,
    double *negf, double *J, double kw, double Kb_NH3)
{
    int i, j;
    int countComplex = 0;

    double conc_OH = pow(10, -1.0*pConcs[indexOH]);
    double conc_NH3 = pow(10, -1.0*pConcs[indexNH3]);
    double conc_NH4 = Kb_NH3*conc_NH3 / conc_OH;

    double f_NH3 = 0.0;
    double J_NH3 = 0.0;
    for (i=0; i<N_METALS; i++)
    {
        int numOfComplex = numOfComplexes[i];

        double conc_i = pow(10, -1.0*pConcs[i]);

        double f_i = 0.0;
        double f_NH3_i = 0.0;
        for (j=1; j<=numOfComplex; j++)
        {
            double temp = kn_NMC[countComplex]*conc_i*pow(conc_NH3, j);
            f_i += temp;
            f_NH3_i += j*temp;
            J_NH3 += pow(j, 2)*temp;

            countComplex ++;
        }

        negf[i] = conc_i + f_i - totalConcs[i];

        J[i*N_COMPS + i] = ln10*(conc_i + f_i);

        /* The index is changed with respect to the python code
        to respect the order required by the dgesv routine of LAPACK */
        J[indexNH3*N_COMPS + i] = ln10*f_NH3_i; 
        J[i*N_COMPS + indexNH3] = ln10*f_i;

        f_NH3 += f_NH3_i;
    }

    negf[indexNH3] = conc_NH3 + conc_NH4 + f_NH3 - totalConcs[indexNH3];

    J[indexNH3*N_COMPS + indexNH3] = ln10*(conc_NH3 + conc_NH4 + J_NH3);
    J[indexOH*N_COMPS + indexNH3] = -1.0*ln10*conc_NH4;

    negf[indexOH] = conc_OH
        - (totalConcs[indexNa] - 2*totalConcs[indexSO4]) - 2*cationTotalConc - conc_NH4 - kw / conc_OH;

    J[indexOH*N_COMPS + indexOH] = ln10*(conc_NH4 + kw / conc_OH + conc_OH);
    J[indexNH3*N_COMPS + indexOH] = -1.0*ln10*conc_NH4;
}


double pureActivity(const double* z_c, const double* z_a, const double I_s,
                    const int j, const int k, double log_gamma_pure)
{
    /* First column SO4 and the second OH, rows has the same order of cations

    The matrices are defined in precNMC_adjust.c

    double B_c[6] = {0.054, 0.037, 0.049, 0.0, -0.042, 0.0875};
    double B_a[2] = {0.0, 0.076};
    double delta_c[6] = {0.21, 0.21, 0.21, 0.028, -0.02, 0.103};
    double delta_a[2] = {-0.4, -1.0};

    double B[6][2] = 
    {
        {0.1056, (B_c[0] + B_a[1] + delta_c[0]*delta_a[1])},
        {0.1226, (B_c[1] + B_a[1] + delta_c[1]*delta_a[1])},
        {0.1244, (B_c[2] + B_a[1] + delta_c[2]*delta_a[1])},
        {-0.0204, 0.0747},
        {-0.0287, (B_c[4] + B_a[1] + delta_c[4]*delta_a[1])},
        {0.0606, (B_c[5] + B_a[1] + delta_c[5]*delta_a[1])}
    }; */

    double nu_c = z_a[k];
    double nu_a = z_c[j];
    double zByZ = (nu_c*pow(z_c[j], 2) + nu_a*pow(z_a[k], 2)) / (nu_c + nu_a);

    double sqrt_I = pow(I_s, 0.5);

    log_gamma_pure = (-A_GAMMA * zByZ * sqrt_I) / (1 + sqrt_I) 
                   + ((0.06 + 0.6 * B[j][k]) * zByZ * I_s) / (pow(1.0 + 1.5 * I_s / zByZ, 2))
                   + B[j][k] * I_s - E[j][k] * ALPHA * sqrt_I * (1 - exp(-ALPHA * sqrt_I));

    return log_gamma_pure;

}


void activityBromley(const double* cationMolalConc, const double* anionMolalConc, double* gamma_ca)
{
    /* cationOrder = Ni, Mn, Co, Na, NH4, H
        anionOrder = SO4, OH */
    int i, j, k;
    double z_c[N_CATIONS] = {2, 2, 2, 1, 1, 1};
    double z_a[N_ANIONS] = {2, 1};
    double I_s = 0.0;

    for (i=0; i<N_CATIONS; i++)
    {
        I_s += cationMolalConc[i]*pow(z_c[i], 2);
    }

    for (i=0; i<N_ANIONS; i++)
    {
        I_s += anionMolalConc[i]*pow(z_a[i], 2);
    }

    I_s *= 0.5;

    for (i=0; i<N_METALS; i++)  /* Loop over the cations Ni, Mn and Co */
    {
        /* To evaluate f_c and f_a, do another loop with j and k,
           first index related to cations and second one to anions */
        double log_gamma_pure_ac = 0.0;
        double log_gamma_pure = 0.0;
        j = i;
        k = 1;
        log_gamma_pure_ac = pureActivity(z_c, z_a, I_s, j, k, log_gamma_pure);

        double f_c = 0.0;
        double f_a = 0.0;
        double sqrt_I = pow(I_s, 0.5);
        double coeff = (A_GAMMA * sqrt_I) / (1 + sqrt_I);
        
        /* Calculate f_c */
        for (k=0; k<N_ANIONS; k++)
        {
            j = i;
            double y_ac = 0.0;

            if (k != 1)
            {
                log_gamma_pure = pureActivity(z_c, z_a, I_s, j, k, log_gamma_pure);

                y_ac = pow((z_c[j] + z_a[k]), 2) * anionMolalConc[k] / (4 * I_s);
                f_c += y_ac * (log_gamma_pure + coeff *z_c[j] * z_a[k]);
            }
            else
            {
                y_ac = pow((z_c[j] + z_a[k]), 2) * anionMolalConc[k] / (4 * I_s);
                f_c += y_ac * (log_gamma_pure_ac + coeff *z_c[j] * z_a[k]);
            }
        }

        /* Calculate f_a */
        for (j=0; j<N_CATIONS; j++)
        {
            k = 1;
            double y_ac = 0.0;

            if (j != i)
            {
                log_gamma_pure = pureActivity(z_c, z_a, I_s, j, k, log_gamma_pure);

                y_ac = pow((z_c[j] + z_a[k]), 2) * cationMolalConc[j] / (4 * I_s);
                f_a += y_ac * (log_gamma_pure + coeff *z_c[j] * z_a[k]);
            }
            else
            {
                y_ac = pow((z_c[j] + z_a[k]), 2) * cationMolalConc[j] / (4 * I_s);
                f_a += y_ac * (log_gamma_pure_ac + coeff *z_c[j] * z_a[k]);
            }
        }

        double nu_c = z_a[1];
        double nu_a = z_c[i];

        double log_gamma_ca = (-coeff * (nu_c*pow(z_c[i], 2) + nu_a*pow(z_a[1], 2))
                            + (nu_c*f_c + nu_a*f_a)) / (nu_c + nu_a);
        
        gamma_ca[i] = pow(10, log_gamma_ca);        
    }
}


void solveEquilibria(const double* totalConcs, double* pConcs, double cationTotalConc,
    const double* cationConcRatios, double* equilConcs, double rhoLiq, double *pH, double* superSat,
    int* interruptFlag)
{
    double negf[N_COMPS];
    double J[N_COMPS*N_COMPS];
    double pKw, kw, pKb_NH3, Kb_NH3;    

    pKw = Aw + Bw/T + Cw/(pow(T, 2)) + Dw/(pow(T, 3)) + (Ew + Fw/T + Gw/(pow(T, 3)))*log10(rhoLiq/1000);
    kw = POW10(pKw);
    pKb_NH3 = -pKw - (0.09018+2729.92/T);
    Kb_NH3 = POW10(-pKb_NH3);

    int i;
    for (i=0; i<N_COMPS; i++)
    {
        negf[i] = 0.0;
    }

    for (i=0; i<N_COMPS*N_COMPS; i++)
    {
        J[i] = 0.0;
    }

    int j;
    int n = N_COMPS, info;
    int ipiv[N_COMPS];
    double error;
    for(i=0; i<MAX_ITER; i++)
    {
        equilibriaEqs(pConcs, totalConcs, cationTotalConc, negf, J, kw, Kb_NH3);
        
        dgesv(&n, &nrhs, J, &lda, ipiv, negf, &ldb, &info);
        
        if (info == 0)
        {
            error = 0.0;
            for (j=0; j<N_COMPS; j++)
            {
                error += fabs(negf[j]);
            }
        
            if (error < TOLERANCE)
            {
                break;
            }
            else
            {
                for (j=0; j<N_COMPS; j++)
                {
                    pConcs[j] += negf[j];
                }
            }
        }
        else if (info > 0) /* Check for the exact singularity */
        {
            printf("The diagonal element of the triangular factor of A,\n");
            printf("U(%i,%i) is zero, so that A is singular;\n", info, info);
            printf("the solution could not be computed.\n");
            *interruptFlag = info;
            break;
        }
        else
        {
            exit(1);
        }
        
    }

    if (info == 0)
    {
        for (j=0; j<N_COMPS; j++)
        {
            equilConcs[j] = pow(10, -1.0*pConcs[j]);
        }

        double conc_OH = equilConcs[indexOH];
        double conc_NH4 = Kb_NH3*equilConcs[indexNH3] / conc_OH;
        double conc_H = kw / conc_OH;
        double gamma_ca[N_METALS];
        double cationMolalConc[N_CATIONS];
        double anionMolalConc[N_ANIONS];

        for (j=0; j<N_METALS; j++)
        {
            cationMolalConc[j] = totalConcs[j];
        }
        
        cationMolalConc[3] = totalConcs[indexNa];
        cationMolalConc[4] = conc_NH4;
        cationMolalConc[5] = conc_H;

        anionMolalConc[0] = totalConcs[indexSO4];
        anionMolalConc[1] = conc_OH;

        activityBromley(cationMolalConc, anionMolalConc, gamma_ca);

        double k_sp_NMC = 1;
        double powConcs_NMC = 1;
        double cationConcRatio;
        for (j=0; j<N_METALS; j++)
        {
            cationConcRatio = cationConcRatios[j];
            k_sp_NMC *= pow(k_sp[j], cationConcRatio);

            powConcs_NMC *= pow(equilConcs[j]*pow(gamma_ca[j], 3), cationConcRatio);
        }

        *pH = 14 - pConcs[indexOH];

        *superSat = pow(powConcs_NMC*conc_OH*conc_OH / k_sp_NMC, 1.0/3.0);
    }
}


