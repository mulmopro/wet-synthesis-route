double nucleation(double superSat)
{
    double J_n = 0;
    
    if (superSat > 1.0)
    {
        /* J_n = K_J*pow((superSat - 1), N_J); */
        J_n = pow(10, K_J_1) * exp(-1.0*B_J_1 / pow(log(superSat), 2))
            + pow(10, K_J_2) * exp(-1.0*B_J_2 / pow(log(superSat), 2));
        /* J_n = 1e10; */
    }

    return J_n;
}

/* Function to calculate the growth rate */
/* The function depends on supersaturation and particle size */
double growth(double superSat, double particleSize)
{
	double G_L = 0.0;

    /* if (particleSize > 1.0e-15 && superSat > 1.0)
    {
        G_L = G0*(superSat-1) / particleSize;
    } */

    if (superSat > 1.0)
    {
        G_L = K_G*pow((superSat - 1), N_G);
        /* G_L = G0; */
    }

    return G_L;
}

/* Function to calculate the aggregation efficiency */
double aggrEfficiency(double L1, double L2, double L_eq, double growthRate,
    double epsilon, double rhoLiq, double nu)
{
    double Db = pow(rhoLiq / A_p, 0.5) * pow(epsilon * nu, 0.25) * L_eq;

    double r_L = L1 / L2;

    if (L2 > L1)
    {
        r_L = L2 / L1;
    }

    double sqrt_r_L = pow(r_L*r_L - 1.0, 0.5);

    double f_lambda = 4 * (1 + r_L - sqrt_r_L)
        / ((1.0/3.0 + r_L - sqrt_r_L)
        - pow(r_L - sqrt_r_L, 2)*(2*r_L/3 + sqrt_r_L/3));

    return exp(-1.0 * pow(epsilon / nu, 0.5) * Db / f_lambda / growthRate);
}

/* Function to calculate the aggregation rate */
double aggregation(double superSat, double L1, double L2, double epsilon,
    double rhoLiq, double mu, double nu)
{
    double aggrRate = 0.0;

    double L_eq = SMALL_SIZE; /* This definition has no effect */
    double growthRate = 0.0;
    
    if (L1*L2 > 0.0)
    {
        L_eq = L1*L2 / pow(pow(L1 - L2, 2) + L1*L2, 0.5);
        growthRate = growth(superSat, L_eq);
    }

    if (L1 < 2e-3 && L2 < 2e-3 && growthRate > 0.0)
    {
        aggrRate = (2*KB*T / mu / 3)*(L1 / L2 + L2 / L1 + 2);

        aggrRate += c_adj_h*2.2943 * pow(epsilon / nu, 0.5) * pow(L1 + L2, 3);

        aggrRate *=
            aggrEfficiency(L1, L2, L_eq, growthRate, epsilon, rhoLiq, nu);
    }

    return aggrRate;
}

/* Function to calculate the breakage rate */
double breakage(double L1, double epsilon, double nu)
{
    double brRate = 0.0;
    
    if (L1 > 0.0)
    {
        brRate = C_BR_1 * pow(epsilon / nu, C_BR_2 / 2.0) * pow(L1, C_BR_3);
    }

    return brRate;
}

/* erosion daughter distribution */
double erosionDD(double L1, double k)
{
    return 1.0 + pow(pow(L1, 3) - 1.0, k/3.0);
}

/* parabolic daughter distribution */
double parabolicDD(double L1, double k)
{
    return pow(L1, k)
        * (
            3 * C_PARABOLIC_DD / (k + 3.0)
            + (1.0 - C_PARABOLIC_DD/2.0)
            * 18 * (6.0 - k) / ((k + 9.0)*(k + 6.0)*(k + 3.0))
        );
}

/* symmetric daughter distribution */
double symmetricDD(double L1, double k)
{
    return pow(2, 1.0 - k/3.0)*pow(L1, k);
}

/* uniform daughter distribution */
double uniformDD(double L1, double k)
{
    return (6.0 / (k + 3.0)) * pow(L1, k);
}

/* Function to calculate the nucleate size */
/* The function depends on supersaturation */
double nucleateSize(double superSat)
{
    /* double xc; */

    /* if (superSat > 1.0)
    {
        xc = 4*GAMMACL*VM / 2 / KB / T / log(superSat);
    }
    else
    {
        xc = 1e-10;
    } */

    return X_C;
}
