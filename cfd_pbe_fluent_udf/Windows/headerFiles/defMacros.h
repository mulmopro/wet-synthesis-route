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


#define SMALL_R 1e-38
#define SMALL_CONC 1e-14
#define EFFECTIVE_CONC 1e-9
#define SMALL_SIZE 1e-20
#define SMALL_M3 2e-20

#define N_NODES 2    /* number of quadrature nodes */
#define N_UDS_E 3    /* number of user-defined scalars (for environments) */ 
#define N_UDS_C 6     /* number of user-defined scalars (for concentrations) */
#define N_COMPS 5    /* number of equilibrium concentrations*/
#define N_METALS 3
#define N_CATIONS 6
#define N_ANIONS 2
#define N_COMPLEXES_NI 6
#define N_COMPLEXES_MN 4
#define N_COMPLEXES_CO 6
#define N_COMPLEXES 16
#define MAX_ITER 200
#define TOLERANCE 1e-6

/* Temperature */
#define T 298.15 /* Kelvin */

/* Turbulent viscosity and turbulent Schmidt number */
#define TURB_VISCOSITY C_MU_T(c, t)
#define SC_TURB 1.0

/* Aggregation model parameters */
/* #define C_ADJ_H 1  Correction coefficient for the hydrodynamic aggregation */
/* #define A_P 1e6  yield stress of crystals */
#define DISS_RATE(i) C_UDMI(c, t, i)
#define MU_LIQ C_MU_L(c, t) /* 1e-3 */
#define TURB_KIN_ENERGY C_K(c, t) /* 0.01011 for 2D simulations */

/* Breakage model parameters */
#define C_BR_1 0.0
#define C_BR_2 0
#define C_BR_3 0
#define BR_DAUGHTER_DIST uniformDD(L_i, momIndex) /* erosionDD, parabolicDD, symmetricDD, uniformDD */
#define C_PARABOLIC_DD 4

/* #define G0 1e-4 constant growth rate */
#define K_G 1e-7 /* growth rate parameter */
#define N_G 1 /* growth rate parameter */

#define K_J_1 30 /* nucleation rate parameter */
#define B_J_1 2000 /* nucleation rate parameter */
#define K_J_2 12 /* nucleation rate parameter */
#define B_J_2 20 /* nucleation rate parameter */
#define X_C 1e-9 /* nucleate size */

/* crystal properties */
#define RHO_CRYST 3953 /* density of crystals in kg/m3 */
#define MW_CRYST 92.3383 /* molecular weight of crystal [kg/kmol] */

/* Bromley's activity coefficients constants */
#define A_GAMMA 0.511 /* valid for T = 25 °C */
#define ALPHA 70.0

/* Micromixing parameters */ 
#define MIX_CORR 2.85
#define N_CP 7 /* Number of correlation parameters */

#define KV 0.523599 /* volume shape factor */
#define KB 1.38064852e-23  /* Boltzmann number m^2 kg s^-2 K^-1 */

#define POW10(a) pow(10, a)

#define SUPERSATURATION C_UDMI(c, t, 5)
#define PH C_UDMI(c, t, 6)
#define NUC_RATE C_UDMI(c, t, 7)
#define NUCLEATE_SIZE C_UDMI(c, t, 8)
#define SMD C_UDMI(c, t, 9)
#define PREC_RATE C_UDMI(c, t, 10)
#define REACT_ENV_P C_UDMI(c, t, 23)
