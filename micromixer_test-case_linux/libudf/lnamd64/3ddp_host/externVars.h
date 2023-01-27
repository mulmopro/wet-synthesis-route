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


extern int startUDMIrp;
extern int startUDMIccr;
extern int startUDMIn;
extern int startUDMIw;
extern int startUDMIg;
extern int indexRegister;
extern int indexDissRate;
extern int indexNH3;
extern int indexOH;
extern int indexNa;
extern int indexSO4;
extern int indexP4;

extern real c_adj_h;
extern real A_p;
extern real T;

extern int env_c_rpIndex[N_UDS_C];
extern real env_conc[N_UDS_C];

extern double nucleation(double superSat);
extern double growth(double superSat, double particleSize);
extern double aggrEfficiency(double L1, double L2, double L_eq,
    double growthRate, double epsilon, double rhoLiq, double nu);
extern double aggregation(double superSat, double L1, double L2,
    double epsilon, double rhoLiq, double mu, double nu);
extern double breakage(double L1, double epsilon, double nu);
extern double erosionDD(double L1, double k);
extern double parabolicDD(double L1, double k);
extern double symmetricDD(double L1, double k);
extern double uniformDD(double L1, double k);
extern double nucleateSize(double superSat);
