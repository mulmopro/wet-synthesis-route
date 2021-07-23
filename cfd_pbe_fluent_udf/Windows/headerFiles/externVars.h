/* Header file that contains the external variable/function declaration */

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
