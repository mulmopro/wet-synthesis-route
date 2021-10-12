/* This file generated automatically. */
/*          Do not modify.            */
#include "udf.h"
#include "prop.h"
#include "dpm.h"
extern DEFINE_EXECUTE_ON_LOADING(on_loading_precNMC, libname);
extern DEFINE_DIFFUSIVITY(conc_diffusivity, c, t, i);
extern DEFINE_DIFFUSIVITY(env_diffusivity, c, t, i);
extern DEFINE_ADJUST(adjust, domain);
extern DEFINE_SOURCE(conc_source, c, t, dS, eqn);
extern DEFINE_SOURCE(env_source, c, t, dS, eqn);
extern DEFINE_SOURCE(mom_source, c, t, dS, eqn);
UDF_Data udf_data[] = {
{"on_loading_precNMC", (void (*)(void))on_loading_precNMC, UDF_TYPE_EXECUTE_ON_LOADING},
{"conc_diffusivity", (void (*)(void))conc_diffusivity, UDF_TYPE_DIFFUSIVITY},
{"env_diffusivity", (void (*)(void))env_diffusivity, UDF_TYPE_DIFFUSIVITY},
{"adjust", (void (*)(void))adjust, UDF_TYPE_ADJUST},
{"conc_source", (void (*)(void))conc_source, UDF_TYPE_SOURCE},
{"env_source", (void (*)(void))env_source, UDF_TYPE_SOURCE},
{"mom_source", (void (*)(void))mom_source, UDF_TYPE_SOURCE},
};
int n_udf_data = sizeof(udf_data)/sizeof(UDF_Data);
#include "version.h"
void UDF_Inquire_Release(int *major, int *minor, int *revision)
{
  *major = RampantReleaseMajor;
  *minor = RampantReleaseMinor;
  *revision = RampantReleaseRevision;
}
