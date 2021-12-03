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


/******************************************************************************
UDF to write the fluxes between different zones of a reactor network model
and the field averages inside each reactor zone to text files.
This udf should not be used with more than one compute node.
******************************************************************************/

#include "udf.h"

/* This macro is used to jump elements of the previous rows in a lower
triangular matrix */
#define JUMP_ELEMENTS(rowIndex) (((rowIndex - 1) * (rowIndex)) / 2)

#define NUM_UDM 3 /* Number of required UDMs in the current UDF */


/************************** Hard-coded definitions ***************************/
static int udm_offset = 0;

DEFINE_EXECUTE_ON_LOADING(on_loading, libname)
{
    if (N_UDM < (udm_offset + NUM_UDM))
    {
        Error("\nYou need to define up to %d extra UDMs in GUI and "
            "then reload the current library '%s'\n",
            udm_offset + NUM_UDM - N_UDM, libname);
    }
    else
    {
        Set_User_Memory_Name(udm_offset, "zoneID");
        Set_User_Memory_Name(udm_offset + 1, "zoneAverageEps");
        Set_User_Memory_Name(udm_offset + 2, "zoneAverageK");

        Message("\nThe name of two UDMs are updated\n");
    }

    Message("\n\n******************************************************************"
        "\nThe loaded library '%s' uses hard-coded variable definitions.\n"
        "Please check the definitions to be valid for the current case.\n"
        "******************************************************************\n\n", libname);
}


int n_InThread = 3; /* Total number of inlet and outlet boundaries */
int threadIn_ID[3] = {84, 82, 83}; /* IDs of the inlet and outlet boundaries assigned by Fluent*/
char inlets[3][7] = {
                    "metals",
                    "nh3",
                    "naoh"
                     };

int n_OutThread = 1; /* Total number of inlet and outlet boundaries */
int threadOut_ID[1] = {86}; /* IDs of the inlet and outlet boundaries assigned by Fluent*/


/*****************************************************************************/


DEFINE_ON_DEMAND(net_react_extr)
{
    if (compute_node_count > 1)
    {
        #if !RP_NODE /* SERIAL or HOST */
            Message("\nError:\n");
            Message("This UDF is not designed to be executed by more than one compute node.");
            Message("\nPlease re-run the Fluent with only one processor.");
        #endif
    }
    else
    {
        #if !RP_HOST /* SERIAL or NODE */

            /**********************************************************************
            Variables declaration and definition
            **********************************************************************/
            Domain *domain; /* declare domain pointer */
            Thread *f_thread;
            face_t f;
            Thread *c_thread0, *c_thread1 = NULL;
            cell_t c0, c1 = -1;

            int i, j = 0;

            int n_react_zone = 0;

            int react_zone_id0, react_zone_id1;

            int n_flux; /* Number of all possible fluxes between the pairs (i, j)
            of the reaction zones without considering connectivity */

            real *iToj_mFlux; /* A mapped 1-D array to keep the positive mass flux
            from reaction zone "i" to reaction zone "j", where "i" is the row index
            and "j" is the column index of elements of the original 2-D array. */

            real *bToi_mFlux; /* A 1-D array to keep the positive mass flux
            from boundaries to reaction zone "i". */

            int n_connectivity; /* Number of elements needed to keep the connectivity
            between the pairs (i, j) of the reaction zones */

            cxboolean *connectivity; /* A mapped 1-D array with boolean elements to
            keep the connectivity of reaction zones "i" and "j", where "i" and "j"
            are respectively the row and column indexes of a lower triangular
            matrix showing the connectivity between the reaction zone "i" and
            "j" (i > j) .*/

            cxboolean *onBoundary; /* A 1-D array with boolean elements to identify
            the zones that exchange mass with boundaries.*/

            cxboolean connected = FALSE;

            int count = 0;
            int mapped_index;

            real f_mFlux;

            real *react_zone_ave; /* An array to keep the average value
            of a selected variable in the reaction zones */
            real *react_zone_ave_k; /* An array to keep the average value of k
            in the reaction zones */
            real *react_zone_var; /* An array to keep the variance of cell values
            of a selected variable in the reaction zones */
            real *react_zone_var_k; /* An array to keep the variance of k
            in the reaction zones */
            real *react_zone_skew; /* An array to keep the skewness of cell values
            of a selected variable in the reaction zones */
            real *react_zone_vol; /* An array to keep the volume of
            the reaction zones */
            real *react_zone_rho; /* An array to keep the density of
            the reaction zones */
            real *react_zone_contErr; /* An array to keep the continuity error of
            the reaction zones */

            real cell_value = 0.0;
            real cell_value_k = 0.0;
            real cell_vol = 0.0;
            real zone_vol = 0.0;
            real globalContErr = 0.0;

            FILE *fp = NULL;
            char *filename_flux = "react_zone_flux.txt";
            char *filename_fluxB = "react_zone_feeds.txt";
            char *filename_ave = "react_zone_ave.txt";
            
            domain = Get_Domain(1); /* Get the domain using ANSYS Fluent utility */

            /*********************************************************************/

            /* Find the number of reactor zones from the maximum of zone ID*/
            thread_loop_c(c_thread0, domain) /* loops over all cell threads*/
            {
                begin_c_loop(c0, c_thread0) /* Loop over cells in a cell thread*/
                {
                    react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);

                    if (react_zone_id0 > n_react_zone)
                        n_react_zone = react_zone_id0;
                }
                end_c_loop(c0, c_thread0)
            }
            n_react_zone += 1;
            
            
            /**********************************************************************
            Calculation of fluxes between zones, fluxes from zones to boundaries
            and finding the connectivity between zones
            **********************************************************************/

            n_flux = n_react_zone * n_react_zone;

            iToj_mFlux = (real *)malloc(n_flux * sizeof(real));

            if (iToj_mFlux == NULL) {
                Message("\nError:\nmalloc of size %u for type \"real\" failed!\n",
                    n_flux);
                Message("Please execute the UDF again.\n");
                /* exit(1); */
            }

            for (i = 0; i < n_flux ; i++) /* initializing fluxes to zero */
            {
                iToj_mFlux[i] = 0.0;
            }

            bToi_mFlux = (real *)malloc(n_react_zone * n_InThread * sizeof(real));
            
            if (bToi_mFlux == NULL) {
                Message("\nError:\nmalloc of size %u for type \"real\" failed!\n",
                    n_react_zone);
                Message("Please execute the UDF again.\n");
                /* exit(1); */
            }

            for (i = 0; i < n_react_zone * n_InThread ; i++) /* initializing to zero */
            {
                bToi_mFlux[i] = 0.0;
            }

            n_connectivity = JUMP_ELEMENTS(n_react_zone);

            connectivity = (cxboolean *)malloc(n_connectivity * sizeof(cxboolean));

            if (connectivity == NULL) {
                Message("\nError:\nmalloc of size %u for type \"boolean\" failed!\n",
                    n_connectivity);
                Message("Please execute the UDF again.\n");
                /* exit(1); */
            }

            for (i = 0; i < n_connectivity ; i++) /* initializing to FALSE */
            {
                connectivity[i] = FALSE;
            }

            onBoundary = (cxboolean *)malloc(n_react_zone * sizeof(cxboolean));

            if (onBoundary == NULL) {
                Message("\nError:\nmalloc of size %u for type \"boolean\" failed!\n",
                    n_react_zone);
                Message("Please execute the UDF again.\n");
                /* exit(1); */
            }

            for (i = 0; i < n_react_zone ; i++) /* initializing to FALSE */
            {
                onBoundary[i] = FALSE;
            }

            /* loops over all face threads in a domain*/
            thread_loop_f(f_thread, domain)
            {
                c_thread0 = THREAD_T0(f_thread);
                c_thread1 = THREAD_T1(f_thread);

                if (c_thread1 != NULL)
                {
                    begin_f_loop(f, f_thread) /* loops over faces in a face thread */
                    {
                        c0 = F_C0(f, f_thread);
                        react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);

                        c1 = F_C1(f, f_thread);
                        react_zone_id1 = C_UDMI(c1, c_thread1, udm_offset);

                        if (react_zone_id0 != react_zone_id1)
                        {
                            f_mFlux = F_FLUX(f, f_thread);

                            if (f_mFlux > 0)
                            {
                                i = react_zone_id0;
                                j = react_zone_id1;

                                iToj_mFlux[i * n_react_zone + j] += f_mFlux;
                            }
                            else if (f_mFlux < 0)
                            {
                                i = react_zone_id1;
                                j = react_zone_id0;
                                
                                iToj_mFlux[i * n_react_zone + j] -= f_mFlux;
                            }

                            if (react_zone_id0 > react_zone_id1)
                            {
                                i = react_zone_id0;
                                j = react_zone_id1;
                            }
                            else
                            {
                                i = react_zone_id1;
                                j = react_zone_id0;
                            }
                            connectivity[JUMP_ELEMENTS(i) + j] = TRUE;
                        }
                    }
                    end_f_loop(f, f_thread)
                }
            }

            for (j = 0; j < n_InThread; j++)
            {
                f_thread = Lookup_Thread(domain, threadIn_ID[j]);

                if (BOUNDARY_FACE_THREAD_P(f_thread))
                {
                    c_thread0 = THREAD_T0(f_thread);

                    begin_f_loop(f, f_thread) /* loops over faces in a face thread */
                    {
                        c0 = F_C0(f, f_thread);
                        react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);

                        i = react_zone_id0;

                        f_mFlux = F_FLUX(f, f_thread);

                        if (f_mFlux < 0)
                        {
                            bToi_mFlux[i + n_react_zone * j] -= f_mFlux;
                            onBoundary[i] = TRUE;
                        }
                    }
                    end_f_loop(f, f_thread)
                }
                else
                {
                    free(iToj_mFlux);
                    free(bToi_mFlux);
                    free(connectivity);
                    free(onBoundary);
                    Error("\nThread ID %u is not a face boundary thread.", threadIn_ID[j]);
                }
            }

            for (j = 0; j < n_OutThread; j++)
            {
                f_thread = Lookup_Thread(domain, threadOut_ID[j]);

                if (BOUNDARY_FACE_THREAD_P(f_thread))
                {
                    c_thread0 = THREAD_T0(f_thread);

                    begin_f_loop(f, f_thread) /* loops over faces in a face thread */
                    {
                        c0 = F_C0(f, f_thread);
                        react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);

                        i = react_zone_id0;

                        f_mFlux = F_FLUX(f, f_thread);

                        if (f_mFlux > 0)
                        {
                            iToj_mFlux[i * n_react_zone + i] += f_mFlux;
                            onBoundary[i] = TRUE;
                        }
                    }
                    end_f_loop(f, f_thread)
                }
                else
                {
                    free(iToj_mFlux);
                    free(bToi_mFlux);
                    free(connectivity);
                    free(onBoundary);
                    Error("\nThread ID %u is not a face boundary thread.", threadOut_ID[j]);
                }
            }

            /**********************************************************************
            Calculation of zone volumes and zone averages of interested variables
            and zone continuity error using the fluxes between the zones
            **********************************************************************/

            react_zone_vol = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_ave = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_ave_k = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_var = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_var_k = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_skew = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_rho = (real *)malloc(n_react_zone * sizeof(real));
            react_zone_contErr = (real *)malloc(n_react_zone * sizeof(real));

            if (react_zone_vol == NULL || react_zone_ave == NULL ||
                react_zone_ave_k == NULL || react_zone_var_k == NULL ||
                react_zone_var == NULL || react_zone_skew == NULL ||
                react_zone_rho == NULL || react_zone_contErr == NULL)
            {
                Message("\nError:\nmalloc of size %u for type \"real\" failed!\n",
                    n_react_zone);
                Message("Please execute the UDF again.\n");
                /* exit(1); */
            }
            
            for (i = 0; i < n_react_zone ; i++) /* initializing the arrays to zero */
            {
                react_zone_vol[i] = 0.0;
                react_zone_ave[i] = 0.0;
                react_zone_ave_k[i] = 0.0;
                react_zone_var[i] = 0.0;
                react_zone_var_k[i] = 0.0;
                react_zone_skew[i] = 0.0;
                react_zone_rho[i] = 0.0;
                react_zone_contErr[i] = 0.0;
            }
            
            /* Loop over all cell threads in the domain */
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0) /* Loop over all cells */
                {
                    react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);
                    
                    cell_vol = C_VOLUME(c0, c_thread0);
                    
                    cell_value = C_D(c0, c_thread0);
                    cell_value_k = C_K(c0, c_thread0);

                    react_zone_vol[react_zone_id0] += cell_vol;
                    react_zone_ave[react_zone_id0] += cell_value * cell_vol;
                    react_zone_ave_k[react_zone_id0] += cell_value_k * cell_vol;
                    react_zone_var[react_zone_id0] += pow(cell_value, 2.0) * cell_vol;
                    react_zone_var_k[react_zone_id0] += pow(cell_value_k, 2.0) * cell_vol;
                    react_zone_skew[react_zone_id0] += pow(cell_value, 3.0) * cell_vol;
                    react_zone_rho[react_zone_id0] += C_R(c0, c_thread0) * cell_vol;
                }
                end_c_loop(c0, c_thread0)
            }
            
            

            for (i = 0; i < n_react_zone ; i++)
            {
                zone_vol = react_zone_vol[i];
                react_zone_ave[i] /= zone_vol;
                react_zone_ave_k[i] /= zone_vol;

                react_zone_var[i] /= zone_vol;
                react_zone_var_k[i] /= zone_vol;

                react_zone_skew[i] /= zone_vol;
                react_zone_skew[i] += 2*pow(react_zone_ave[i], 3.0)
                                    - 3*react_zone_ave[i]*react_zone_var[i];

                react_zone_var[i] -= pow(react_zone_ave[i], 2.0);

                react_zone_skew[i] /= pow(react_zone_var[i], 1.5);

                react_zone_rho[i] /= zone_vol;
            }

            for (i = 0; i < n_react_zone; i++)
            {
                for (j = 0; j < i; j++)
                {
                    connected = connectivity[JUMP_ELEMENTS(i) + j];

                    if (connected)
                    {
                        mapped_index = i * n_react_zone + j;

                        react_zone_contErr[i] -= iToj_mFlux[mapped_index];
                        react_zone_contErr[j] += iToj_mFlux[mapped_index];
                    }
                }

                for (j = i + 1; j < n_react_zone; j++)
                {
                    connected = connectivity[JUMP_ELEMENTS(j) + i];

                    if (connected)
                    {
                        mapped_index = i * n_react_zone + j;

                        react_zone_contErr[i] -= iToj_mFlux[mapped_index];
                        react_zone_contErr[j] += iToj_mFlux[mapped_index];
                    }
                }

                if (onBoundary[i])
                {
                    mapped_index = i * n_react_zone + i;
                    react_zone_contErr[i] -= iToj_mFlux[mapped_index];
                }
            }

            for (j = 0; j < n_InThread; j++)
            {
                for (i = 0; i < n_react_zone; i++)
                {
                    react_zone_contErr[i] += bToi_mFlux[i + n_react_zone * j];
                    
                }
            }

            for (i = 0; i < n_react_zone ; i++)
            {
                globalContErr += react_zone_contErr[i];
            }

            /* Fill the UDM with the average value of the zone. */
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {
                    react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);
                    C_UDMI(c0, c_thread0, udm_offset + 1) = react_zone_ave[react_zone_id0];
                    C_UDMI(c0, c_thread0, udm_offset + 2) = react_zone_ave_k[react_zone_id0];
                }
                end_c_loop(c0, c_thread0)
            }


            /**********************************************************************
            Write the zone information to text files
            **********************************************************************/

            if ((fp = fopen(filename_flux, "w")) == NULL)
            {
                Message("\nError: Unable to open %s for writing\n",
                    filename_flux);
            }
            else
            {
                Message("\nWriting flux between reactor zones to %s ...\n",
                    filename_flux);
                
                fprintf(fp, "%-11s %-13s %-8s %-19s %-15s\n",
                    "#", "From", "To", "Mass Flux", "Volumetric Flux");
                
                fprintf(fp, "%-11s %-13s %-8s %-19s %-6s\n\n",
                    "", "", "", "(kg/s)", "(m3/s)");

                for (i = 0; i < n_react_zone; i++)
                {
                    for (j = 0; j < i; j++)
                    {
                        connected = connectivity[JUMP_ELEMENTS(i) + j];

                        if (connected)
                        {
                            mapped_index = i * n_react_zone + j;

                            fprintf(fp, "%-8u    %-4u  %-4s    %-4u    %- 16.10f    %- 16.10f\n",
                                count, i, "--->", j, iToj_mFlux[mapped_index],
                                iToj_mFlux[mapped_index] / react_zone_rho[i]);
                            
                            count += 1;
                        }
                    }

                    if (onBoundary[i])
                    {
                        mapped_index = i * n_react_zone + i;

                        fprintf(fp, "%-8u    %-4u  %-4s    %-4u    %- 16.10f    %- 16.10f\n",
                            count, i, "--->", i, iToj_mFlux[mapped_index],
                            iToj_mFlux[mapped_index] / react_zone_rho[i]);
                            
                        count += 1;
                    }

                    for (j = i + 1; j < n_react_zone; j++)
                    {
                        connected = connectivity[JUMP_ELEMENTS(j) + i];

                        if (connected)
                        {
                            mapped_index = i * n_react_zone + j;

                            fprintf(fp, "%-8u    %-4u  %-4s    %-4u    %- 16.10f    %- 16.10f\n",
                                count, i, "--->", j, iToj_mFlux[mapped_index],
                                iToj_mFlux[mapped_index] / react_zone_rho[i]);

                            count += 1;
                        }
                    }
                }
            }
            fclose(fp);

            /*********************************************************************/

            count = 0;
            if ((fp = fopen(filename_fluxB, "w")) == NULL)
            {
                Message("\nError: Unable to open %s for writing\n",
                    filename_fluxB);
            }
            else
            {
                Message("\nWriting flux from boundaries to %s ...\n",
                    filename_fluxB);
                
                fprintf(fp, "%-11s %-13s %-8s %-19s %-19s %-7s\n",
                    "#", "From", "To", "Mass Flux", "Volumetric Flux", "Average");
                
                fprintf(fp, "%-11s %-13s %-8s %-19s %-19s %-7s\n\n",
                    "", "", "", "(kg/s)", "(m3/s)", "(m2/s3)");

                for (j = 0; j < n_InThread; j++)
                {
                    for (i = 0; i < n_react_zone; i++)
                    {
                        real bFlux = bToi_mFlux[i + n_react_zone * j];
                        if (bFlux > 0)
                        {
                            fprintf(fp, "%-8u    %-8s  %-4s    %-4u    %- 16.10f    %- 16.10f\n",
                                count, inlets[j], "--->", i, bFlux,
                                bFlux / react_zone_rho[i]);
                                
                            count += 1;
                        }
                    }
                }
            }
            fclose(fp);

            free(iToj_mFlux);
            free(bToi_mFlux);
            free(connectivity);
            free(onBoundary);

            /*********************************************************************/

            fp = NULL;

            if ((fp = fopen(filename_ave, "w")) == NULL)
            {
                Message("\nError: Unable to open %s for writing\n", filename_ave);
            }
            else
            {
                Message("\nWriting field averages of reactor zones to %s ...\n",
                    filename_ave);
                
                fprintf(fp, "%-9s %-19s %-19s %-19s %-19s %-19s %-19s %-19s %-19s %-16s\n", "Zone",
                "Average Eps", "Average K", "sdev Eps", "sdev K", "Skewness", "Volume", "Density",
                "Continuity Error", "Global Cont Error");
                    
                fprintf(fp, "%-9s %-19s %-19s %-19s %-19s %-19s %-19s %-19s %-19s %-8s\n\n", "ID",
                "(m^2/s^3)", "(m^2/s^2)","(m^2/s^3)", "(m^2/s^2)", " ","(m^3)", "(kg/m^3)",
                "(kg/s)", "(kg/s)");

                fprintf(fp, "%-5u    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f\n",
                    0, react_zone_ave[0], react_zone_ave_k[0], pow(react_zone_var[0], 0.5), pow(react_zone_var_k[0], 0.5),
                    react_zone_skew[0], react_zone_vol[0], react_zone_rho[0],
                    react_zone_contErr[0], globalContErr);

                for (i = 1; i < n_react_zone; i++)
                {
                    fprintf(fp, "%-5u    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f\n",
                    i, react_zone_ave[i], react_zone_ave_k[i], pow(react_zone_var[i], 0.5), pow(react_zone_var_k[i], 0.5), react_zone_skew[i],
                    react_zone_vol[i], react_zone_rho[i], react_zone_contErr[i]);
                }
            }
            fclose(fp);

            free(react_zone_ave);
            free(react_zone_var);
            free(react_zone_ave_k);
            free(react_zone_var_k);
            free(react_zone_vol);
            free(react_zone_skew);
            free(react_zone_rho);
            free(react_zone_contErr);

        #endif
    }
}
