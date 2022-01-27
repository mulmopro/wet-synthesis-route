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
    #if !RP_NODE  /* SERIAL or HOST */
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
    #endif  /* !RP_NODE  */
}


int n_InThread = 3; /* Total number of inlet and outlet boundaries */
int threadIn_ID[3] = {15, 17, 16}; /* IDs of the inlet and outlet boundaries assigned by Fluent*/
char inlets[3][7] = {
                    "metals",
                    "nh3",
                    "naoh"
                     };

int n_OutThread = 1; /* Total number of inlet and outlet boundaries */
int threadOut_ID[1] = {24}; /* IDs of the inlet and outlet boundaries assigned by Fluent*/


/*****************************************************************************/


DEFINE_ON_DEMAND(net_react_extr)
{

    int n_react_zone = 0;

    /*************************************************************************/

    #if !RP_HOST  /* SERIAL or NODE */

        Domain *domain; /* declare domain pointer */
        Thread *c_thread0 = NULL;
        cell_t c0 = -1;

        int react_zone_id0;

        domain = Get_Domain(1); /* Get the domain using ANSYS Fluent utility */

        /* Find the number of reactor zones from the maximum of zone ID*/
        thread_loop_c(c_thread0, domain) /* loops over all cell threads*/
        {
            begin_c_loop_int(c0, c_thread0) /* Loop over cells in a cell thread*/
            {
                react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);

                if (react_zone_id0 > n_react_zone)
                    n_react_zone = react_zone_id0;
            }
            end_c_loop_int(c0, c_thread0)
        }
        n_react_zone += 1;

        n_react_zone = PRF_GIHIGH1(n_react_zone);  /* global maximum */

    #endif  /* !RP_HOST */

    node_to_host_int_1(n_react_zone);  /* copy number of zones to host */

    /*************************************************************************/

    int n_flux; /* Number of all possible fluxes between the pairs (i, j)
        of the reaction zones without considering connectivity */

    n_flux = n_react_zone * n_react_zone;

    double *iToj_mFlux; /* A mapped 1-D array to keep the positive mass flux
        from reaction zone "i" to reaction zone "j", where "i" is the row index
        and "j" is the column index of elements of the original 2-D array. */

    iToj_mFlux = (double *)malloc(n_flux * sizeof(double));

    if (iToj_mFlux == NULL) {
        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_flux);
    }

    for (int i = 0; i < n_flux ; i++) /* initializing fluxes to zero */
    {
        iToj_mFlux[i] = 0.0;
    }

    double *bToi_mFlux; /* A 1-D array to keep the positive mass flux
        from boundaries to reaction zone "i". */

    bToi_mFlux = (double *)malloc(n_react_zone * n_InThread * sizeof(double));

    if (bToi_mFlux == NULL) {
        free(iToj_mFlux);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone * n_InThread);
    }

    for (int i = 0; i < n_react_zone * n_InThread ; i++) /* initializing to zero */
    {
        bToi_mFlux[i] = 0.0;
    }

    int n_connectivity; /* Number of elements needed to keep the connectivity
        between the pairs (i, j) of the reaction zones */

    n_connectivity = JUMP_ELEMENTS(n_react_zone);

    cxboolean *connectivity; /* A mapped 1-D array with boolean elements to
        keep the connectivity of reaction zones "i" and "j", where "i" and "j"
        are respectively the row and column indexes of a lower triangular
        matrix showing the connectivity between the reaction zone "i" and
        "j" (i > j) .*/

    connectivity = (cxboolean *)malloc(n_connectivity * sizeof(cxboolean));

    if (connectivity == NULL) {
        free(iToj_mFlux); free(bToi_mFlux);

        Error("\nError:\nmalloc of size %u for type \"boolean\" failed!\n",
            "Please execute the UDF again.\n", n_connectivity);
    }

    for (int i = 0; i < n_connectivity ; i++) /* initializing to FALSE */
    {
        connectivity[i] = FALSE;
    }

    cxboolean *onBoundary; /* A 1-D array with boolean elements to identify
        the zones that exchange mass with boundaries.*/

    onBoundary = (cxboolean *)malloc(n_react_zone * sizeof(cxboolean));

    if (onBoundary == NULL) {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity);

        Error("\nError:\nmalloc of size %u for type \"boolean\" failed!\n",
            "Please execute the UDF again.\n", n_react_zone);
    }

    for (int i = 0; i < n_react_zone ; i++) /* initializing to FALSE */
    {
        onBoundary[i] = FALSE;
    }

    /*************************************************************************
        Calculation of fluxes between zones, fluxes from zones to boundaries
        and finding the connectivity between zones
    **************************************************************************/

    #if !RP_HOST  /* SERIAL or NODE */

        Thread *f_thread;
        Thread *c_thread1 = NULL;

        face_t f;
        cell_t c1 = -1;

        int react_zone_id1;

        double f_mFlux;

        /* loops over all face threads in a domain*/
        thread_loop_f(f_thread, domain)
        {
            c_thread0 = THREAD_T0(f_thread);
            c_thread1 = THREAD_T1(f_thread);

            if (c_thread1 != NULL)
            {
                int i, j;

                begin_f_loop(f, f_thread) /* loops over faces in a face thread */
                {
                    if (PRINCIPAL_FACE_P(f, f_thread))
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
                }
                end_f_loop(f, f_thread)
            }
        }

        for (int j = 0; j < n_InThread; j++)
        {
            f_thread = Lookup_Thread(domain, threadIn_ID[j]);

            if (BOUNDARY_FACE_THREAD_P(f_thread))
            {
                c_thread0 = THREAD_T0(f_thread);

                int i;

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
                free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
                Error("\nThread ID %u is not a face boundary thread.", threadIn_ID[j]);
            }
        }

        for (int j = 0; j < n_OutThread; j++)
        {
            f_thread = Lookup_Thread(domain, threadOut_ID[j]);

            int i;

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
                free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
                Error("\nThread ID %u is not a face boundary thread.", threadOut_ID[j]);
            }
        }

        double *iwork_r;
        iwork_r = (double *)malloc(n_flux * sizeof(double));
        if (iwork_r == NULL) {
            Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
                "Please execute the UDF again.\n", n_flux);
        }

        PRF_GDSUM(iToj_mFlux, n_flux, iwork_r);  /* global sum for fluxes */
        free(iwork_r);

        iwork_r = (double *)malloc(n_react_zone * n_InThread * sizeof(double));
        if (iwork_r == NULL) {
            Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
                "Please execute the UDF again.\n", n_react_zone * n_InThread);
        }

        PRF_GDSUM(bToi_mFlux, n_react_zone * n_InThread, iwork_r);  /* global sum for fluxes */
        free(iwork_r);

        cxboolean *iwork_b;
        iwork_b = (cxboolean *)malloc(n_connectivity * sizeof(cxboolean));
        if (iwork_b == NULL) {
            Error("\nError:\nmalloc of size %u for type \"boolean\" failed!\n",
                "Please execute the UDF again.\n", n_connectivity);
        }

        PRF_GBOR(connectivity, n_connectivity, iwork_b);  /* global sum for fluxes */
        free(iwork_b);

        iwork_b = (cxboolean *)malloc(n_react_zone * sizeof(cxboolean));
        if (iwork_b == NULL) {
            Error("\nError:\nmalloc of size %u for type \"boolean\" failed!\n",
                "Please execute the UDF again.\n", n_react_zone);
        }

        PRF_GBOR(onBoundary, n_react_zone, iwork_b);  /* global sum for fluxes */
        free(iwork_b);

    #endif  /* !RP_HOST */

    /* copy arrays to host */
    node_to_host_double(iToj_mFlux, n_flux);  
    node_to_host_double(bToi_mFlux, n_react_zone * n_InThread);

    node_to_host_boolean(connectivity, n_connectivity);
    node_to_host_boolean(onBoundary, n_react_zone);

    /*************************************************************************/

    double *react_zone_vol; /* An array to keep the volume of
    the reaction zones */
    react_zone_vol = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_vol == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_ave; /* An array to keep the average value
    of a selected variable in the reaction zones */
    react_zone_ave = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_ave == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_ave_k; /* An array to keep the average value of k
    in the reaction zones */
    react_zone_ave_k = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_ave_k == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol); free(react_zone_ave);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_var; /* An array to keep the variance of cell values
    of a selected variable in the reaction zones */
    react_zone_var = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_var == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol); free(react_zone_ave); free(react_zone_ave_k);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_var_k; /* An array to keep the variance of k
    in the reaction zones */
    react_zone_var_k = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_var_k == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol); free(react_zone_ave); free(react_zone_ave_k);
        free(react_zone_var);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_skew; /* An array to keep the skewness of cell values
    of a selected variable in the reaction zones */
    react_zone_skew = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_skew == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol); free(react_zone_ave); free(react_zone_ave_k);
        free(react_zone_var); free(react_zone_var_k);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_rho; /* An array to keep the density of
    the reaction zones */
    react_zone_rho = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_rho == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol); free(react_zone_ave); free(react_zone_ave_k);
        free(react_zone_var); free(react_zone_var_k); free(react_zone_skew);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    double *react_zone_contErr; /* An array to keep the continuity error of
    the reaction zones */
    react_zone_contErr = (double *)malloc(n_react_zone * sizeof(double));
    if (react_zone_contErr == NULL)
    {
        free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);
        free(react_zone_vol); free(react_zone_ave); free(react_zone_ave_k);
        free(react_zone_var); free(react_zone_var_k); free(react_zone_skew);
        free(react_zone_rho);

        Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
            "Please execute the UDF again.\n", n_react_zone);
    }

    for (int i = 0; i < n_react_zone ; i++) /* initializing the arrays to zero */
    {
        react_zone_vol[i] = 0.0; react_zone_ave[i] = 0.0; react_zone_ave_k[i] = 0.0;
        react_zone_var[i] = 0.0; react_zone_var_k[i] = 0.0; react_zone_skew[i] = 0.0;
        react_zone_rho[i] = 0.0; react_zone_contErr[i] = 0.0;
    }

    #if !RP_HOST  /* SERIAL or NODE */

        double cell_vol = 0.0;
        double cell_value = 0.0;
        double cell_value_k = 0.0;

        /* Loop over all cell threads in the domain */
        thread_loop_c(c_thread0, domain)
        {
            begin_c_loop_int(c0, c_thread0) /* Loop over all cells */
            {
                react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);

                cell_vol = C_VOLUME(c0, c_thread0);

                cell_value = C_D(c0, c_thread0);
                cell_value_k = C_K(c0, c_thread0);

                react_zone_vol[react_zone_id0] += cell_vol;

                react_zone_ave[react_zone_id0] += cell_value * cell_vol;
                react_zone_var[react_zone_id0] += pow(cell_value, 2.0) * cell_vol;
                react_zone_skew[react_zone_id0] += pow(cell_value, 3.0) * cell_vol;

                react_zone_ave_k[react_zone_id0] += cell_value_k * cell_vol;
                react_zone_var_k[react_zone_id0] += pow(cell_value_k, 2.0) * cell_vol;

                react_zone_rho[react_zone_id0] += C_R(c0, c_thread0) * cell_vol;
            }
            end_c_loop_int(c0, c_thread0)
        }

        iwork_r = (double *)malloc(n_react_zone * sizeof(double));
        if (iwork_r == NULL) {
            Error("\nError:\nmalloc of size %u for type \"double\" failed!\n"
                "Please execute the UDF again.\n", n_react_zone);
        }

        PRF_GDSUM(react_zone_vol, n_react_zone, iwork_r);  /* global sum */
        PRF_GDSUM(react_zone_ave, n_react_zone, iwork_r);  /* global sum */
        PRF_GDSUM(react_zone_ave_k, n_react_zone, iwork_r);  /* global sum */
        PRF_GDSUM(react_zone_var, n_react_zone, iwork_r);  /* global sum */
        PRF_GDSUM(react_zone_var_k, n_react_zone, iwork_r);  /* global sum */
        PRF_GDSUM(react_zone_skew, n_react_zone, iwork_r);  /* global sum */
        PRF_GDSUM(react_zone_rho, n_react_zone, iwork_r);  /* global sum */
        free(iwork_r);

    #endif  /* !RP_HOST  */

    node_to_host_double(react_zone_vol, n_react_zone);
    node_to_host_double(react_zone_ave, n_react_zone);
    node_to_host_double(react_zone_ave_k, n_react_zone);
    node_to_host_double(react_zone_var, n_react_zone);
    node_to_host_double(react_zone_var_k, n_react_zone);
    node_to_host_double(react_zone_skew, n_react_zone);
    node_to_host_double(react_zone_rho, n_react_zone);

    /*************************************************************************
        Calculation of zone volumes and zone averages of interested variables
        and zone continuity error using the fluxes between the zones
    **************************************************************************/

    #if !RP_NODE  /* SERIAL or HOST */

        double zone_vol = 0.0;

        for (int i = 0; i < n_react_zone ; i++)
        {
            zone_vol = react_zone_vol[i];

            react_zone_ave[i] /= zone_vol;  /* 1st moment */

            react_zone_var[i] /= zone_vol;  /* 2nd moment */

            react_zone_skew[i] /= zone_vol;  /* 3rd moment */
            react_zone_skew[i] += 2*pow(react_zone_ave[i], 3.0)
                                - 3*react_zone_ave[i]*react_zone_var[i];

            /* This final update of variance should come after using
            the second moment for the calculation of skewness (see above lines)*/
            react_zone_var[i] -= pow(react_zone_ave[i], 2.0);

            /* This final update of skewness should come after updating variance */
            react_zone_skew[i] /= pow(react_zone_var[i], 1.5);

            react_zone_ave_k[i] /= zone_vol;  /* 1st moment */

            react_zone_var_k[i] /= zone_vol;  /* 2nd moment */
            react_zone_var_k[i] -= pow(react_zone_ave_k[i], 2.0);

            react_zone_rho[i] /= zone_vol;
        }

    #endif  /* !RP_NODE  */

    host_to_node_double(react_zone_ave, n_react_zone);
    host_to_node_double(react_zone_ave_k, n_react_zone);

    #if !RP_HOST  /* SERIAL or NODE */

        /* Fill the UDM with the average value of the zone. */
        thread_loop_c(c_thread0, domain)
        {
            begin_c_loop_int(c0, c_thread0)
            {
                react_zone_id0 = C_UDMI(c0, c_thread0, udm_offset);
                C_UDMI(c0, c_thread0, udm_offset + 1) = react_zone_ave[react_zone_id0];
                C_UDMI(c0, c_thread0, udm_offset + 2) = react_zone_ave_k[react_zone_id0];
            }
            end_c_loop_int(c0, c_thread0)
        }

    #endif  /* !RP_HOST  */

    #if !RP_NODE  /* SERIAL or HOST */

        int mapped_index;
        cxboolean connected = FALSE;

        for (int i = 0; i < n_react_zone; i++)
        {
            for (int j = 0; j < i; j++)
            {
                connected = connectivity[JUMP_ELEMENTS(i) + j];

                if (connected)
                {
                    mapped_index = i * n_react_zone + j;

                    react_zone_contErr[i] -= iToj_mFlux[mapped_index];
                    react_zone_contErr[j] += iToj_mFlux[mapped_index];
                }
            }

            for (int j = i + 1; j < n_react_zone; j++)
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

        for (int j = 0; j < n_InThread; j++)
        {
            for (int i = 0; i < n_react_zone; i++)
            {
                react_zone_contErr[i] += bToi_mFlux[i + n_react_zone * j];
                
            }
        }

        double globalContErr = 0.0;

        for (int i = 0; i < n_react_zone ; i++)
        {
            globalContErr += react_zone_contErr[i];
        }

        /**********************************************************************
            Write the zone information to text files
        **********************************************************************/

        int count = 0;

        FILE *fp = NULL;
        char *filename_flux = "react_zone_flux.txt";
        char *filename_fluxB = "react_zone_feeds.txt";
        char *filename_ave = "react_zone_ave.txt";

        if ((fp = fopen(filename_flux, "w")) == NULL)
        {
            Message("\nError: Unable to open %s for writing\n", filename_flux);
        }
        else
        {
            Message("\nWriting flux between reactor zones to %s ...\n",
                filename_flux);
            
            fprintf(fp, "%-11s %-13s %-8s %-19s %-15s\n",
                "#", "From", "To", "Mass Flux", "Volumetric Flux");
            
            fprintf(fp, "%-11s %-13s %-8s %-19s %-6s\n\n",
                "", "", "", "(kg/s)", "(m3/s)");

            for (int i = 0; i < n_react_zone; i++)
            {
                for (int j = 0; j < i; j++)
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

                for (int j = i + 1; j < n_react_zone; j++)
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
            Message("\nError: Unable to open %s for writing\n", filename_fluxB);
        }
        else
        {
            Message("\nWriting flux from boundaries to %s ...\n",
                filename_fluxB);
            
            fprintf(fp, "%-11s %-13s %-8s %-19s %-19s %-7s\n",
                "#", "From", "To", "Mass Flux", "Volumetric Flux", "Average");
            
            fprintf(fp, "%-11s %-13s %-8s %-19s %-19s %-7s\n\n",
                "", "", "", "(kg/s)", "(m3/s)", "(m2/s3)");

            for (int j = 0; j < n_InThread; j++)
            {
                for (int i = 0; i < n_react_zone; i++)
                {
                    double bFlux = bToi_mFlux[i + n_react_zone * j];
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

            for (int i = 1; i < n_react_zone; i++)
            {
                fprintf(fp, "%-5u    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f    %- 16.10f\n",
                i, react_zone_ave[i], react_zone_ave_k[i], pow(react_zone_var[i], 0.5), pow(react_zone_var_k[i], 0.5), react_zone_skew[i],
                react_zone_vol[i], react_zone_rho[i], react_zone_contErr[i]);
            }
        }
        fclose(fp);

    #endif  /* !RP_NODE  */

    free(iToj_mFlux); free(bToi_mFlux); free(connectivity); free(onBoundary);

    free(react_zone_vol); free(react_zone_ave); free(react_zone_ave_k);
    free(react_zone_var); free(react_zone_var_k); free(react_zone_skew);
    free(react_zone_rho); free(react_zone_contErr);
}
