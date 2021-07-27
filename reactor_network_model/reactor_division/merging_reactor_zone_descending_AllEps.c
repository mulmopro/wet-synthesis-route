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
#include <stdbool.h>

#define NUM_UDM 2

static int udm_offset = 0;

static int reactID_supersatOffset = 22;    /*It is necessary to patch the reactID_supersatOffset********************/
static int reactID_dissipRateOffset = 23;  /*and reactID_dissipRateOffset in fluent*********************************/
static int supersatAveOffset = 24;   /*The supersat and dissipRate average will be computed in this UDF*/
static int dissipRateAveOffset = 25; /******************************************************************/
static int supersatOffset = 5;


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
        Set_User_Memory_Name(udm_offset + 1, "zoneAverage");

        Message("\nThe name of two UDMs are updated\n");
    }

    Message("\n\n******************************************************************"
        "\nThe loaded library '%s' uses hard-coded variable definitions.\n"
        "Please check the definitions to be valid for the current case.\n"
        "******************************************************************\n\n", libname);
}

DEFINE_ON_DEMAND(merging_reactor_zone_descending)
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
            Thread *c_thread0 = NULL;
            cell_t c0 = -1;

            int i, j;
            real temp = 0;

            /*variables to count the numbers of comparments based on **
            **supersaturaton and dissiRate*/ 
            int n_reactDissipationRateZones = 0;
            int n_reactSupersatZones = 0;

            int react_id;
            int reactZone_id0Supersat;
            int reactZone_id0DissipRate;

            real logSupersaturation = 0.0;
            real logDissipationRate = 0.0;
            real min_logSupersaturation = 0.39794; /*The minimum value of supersaturation to use as a base to the supersaturation comparments*/

            real *reactZone_dissAve; /* An array to keep the average value
            of a dissipation rate in the reaction zones */
            real *reactZone_supersatAve; /* An array to keep the average value
            of supersaturation in the reaction zones */
            real *reactZone_volSupersat; /* An array to keep the volume of
            the supersaturation reaction zones */
            real *reactZone_volDissipRate; /* An array to keep the volume of the 
            dissipation reaction zones */
            
            real cell_vol = 0.0;
            real zone_vol = 0.0;

            domain = Get_Domain(1);

            /* Find the number of reactor zones from the maximum of zone ID*/
            thread_loop_c(c_thread0, domain) /* loops over all cell threads*/
            {
                begin_c_loop(c0, c_thread0) /* Loop over cells in a cell thread*/
                {
                    reactZone_id0Supersat    = C_UDMI(c0, c_thread0, reactID_supersatOffset);
                    reactZone_id0DissipRate = C_UDMI(c0, c_thread0, reactID_dissipRateOffset);

                    if (reactZone_id0Supersat > n_reactSupersatZones)
                        n_reactSupersatZones = reactZone_id0Supersat;

                    if (reactZone_id0DissipRate > n_reactDissipationRateZones)
                        n_reactDissipationRateZones = reactZone_id0DissipRate;
                }
                end_c_loop(c0, c_thread0)
            }
            n_reactSupersatZones++; /*Numbers of compartments based on supersaturation value*/
            n_reactDissipationRateZones++; /*Numbers of compartments based on dissipRate value*/

            Message("\nDissipRate comparments %d:\n",
                    n_reactDissipationRateZones);
            Message("\nSupersatiration  comparments %d:\n",
                    n_reactSupersatZones);


            /*************************************************************************/
            /* Calculation of zone volumes and zone averages of interested variables */
            /*************************************************************************/

            reactZone_volSupersat = (real *)malloc(n_reactSupersatZones * sizeof(real));
            reactZone_volDissipRate = (real *)malloc(n_reactDissipationRateZones * sizeof(real));
            reactZone_dissAve = (real *)malloc(n_reactDissipationRateZones * sizeof(real));
            reactZone_supersatAve = (real *)malloc(n_reactSupersatZones * sizeof(real));
            
            if (reactZone_volDissipRate == NULL || reactZone_dissAve == NULL)
            {
                Message("\nError:\nDissipationComparments: malloc of size %u for type \"real\" failed!\n",
                    n_reactDissipationRateZones);
                Message("Please execute the UDF again.\n");
            }
            if (reactZone_volSupersat == NULL || reactZone_supersatAve == NULL)
            {
                Message("\nError:\nSupersaturationCompartments: malloc of size %u for type \"real\" failed!\n",
                    n_reactSupersatZones);
                Message("Please execute the UDF again.\n");
            }

            /* initializing the arrays to zero */
            for (i = 0; i < n_reactDissipationRateZones ; i++) 
            {
                reactZone_volDissipRate[i] = 0.0;
                reactZone_dissAve[i] = 0.0;
            }       
            /* initializing the arrays to zero */
            for (i = 0; i < n_reactSupersatZones ; i++)
            {
                reactZone_volSupersat[i] = 0.0;
                reactZone_supersatAve[i] = 0.0;
            }

             /* Loop over all cell threads in the domain */
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0) /* Loop over all cells */
                {
                    reactZone_id0Supersat = C_UDMI(c0, c_thread0, reactID_supersatOffset);
                    reactZone_id0DissipRate = C_UDMI(c0, c_thread0, reactID_dissipRateOffset);
                    
                    cell_vol      = C_VOLUME(c0, c_thread0);
                    
                    logDissipationRate     = C_D(c0, c_thread0);
                    logSupersaturation = C_UDMI(c0, c_thread0, supersatOffset);
                    
                    if(logDissipationRate!=0)
                        logDissipationRate = log10(logDissipationRate);
                                        
                    if(logSupersaturation!=0)
                        logSupersaturation = log10(logSupersaturation);


                    reactZone_volSupersat[reactZone_id0Supersat] += cell_vol;
                    reactZone_volDissipRate[reactZone_id0DissipRate] += cell_vol;
                    reactZone_dissAve[reactZone_id0DissipRate] += logDissipationRate * cell_vol;
                    reactZone_supersatAve[reactZone_id0Supersat] += logSupersaturation * cell_vol;
                }
                end_c_loop(c0, c_thread0)
            }
 
            for (i = 0; i < n_reactDissipationRateZones ; i++)
            {
                zone_vol = reactZone_volDissipRate[i];
                reactZone_dissAve[i] /= zone_vol;
            }

            for (i = 0; i < n_reactSupersatZones ; i++)
            {
                zone_vol = reactZone_volSupersat[i];
                reactZone_supersatAve[i] /= zone_vol;
            }

            /* Fill the UDMs with the average value of the zone. */
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {
                    reactZone_id0Supersat = C_UDMI(c0, c_thread0, reactID_supersatOffset);
                    reactZone_id0DissipRate = C_UDMI(c0, c_thread0, reactID_dissipRateOffset);

                    C_UDMI(c0, c_thread0, supersatAveOffset) = reactZone_supersatAve[reactZone_id0Supersat];
                    C_UDMI(c0, c_thread0, dissipRateAveOffset) = reactZone_dissAve[reactZone_id0DissipRate];
                }
                end_c_loop(c0, c_thread0)
            }

            /* rearranges the supersaturations in descending order */          
            for(i = 0; i < n_reactSupersatZones ; i ++)
	        {
		
		        for(j = i + 1; j < n_reactSupersatZones; j ++)
		        {
			
			        if(reactZone_supersatAve[j] > reactZone_supersatAve[i])
		            {
			            temp = reactZone_supersatAve[i];
				        reactZone_supersatAve[i] = reactZone_supersatAve[j];
				        reactZone_supersatAve[j] = temp;
			        }
		        }
	        }

            /* rearranges the dissipation rates in descending order */          
            for(i = 0; i < n_reactDissipationRateZones; i ++)
	        {
	    
                for(j = i + 1; j < n_reactDissipationRateZones; j ++)
		        {
			        			
			        if(reactZone_dissAve[j] > reactZone_dissAve[i])
		            {
			            temp = reactZone_dissAve[i];
				        reactZone_dissAve[i] = reactZone_dissAve[j];
				        reactZone_dissAve[j] = temp;
			        }
		        }
	        }

            
            Message("Ordered average dissipRate (id | averageValue):\n");
            for(i=0; i < n_reactDissipationRateZones; i++){
                Message("\n%d",i);
                Message(" | %lf\n", reactZone_dissAve[i]);
            }

            Message("Ordered average superSat (id | averageValue):\n");
            for(i=0; i < n_reactSupersatZones; i++){
                Message("\n%d",i);
                Message(" | %lf\n", reactZone_supersatAve[i]);
            }
            

            /**************************************/
            /* Assing the ordered id to the cells */
            /**************************************/
            
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {
                    logSupersaturation = C_UDMI(c0, c_thread0, supersatAveOffset);
                    logDissipationRate = C_UDMI(c0, c_thread0, dissipRateAveOffset);

                    for(i=0 ; i < n_reactDissipationRateZones; i++){
                        if(logDissipationRate == reactZone_dissAve[i])
                            C_UDMI(c0, c_thread0, reactID_dissipRateOffset) = i;

                    }

                    for(i=0 ; i < n_reactSupersatZones; i++){
                        if(logSupersaturation == reactZone_supersatAve[i])
                            C_UDMI(c0, c_thread0, reactID_supersatOffset) = i;

                    }

                    
                }
                end_c_loop(c0, c_thread0)
            }     

            /*count the supersat comparments to be neglected*/
            int countSupersatReactNeglected = 0;

            for(i = 0; i < n_reactSupersatZones; i ++)
	        {
                logSupersaturation = reactZone_supersatAve[i];
		        if(logSupersaturation <= min_logSupersaturation)
                    countSupersatReactNeglected++;
		
		    }
            n_reactSupersatZones = n_reactSupersatZones - countSupersatReactNeglected;
            

            
            /**************************************************************/
            /* mark the supersat and dissipate pairs to be joined    ******/
            /**************************************************************/

            /*array vector to alloc dynamically a matrix of couples*/
            int *possibles_comparmentsCouples[2];
            for (i=0; i<2; i++){
                possibles_comparmentsCouples[i]=(int *)malloc(n_reactDissipationRateZones * 
                                         n_reactSupersatZones * sizeof(int));
            }        

            /*feel the matrix with -1 since no compartment can have -1 as id*/
            for (i = 0; i < (n_reactDissipationRateZones * n_reactSupersatZones) ; i++){
                possibles_comparmentsCouples[0][i] = -1;
                possibles_comparmentsCouples[1][i] = -1;

            }    
            
            /*count the new compartments born from a merge*/
            int n_MaxMergedComparments = 0;

            for(i = 0; i < n_reactDissipationRateZones; i++){
                for(j = 0; j < n_reactSupersatZones; j++){
                           possibles_comparmentsCouples[0][n_MaxMergedComparments] = i;
                           possibles_comparmentsCouples[1][n_MaxMergedComparments] = j;
                           n_MaxMergedComparments++;
                       
                }
            }

            Message("\nMaximum number of possibles couples: %d\n",
                    n_MaxMergedComparments);
             
            /* Assign the new ids obtained with the hybrid compartments logic */
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {
                    logSupersaturation = C_UDMI(c0, c_thread0, supersatAveOffset);
                    logDissipationRate = C_UDMI(c0, c_thread0, dissipRateAveOffset);

                    if(logSupersaturation <= min_logSupersaturation){
                        C_UDMI(c0, c_thread0, udm_offset) = C_UDMI(c0, c_thread0, reactID_dissipRateOffset);

                    } else {
                        reactZone_id0DissipRate = C_UDMI(c0, c_thread0, reactID_dissipRateOffset);
                        reactZone_id0Supersat = C_UDMI(c0, c_thread0, reactID_supersatOffset);

                        for(i=0; i < n_MaxMergedComparments; i++){
                            if(reactZone_id0DissipRate == possibles_comparmentsCouples[0][i] &&
                             reactZone_id0Supersat == possibles_comparmentsCouples[1][i]){                                    
                                C_UDMI(c0, c_thread0, udm_offset) = n_reactDissipationRateZones + i;

                            }

                        }

                    } 
                }
                end_c_loop(c0, c_thread0)
            }

            /****************************************************************************/
            /*  Check if some old epsilon or supersaturaton comparments                 */
            /*  does not exists anymore                                                 */
            /****************************************************************************/

            int maxFinalCompartmentsid = 0;
            thread_loop_c(c_thread0, domain) /* loops over all cell threads*/
            {
                begin_c_loop(c0, c_thread0) /* Loop over cells in a cell thread*/
                {
                    react_id    = C_UDMI(c0, c_thread0, udm_offset);
                    
                    if (react_id > maxFinalCompartmentsid)
                        maxFinalCompartmentsid = react_id;
                }
                end_c_loop(c0, c_thread0)
            }

            /***********************************************************************************/
            /*  Check if there are compartments with only one cell and merge them with  */
            /*  the previous one available                                                     */
            /***********************************************************************************/
            
            real *numberCellsInCompartments = (real *)malloc(maxFinalCompartmentsid * sizeof(real));    
            
            for (i = 0; i < maxFinalCompartmentsid ; i++){
                numberCellsInCompartments[i] = 0;
            }  
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {   
                    react_id = C_UDMI(c0, c_thread0, udm_offset);
                    numberCellsInCompartments[react_id]++; 
                    
                }end_c_loop(c0, c_thread0)
            }

            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {   
                    react_id = C_UDMI(c0, c_thread0, udm_offset);
                    if (numberCellsInCompartments[react_id] == 1){
                        for (i = react_id -1; i >= 0 ; i--){
                            if (numberCellsInCompartments[i] > 1){
                                C_UDMI(c0, c_thread0, udm_offset) = i;
                                break;
                            } else{
                                continue;
                            }
                        }
                        
                    }
                }end_c_loop(c0, c_thread0)
            }
            
            /*********************************************/
            /*  Remove the gaps in the compartments IDs  */
            /*********************************************/

            /*array to count gaps in the numeration IDs: */
            int *gapsID;
            gapsID = (int *)malloc(maxFinalCompartmentsid* sizeof(int));
            for(i=0; i<maxFinalCompartmentsid; i++){
                gapsID[i]=-1;
            }

            bool isA_Gap;
            int n_gaps=0;

            for(i=0; i<maxFinalCompartmentsid; i++){
                isA_Gap = true;
                thread_loop_c(c_thread0, domain)
                {
                    begin_c_loop(c0, c_thread0)
                    { 
                        react_id = C_UDMI(c0, c_thread0, udm_offset);
                        if(react_id == i)
                            isA_Gap = false;

                    }end_c_loop(c0, c_thread0)
                }

                if(isA_Gap == true){
                    gapsID[n_gaps]=i;
                    n_gaps++;
                }
            }
            
            Message("\nNumber of gaps: %d\n",
                    n_gaps);
            
            /*check if there are compartments with only one cell*/
            int n_major;
            thread_loop_c(c_thread0, domain)
            {
                begin_c_loop(c0, c_thread0)
                {   n_major=0;
                    react_id = C_UDMI(c0, c_thread0, udm_offset);
                    for(i=0; i<n_gaps; i++){
                        if(react_id > gapsID[i]){
                            n_major++;
                        }
                    }
                     
                    C_UDMI(c0, c_thread0, udm_offset) -= n_major;
                    
                }end_c_loop(c0, c_thread0)
            }  

            int TotalFinalCompartments = 0;
            thread_loop_c(c_thread0, domain) /* loops over all cell threads*/
            {
                begin_c_loop(c0, c_thread0) /* Loop over cells in a cell thread*/
                {
                    react_id    = C_UDMI(c0, c_thread0, udm_offset);
                    
                    if (react_id > TotalFinalCompartments)
                        TotalFinalCompartments = react_id;
                }
                end_c_loop(c0, c_thread0)
            }
            TotalFinalCompartments++;
            
            int nFinalGaps=0;
            for(i=0; i<TotalFinalCompartments; i++){
                isA_Gap = true;
                thread_loop_c(c_thread0, domain)
                {
                    begin_c_loop(c0, c_thread0)
                    { 
                        react_id = C_UDMI(c0, c_thread0, udm_offset);
                        if(react_id == i)
                            isA_Gap = false;

                    }end_c_loop(c0, c_thread0)
                }

                if(isA_Gap == true){
                    nFinalGaps++;
                }
            }


            Message("\nDissipRate compartments %d: \n",
                    n_reactDissipationRateZones);

            Message("\nTotal supersaturation compartments neglected %d: \n",
                    countSupersatReactNeglected);
            Message("\nFinal supersaturation compartments %d: \n",
                    n_reactSupersatZones);

            Message("\nTotal number of compartments: %d\n",
                    TotalFinalCompartments);

            Message("\nFinal gaps: %d\n",
                    nFinalGaps);

            free(reactZone_volSupersat);
            free(reactZone_volDissipRate);
            free(reactZone_dissAve);
            free(reactZone_supersatAve);
            free(possibles_comparmentsCouples[0]);
            free(possibles_comparmentsCouples[1]);
            free(gapsID);

        #endif
    }
}
