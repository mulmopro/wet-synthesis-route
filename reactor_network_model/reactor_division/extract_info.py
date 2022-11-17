import inspect
import os
import sys
import traceback
# import math
import time  # argparse, importlib, inspect, traceback
# import pickle
import numpy as np

from scipy.sparse import dok_matrix
from sklearn.cluster import AgglomerativeClustering


# Entry to the main function
def main(argv):

    # Find location of the executed script
    frame = inspect.currentframe()
    filePath = inspect.getfile(frame)
    caseDir = os.path.realpath(os.path.abspath(os.path.dirname(filePath)))
    # Path to the import files
    importDir = os.path.join(caseDir, 'importScripts')

    # Start the calculations
    try:
        from importScripts.read_files import (read_field, read_owner_nbr,
            read_boundary)

        OF_data_path = os.path.join(caseDir, '../cfdSimulation')

        time_dir = '0'

        density = 998.2

        inOutPatches = {
            'inlets':{
                0:{'patchName': 'inlet-metal', 'name': 'metals'},
                1:{'patchName': 'inlet-nh3', 'name': 'nh3'},
                2:{'patchName': 'inlet-naoh', 'name': 'naoh'}},
            'outlets':{
                0:{'patchName': 'outlet', 'name': 'outlet'}}}

        fieldName = 'zone_id'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        zone_id, _ = read_field(field_path, read_boundary=False)

        fieldName = 'phi'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        phi, phi_b = read_field(field_path)

        fieldName = 'epsilon'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        epsilon, _ = read_field(field_path, read_boundary=False)

        fieldName = 'k'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        kappa, _ = read_field(field_path, read_boundary=False)

        fieldName = 'V'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        volume, _ = read_field(field_path, read_boundary=False)

        mesh_path = os.path.join(OF_data_path, 'constant', 'polyMesh')

        owner_path = os.path.join(mesh_path, 'owner')
        owner_cells = read_owner_nbr(owner_path)

        neighbor_path = os.path.join(mesh_path, 'neighbour')
        neighbor_cells = read_owner_nbr(neighbor_path)

        boundary_path = os.path.join(mesh_path, 'boundary')
        boundary_faces = read_boundary(boundary_path)

        num_cells = max(np.amax(neighbor_cells), np.amax(owner_cells)) + 1

        cells = np.linspace(
            0, num_cells, num=num_cells, endpoint=False, dtype='int32')

        num_zones = int(np.amax(zone_id) + 1)

        num_inlets = int(3)

        iToj_flux = dok_matrix((num_zones, num_zones), dtype=np.float64)
        bToi_flux = dok_matrix((num_inlets, num_zones), dtype=np.float64)

        ''' calculate fluxes between the zones '''
        for face, (neighbor, owner) in enumerate(
                zip(neighbor_cells, owner_cells)):

            owner_zone_id = zone_id[owner]
            nbr_zone_id = zone_id[neighbor]

            if (owner_zone_id != nbr_zone_id):
                face_flux = phi[face]

                if face_flux > 0.0:  # from owner to neighbor
                    iToj_flux[owner_zone_id, nbr_zone_id] += face_flux
                elif face_flux < 0.0:  # from neighbor to owner
                    iToj_flux[nbr_zone_id, owner_zone_id] -= face_flux

        ''' calculate fluxes from the zones to the outlets and
            from the inlets to the zones.
            The outgoing fluxes are saved as diagonal elements! '''
        for key in inOutPatches:
            patches = inOutPatches[key]

            for index in patches.keys():
                patch = patches[index]
                patchName = patch['patchName']

                nFaces = boundary_faces[patchName]['nFaces']
                startFace = boundary_faces[patchName]['startFace']

                phi_patch = phi_b[patchName]['value']

                for face in range(nFaces):

                    owner = owner_cells[startFace + face]

                    owner_zone = zone_id[owner]

                    face_flux = phi_patch[face]

                    if key == 'outlets':
                        if face_flux > 0.0:  # from owner to outlet
                            iToj_flux[owner_zone, owner_zone] += face_flux
                    elif key == 'inlets':
                        if face_flux < 0.0:  # from owner to inlet
                            bToi_flux[index, owner_zone] -= face_flux

        ''' calculate zone volumes and field averages '''
        zone_vol = np.zeros(num_zones)
        zone_epsilon_volAve = np.zeros(num_zones)
        zone_k_volAve = np.zeros(num_zones)

        zone_numCells = np.zeros(num_zones)
        zone_epsilon_numAve = np.zeros(num_zones)
        zone_k_numAve = np.zeros(num_zones)
        zone_epsilon_numVar = np.zeros(num_zones)
        zone_k_numVar = np.zeros(num_zones)

        # Volume average of epsilon and kappa are calculated in zones
        # For statistics of the division, variances are also calculated
        # without volume averaging
        for cell in range(num_cells):
            cell_zone_id = int(zone_id[cell])
            cell_vol = volume[cell]
            cell_epsilon = epsilon[cell]
            cell_k = kappa[cell]

            zone_vol[cell_zone_id] += cell_vol

            zone_epsilon_volAve[cell_zone_id] += cell_vol * cell_epsilon
            zone_k_volAve[cell_zone_id] += cell_vol * cell_k

            zone_numCells[cell_zone_id] += 1

            zone_epsilon_numAve[cell_zone_id] += cell_epsilon
            zone_epsilon_numVar[cell_zone_id] += cell_epsilon**2

            zone_k_numAve[cell_zone_id] += cell_k
            zone_k_numVar[cell_zone_id] += cell_k**2

        # complete the calculation of the averages, variances and skewness
        # the order of calcualtions matters!
        for i in range(num_zones):
            zone_i_vol = zone_vol[i]

            zone_epsilon_volAve[i] /= zone_i_vol
            zone_k_volAve[i] /= zone_i_vol

            zone_i_numCells = zone_numCells[i]

            zone_epsilon_numAve[i] /= zone_i_numCells  # 1st moment

            zone_epsilon_numVar[i] /= zone_i_numCells  # 2nd moment
            zone_epsilon_numVar[i] -= zone_epsilon_numAve[i]**2.0

            zone_k_numAve[i] /= zone_i_numCells  # 1st moment

            zone_k_numVar[i] /= zone_i_numCells  # 2nd moment
            zone_k_numVar[i] -= zone_k_numAve[i]**2.0

        ''' sort sparse matrices'''
        iToj_flux = iToj_flux.tocsr()
        iToj_flux.sort_indices()
        iToj_flux = iToj_flux.tocoo()

        bToi_flux = bToi_flux.tocsr()
        bToi_flux.sort_indices()
        bToi_flux = bToi_flux.tocoo()

        zone_contErr = np.zeros(num_zones)
        for (i, j, flux) in zip(iToj_flux.row, iToj_flux.col, iToj_flux.data):

            if i == j:
                zone_contErr[i] -= flux
            else:
                zone_contErr[i] -= flux
                zone_contErr[j] += flux

        for (j, flux) in zip(bToi_flux.col, bToi_flux.data):
            zone_contErr[j] += flux

        globalContErr = 0.0
        for i in range(num_zones):
            globalContErr += zone_contErr[i]

        globalContErr *= density

        ''' write zone fluxes '''
        fileName = 'react_zone_flux.txt'
        filePath = os.path.join(caseDir, fileName)

        print("\nWriting fluxes between reactor zones to \"{}\" ...\n".format(
            fileName))
        with open(filePath, 'w') as f:

            f.write("{:<11s} {:<13s} {:<7s} {:<22s} {:<15s}\n".format(
                "#", "From", "To", "Mass Flux", "Volumetric Flux"))

            f.write("{:<11s} {:<13s} {:<7s} {:<22s} {:<6s}\n\n".format(
                "", "", "", "(kg/s)", "(m3/s)"))

            for count, (i, j, flux) in enumerate(
                    zip(iToj_flux.row, iToj_flux.col, iToj_flux.data)):

                f.write(
                    "{:<8d}    {:<4d}  {:<4s}    {:<4d}    {:<.13e}    "
                    "{:<.13e}\n".format(
                        count + 1, i, "--->", j, flux*998.2, flux))

        ''' write inlet boundary fluxes '''
        fileName = 'react_zone_feeds.txt'
        filePath = os.path.join(caseDir, fileName)

        print("\nWriting inlet fluxes from boundaries to \"{}\" ...\n".format(
            fileName))
        with open(filePath, 'w') as f:

            f.write("{:<11s} {:<17s} {:<7s} {:<22s} {:<15s}\n".format(
                "#", "From", "To", "Mass Flux", "Volumetric Flux"))

            f.write("{:<11s} {:<17s} {:<7s} {:<22s} {:<6s}\n\n".format(
                "", "", "", "(kg/s)", "(m3/s)"))

            for count, (i, j, flux) in enumerate(
                    zip(bToi_flux.row, bToi_flux.col, bToi_flux.data)):

                inlet_name = inOutPatches['inlets'][i]['name']

                f.write(
                    "{:<8d}    {:<8s}  {:<4s}    {:<4d}    {:<.13e}    "
                    "{:<.13e}\n".format(
                        count + 1, inlet_name, "--->", j, flux*density, flux))

        ''' write zone volumes and averages '''
        fileName = 'react_zone_ave.txt'
        filePath = os.path.join(caseDir, fileName)

        print("\nWriting field averages of zones to \"{}\" ...\n".format(
            fileName))
        with open(filePath, 'w') as f:

            f.write(
                "{:<10s} {:<22s} {:<22s} {:<22s} "
                "{:<17s} {:<19s} {:<16s}\n".format(
                    "Zone ID", "Volume Average eps", "Volume Average K",
                    "Volume", "Density", "Continuity Error",
                    "Global Cont Error"))

            f.write(
                "{:<10s} {:<22s} {:<22s} {:<22s} "
                "{:<17s} {:<19s} {:<6s}\n".format(
                    "#", "(m^2/s^3)", "(m^2/s^2)",
                     "(m^3)", "(kg/m^3)", "(kg/s)", "(kg/s)"))

            for i in range(num_zones):

                line = \
                    "{:<7d}    {:<.13e}    {:<.13e}    "\
                    "{:<.13e}    {:<.8e}    {:<+.8e}".format(
                        i, zone_epsilon_volAve[i], zone_k_volAve[i],
                        zone_vol[i], density, zone_contErr[i]*density)

                if i == 0:
                    line = line + '    {:<+.8e}\n'.format(globalContErr)
                else:
                    line = line + '\n'

                f.write(line)

        import warnings

        with warnings.catch_warnings():
            warnings.simplefilter("ignore")

            pass


    # If the necessary files to be imported are not found in the
    # .\caseDir\importScripts
    except ImportError as e:
        print(str(e))
        print('\nThe above python script is not found in the following'
              + 'address:\n' + importDir)
    # All other exceptions raised during calculation will be processed here
    except Exception as e:

            print(str(e))
            traceback.print_exc()
            sys.exit(1)
    finally:
        print('\n')


# Python starts the program from here because it is the only executable code
# at indentation level 0. Python checks if the code is executed as a script
# and not imported by other scripts
if __name__ == "__main__":
    # If it is executed as script, it calls function "main" with arguments
    # entered by the user
    main(sys.argv[1:])
