import inspect
import os
import sys
import traceback
# import math
import time  # argparse, importlib, inspect, traceback
# import pickle
import re
import struct
import numpy as np

from scipy.sparse import csr_matrix
from sklearn.cluster import AgglomerativeClustering
from sklearn.preprocessing import StandardScaler


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
        from importScripts.read_files import read_field, read_owner_nbr

        OF_data_path = os.path.join(caseDir, '../cfdSimulation')

        time_dir = '0'

        fieldName = 'p_env1'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        p_env1, _ = read_field(field_path, read_boundary=False)

        fieldName = 'p_env3'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        p_env3, _ = read_field(field_path, read_boundary=False)

        p_Metal_NaOH = p_env1 + p_env3

        fieldName = 'epsilon'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)
        epsilon, _ = read_field(field_path, read_boundary=False)

        mesh_path = os.path.join(OF_data_path, 'constant', 'polyMesh')

        owner_path = os.path.join(mesh_path, 'owner')
        owner_cells = read_owner_nbr(owner_path)

        neighbor_path = os.path.join(mesh_path, 'neighbour')
        neighbor_cells = read_owner_nbr(neighbor_path)

        num_cells = max(np.amax(neighbor_cells), np.amax(owner_cells)) + 1

        cells = np.linspace(
            0, num_cells, num=num_cells, endpoint=False, dtype='int32')

        connectivity = construct_connectivity(
            cells, neighbor_cells, owner_cells)

        ''' specification data '''

        low_limit_p = 0.0001
        high_limit_p = 0.9999

        redo_regions = True

        low_limit_epsilon = 1e-5

        num_zones_out = 22

        num_zones_in = 24

        ''' Mark the region of interest '''

        target_condition = ~np.logical_or(
            p_Metal_NaOH < low_limit_p, p_Metal_NaOH > high_limit_p)

        target_values = np.zeros(num_cells, dtype=np.float16)
        target_values[target_condition] = np.float16(1.0)

        num_target_inlets = 2
        num_regions = 2*num_target_inlets + 1

        regionLabels_path = os.path.join(
            caseDir,
            "regionLabels_L{}_H{}.npy".format(low_limit_p, high_limit_p))

        if redo_regions:
            print("\nDividing the domain into", num_regions, "regions ...\n")
            model_region = AgglomerativeClustering(
                n_clusters=num_regions, connectivity=connectivity,
                linkage="ward", compute_full_tree=False)

            model_region.fit(target_values.reshape(-1, 1))

            regionLabels = model_region.labels_

            np.save(regionLabels_path, regionLabels)

        elif os.path.isfile(regionLabels_path):
            print("\nLoad regions from the disk ...\n")
            regionLabels = np.load(regionLabels_path)

        else:
            print("\nRegion labels not found! ...\n")
            exit()

        regions_in = []
        regions_in_numCells = []
        regions_out = []
        regions_out_meanEpsilon = []
        for i in range(num_regions):
            region_condition = (regionLabels == i)
            region_cells = region_condition.nonzero()[0]

            region_p = p_Metal_NaOH[region_condition]
            region_p_mean = np.mean(region_p)

            region_epsilon = epsilon[region_condition]

            if (region_p_mean < low_limit_p or region_p_mean > high_limit_p):
                regions_out.append(i)
                regions_out_meanEpsilon.append(np.mean(region_epsilon))

            else:
                regions_in.append(i)
                regions_in_numCells.append(region_cells.size)

        sort_index = list(range(len(regions_out_meanEpsilon)))
        sort_index.sort(key=regions_out_meanEpsilon.__getitem__)
        regions_out = [regions_out[i] for i in sort_index]

        sort_index = list(range(len(regions_in_numCells)))
        sort_index.sort(key=regions_in_numCells.__getitem__, reverse=True)
        regions_in = [regions_in[i] for i in sort_index]

        # Array to save the ID of the zones
        zone_id_field = np.zeros(num_cells, dtype=np.int32)

        ''' Dividing the region far from the inlets '''

        count = np.int32(0)
        for i in regions_out:
            region_condition = (regionLabels == i)
            region_cells = region_condition.nonzero()[0]

            connectivity_sliced = connectivity[region_cells, :]
            region_connectivity = connectivity_sliced[:, region_cells]

            region_epsilon = epsilon[region_condition]
            region_epsilon[region_epsilon < low_limit_epsilon] = low_limit_epsilon
            region_epsilon_scaled = np.log10(region_epsilon)

            if i == regions_out[-1]:

                print("\nDividing region", i, "into", num_zones_out, "zones ...\n")

                model_zone = AgglomerativeClustering(
                    linkage="ward", connectivity=region_connectivity,
                    n_clusters=num_zones_out)

                features = region_epsilon_scaled.reshape(-1, 1)

                scaler = StandardScaler()

                features_std = scaler.fit_transform(features)

                model_zone.fit(
                    features_std.astype(np.float32, casting='same_kind'))

                zone_id_field[region_condition] = np.int32(
                    model_zone.labels_ + num_target_inlets)

            else:
                zone_id_field[region_condition] = count
                count += np.int32(1)

        ''' Dividing the regions close to the inlets '''

        zone_ID_offset = num_target_inlets + num_zones_out
        for i in regions_in:
            region_condition = (regionLabels == i)
            region_cells = region_condition.nonzero()[0]

            connectivity_sliced = connectivity[region_cells, :]
            region_connectivity = connectivity_sliced[:, region_cells]

            region_p = p_Metal_NaOH[region_condition]
            region_mean = np.mean(region_p)

            region_epsilon = epsilon[region_condition]
            region_epsilon[region_epsilon < low_limit_epsilon] = low_limit_epsilon
            region_epsilon_scaled = np.log10(region_epsilon)

            print("\nDividing region", i, "into", num_zones_in, "zones ...\n")

            model_zone = AgglomerativeClustering(
                linkage="ward", connectivity=region_connectivity,
                n_clusters=num_zones_in)

            # Avoid cells with a very small values (small order of magnitude)
            # that are introduced in the previous regionalization due to
            # connectivity constraint
            region_p[region_p < low_limit_p] = low_limit_p

            region_p_scaled = scaling_func(region_p)

            # features = region_epsilon_scaled.reshape(-1, 1)
            # features = region_p_scaled.reshape(-1, 1)
            features = np.hstack((
                region_p_scaled.reshape(-1, 1),
                region_epsilon_scaled.reshape(-1, 1)))

            scaler = StandardScaler()

            features_std = scaler.fit_transform(features)

            model_zone.fit(
                features_std.astype(np.float32, casting='same_kind'))

            zone_id_field[region_condition] = np.int32(
                model_zone.labels_ + zone_ID_offset)

            zone_ID_offset += num_zones_in

        # print(max(zone_id_field))

        zone_id_path = os.path.join(OF_data_path, time_dir, 'zone_id')

        fieldName = 'p_env1'
        field_path = os.path.join(OF_data_path, time_dir, fieldName)

        num_header_rows = 16

        print("\nWriting zone information on the disk ...\n")
        with open(zone_id_path, 'wb') as f, open(field_path, 'rb') as headerSource:
            
            is_binary = False
            headerLines = list()
            for i in range(num_header_rows):
                line = headerSource.readline()

                if b'binary' in line:
                    is_binary = True

                if b'class' in line:
                    line = re.sub(
                        re.sub(b'\s', b'', line.split(b'class')[-1]),
                        b'volScalarField::Internal;', line)

                if b'location' in line:
                    line = re.sub(
                        re.sub(b'\s', b'', line.split(b'location')[-1]),
                        '\"{:s}\";'.format(time_dir).encode(), line)

                if b'object' in line:
                    line = re.sub(
                        re.sub(b'\s', b'', line.split(b'object')[-1]),
                        b'zone_id;', line)

                headerLines.append(line)

            f.writelines(headerLines)
            f.write(b'\ndimensions      [0 0 0 0 0 0 0];\n')
            f.write(b'\nvalue           nonuniform List<scalar>\n')
            f.write('{:d}\n'.format(num_cells).encode())

            if is_binary:
                f.write(b'(')
                f.write(struct.pack('{}d'.format(num_cells), *zone_id_field))
            else:
                f.write(b'(\n')
                np.savetxt(f, zone_id_field, fmt='%d')

            f.write(b');')
            f.write(b'\n\n\n')
            f.write(b'// ************************************************************************* //\n')

        # with open(zone_id_path, 'a') as f:
        #     np.savetxt(f, zone_id_field, fmt='%d')

            # f.write(')\n;\n\n\n')
            # f.write('// ************************************************************************* //\n')

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


def scaling_func(field):

    # scaled_field = np.log10(field / (1 - field))
    # scaled_field = (1.0 / np.sqrt(1 - field)) - (1.0 / np.sqrt(field))
    scaled_field = (1.0 / np.cbrt(1 - field)) - (1.0 / np.cbrt(field))

    return scaled_field


def construct_connectivity(
        cells, neighbors, owners, dataType=np.float32):
    cell_dict = dict()

    for index, cell in enumerate(cells):
        cell_dict[cell] = {'index': index}

    row = list()
    col = list()

    for (neighbor, owner) in zip(neighbors, owners):

        if neighbor in cell_dict and owner in cell_dict:
            neighbor_index = cell_dict[neighbor]['index']
            owner_index = cell_dict[owner]['index']

            if neighbor_index < owner_index:
                row.append(neighbor_index)
                col.append(owner_index)
            else:
                row.append(owner_index)
                col.append(neighbor_index)

    data = np.ones(len(row), dtype=dataType)

    target_size = cells.size

    connectivity_U = csr_matrix(
        (data,
         (np.asarray(row, dtype='int32'), np.asarray(col, dtype='int32'))),
        shape=(target_size, target_size))

    connectivity = (
        connectivity_U + connectivity_U.transpose()
        ).astype(dataType, casting='same_kind')

    if np.amin(connectivity.sum(axis=0)) < 1:
        print("\nError:\nDisconnected cells are detected\n")
        exit(1)

    return connectivity


# Python starts the program from here because it is the only executable code
# at indentation level 0. Python checks if the code is executed as a script
# and not imported by other scripts
if __name__ == "__main__":
    # If it is executed as script, it calls function "main" with arguments
    # entered by the user
    main(sys.argv[1:])
