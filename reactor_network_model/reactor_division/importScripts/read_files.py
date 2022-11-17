import os
import re
import sys
import traceback
import struct
import re
from unittest.mock import patch
import numpy as np

from enum import Enum, auto


class Format(Enum):
    ASCII = auto()
    BINARY = auto()


def read_field(filePath, read_boundary=True, num_headerLines=16):
    try:
        with open(filePath, "rb") as f:
            lines = f.readlines()

        num_headerLines = max(num_headerLines, 21)

        headerLines = lines[:num_headerLines]

        data_format = find_format(headerLines)

        if has_boundary(headerLines):

            lines_boundary = split_bField(lines)

            if read_boundary:
                return read_field_internal(lines, data_format), \
                    read_field_boundary(lines_boundary, data_format)
            else:
                return read_field_internal(lines, data_format), {}

        else:
            return read_field_internal(lines, data_format), {}

    except Exception as e:
        print(str(e))
        traceback.print_exc()
        sys.exit(1)


def read_field_internal(lines, data_format):

    for row, line in enumerate(lines):
        if (b'internalField' in line) or (b'value' in line):
            break
    
    internal_field = None

    if b'nonuniform' in line:
        field_size = int(re.sub(b'\W', b'', lines[row + 1]))

        if data_format is Format.ASCII:
            offset = row + 3
            internal_field = np.array(
                lines[offset : offset + field_size], dtype=float)

        else:
            offset = row + 2
            buffer = b''.join(lines[offset: ])

            internal_field = np.frombuffer(
                buffer, dtype='d', count=field_size, offset=1)

    elif b'uniform' in line:
        return float(re.sub(b'\W', b'', line.split(b'uniform')[-1]))

    else:
        raise Exception("wrong data format")

    return internal_field


def read_field_boundary(lines, data_format):

    boundary_data = dict()

    for i, line in enumerate(lines):
        if b'{\n' in line:
            boundary_data = read_dict_bField(
                lines[i + 1: ], 'boundaryField', data_format)
            break

    return boundary_data


def read_dict_bField(lines, dictName, data_format, start=0):

    dict_ = dict()

    for i, _ in enumerate(lines, start):
        line = lines[i]
        # print(dictName, i, line)
        if b'{\n' in line:
            key = re.sub(b'\s', b'', lines[i - 1]).decode()
            lines[i - 1] = b''
            lines[i] = b''
            dict_[key] = read_dict_bField(lines, key, data_format, i+1)

        elif b'type' in line:
            value = re.sub(b'\W', b'', line.split(b'type')[-1]).decode()
            dict_['type'] = value
            lines[i] = b''

        elif b'value' in line:
            if b'nonuniform' in line:
                bField_size = int(re.sub(b'\W', b'', lines[i + 1]))
                # print(bField_size)

                if data_format is Format.ASCII:
                    offset = i + 3
                    value = np.array(
                        lines[offset : offset + bField_size], dtype=float)

                else:
                    offset = i + 2
                    buffer = b''.join(lines[offset: ])

                    value = np.frombuffer(
                        buffer, dtype='d', count=bField_size, offset=1)

                dict_['value'] = value

            elif b'uniform' in line:
                value = float(re.sub(b'\W', b'', line.split(b'uniform')[-1]))
                dict_['value'] = value

            lines[i] = b''

        elif b'}\n' in line:
            lines[i] = b''
            return dict_

    return dict_

def split_bField(lines):
    for row, line in enumerate(lines):
        if b'boundaryField' in line:
            break

    if row == len(lines) - 1:
        raise Exception("Boundary data is missing")

    begin = row

    lines_sliced = lines[row: ]

    for row, line in enumerate(lines_sliced):
        if b'{' in line:
            break
    
    if row == len(lines_sliced) - 1:
        raise Exception("\"{\" is expected after boundaryField}")
    
    count_left_curl = 1
    count_right_curl = 0

    while count_left_curl != count_right_curl:
        row += 1
        try:
            line = lines_sliced[row]

        except:
            Exception("boundaryField is corrupted")

        count_left_curl += line.count(b'{\n')
        count_right_curl += line.count(b'}\n')

    end = begin + row

    boundary_lines = lines[begin: end]

    del lines[begin: end]

    return boundary_lines


def find_format(headerLines):
    for line in headerLines:
        if b'binary' in line:
            return Format.BINARY

    return Format.ASCII


def has_boundary(headerLines):
    for line in headerLines:
        if b'Internal' in line:
            return False

    return True


def read_owner_nbr(filePath, num_headerLines=16):
    cell_indexes = None
    try:
        lines = None
        with open(filePath, 'rb') as f:
            lines = f.readlines()

        headerLines = lines[0: num_headerLines]

        data_format = find_format(headerLines)

        match_word = None
        for line in headerLines:
            if b'object' in line:
                cell_type = re.sub(b'\W', b'', line.split(b'object')[-1])

                if cell_type == b'neighbour':
                    match_word = b'nInternalFaces'

                elif cell_type == b'owner':
                    match_word = b'nFaces'

                else:
                    raise Exception(
                        '\nUnknown object: {} \n expected \"owner\" or '
                        '\"neighbour\"\n'.format(cell_type))

                break
        
        if match_word == None:
            raise Exception('\nObject not found in header of file {}\n'.format(
                os.path.basename(filePath)))

        num_faces = None
        for line in headerLines:
            if match_word in line:
                notes = re.sub(b':', b' ', line).split()
                for i, note in enumerate(notes):
                    if match_word in note:
                        num_faces = int(re.sub(b'\W', b'', notes[i + 1]))
                        break
                break

        for i, line in enumerate(lines):
            if bytes(str(num_faces), 'ascii') in line and b'(' in lines[i + 1]:
                if data_format is Format.ASCII:
                    lines_to_convert = lines[i + 2: i + 2 + num_faces]

                    cell_indexes = np.array(
                        [line.strip() for line in lines_to_convert], np.int32)
                else:
                    buffer = b''.join(lines[i + 1: ])

                    cell_indexes = np.frombuffer(
                        buffer, dtype=np.int32, count=num_faces, offset=1)

                break

    except Exception as e:
        print(str(e))
        traceback.print_exc()
        sys.exit(1)

    return cell_indexes


def read_boundary(filePath):
    boundary = dict()
    try:
        lines = None
        with open(filePath, 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if '(' in line:
                boundary = read_dict_boundary(lines[i + 1: ], 'boundary')
                break

    except Exception as e:
        print(str(e))
        traceback.print_exc()
        sys.exit(1)

    return boundary


def read_dict_boundary(lines, dictName, start=0):
    dict_ = dict()

    for i, _ in enumerate(lines, start):
        line = lines[i]
        if '{' in line:
            key = re.sub('\s', '', lines[i - 1])
            lines[i - 1] = ''
            lines[i] = ''
            dict_[key] = read_dict_boundary(lines, key, i+1)

        elif '}' in line:
            lines[i] = ''
            return dict_

        else:
            for key in ['nFaces', 'startFace']:
                if key in line:
                    value = int(re.sub('\W', '', line.split(key)[-1]))
                    dict_[key] = value
                    lines[i] = ''

    return dict_


if __name__ == "__main__":

    import inspect
    import os

    # Find location of the executed script
    frame = inspect.currentframe()
    filePath = inspect.getfile(frame)
    caseDir = os.path.realpath(os.path.abspath(os.path.dirname(filePath)))

    # Path to the import files
    OF_data_path = os.path.join(caseDir, '..', '../cfdSimulation')

    time_dir = '0'

    field = 'phi'
    # field = 'supersaturation'

    file_path = os.path.join(OF_data_path, time_dir, field)

    internal_field, boundary_field = read_field(file_path, read_boundary=True)

    print(internal_field[:2], internal_field.size)
    print(boundary_field)
