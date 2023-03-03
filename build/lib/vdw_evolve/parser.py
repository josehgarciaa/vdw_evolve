import json
import numpy as np
from .physics import num2chem

def cell_fromJSON(filepath, format="c2db-json"):
    """ Read the cell out of a json file in different formats

    Parameters
    ----------
    filepath : str
        Path to the json file
    format : str
        Format of the json file. 
        Default is c2db-json

    Returns
    -------
    array-like
        The unit cell of structure, with vectors as columns.
    """
    if format == "c2db-json":
        with open(filepath, 'r') as file:
            json_data = json.load(file)
        arr_type = "__ndarray__"
        try:
            data_struct = json_data['1']["cell"][arr_type]
        except ValueError:
            print("The array type is not the proper one")
        shape = data_struct[0]
        d_type = data_struct[1]
        lattice = np.array(data_struct[2], dtype=float).reshape(shape)
        return np.transpose(lattice.reshape(shape))
    return None


def atoms_fromJSON(filepath, format="c2db-json"):
    """ Read the atoms from a json file in different formats

    Parameters
    ----------
    filepath : str
        Path to the json file
    format : str
        Format of the json file. 
        Default is c2db-json

    Returns
    -------
    array-like
        The unit cell of structure, with vectors as columns.
    """

    if format == "c2db-json":
        with open(filepath, 'r') as file:
            json_data = json.load(file)
        arr_type = "__ndarray__"
        
        try:
            positions = json_data['1']["positions"][arr_type]
            positions = np.reshape(np.array(positions[2]), positions[0])
        except ValueError:
            print("The array type is not the proper one")
        
        try:
            type = json_data['1']["numbers"][arr_type]
            type = np.reshape(np.array(type[2]), type[0])
        except ValueError:
            print("The array type is not the proper one")
        atoms = []
        for i, atom_n in enumerate(type):
            atom = [num2chem[atom_n], [positions[i][0], positions[i][1], positions[i][2]]]
            atoms.append(atom)
        return atoms


def cell_fromXYZ(filepath, format):
    """ Read the cell from a xyz when allowed

    Parameters
    ----------
    filepath : str
        Path to the xyz file
    format : str
        Format of the xyz file. 
        Default is c2db-xyz

    Returns
    -------
    array-like
        The unit cell of structure, with vectors as columns.
    or None
        When the lattice vector cannot be extracted
    """
    if format == "c2db-xyz":
        with open(filepath) as f:
            npoints = f.readline()
            lat = f.readline()
            lat = lat[lat.find("\"")+1:lat.find("Properties")-2]
            lat = lat.split(" ")
            try:
                cell = np.transpose(np.array(list(map(float, lat))).reshape(3, 3))
                return cell
            except ValueError:
                print("Could not properly parse input_file from cells")
                return None


def atoms_fromXYZ(filepath, format):
    """ Read the atoms from a xyz when allowed

    Parameters
    ----------
    filepath : str
        Path to the xyz file
    format : str
        Format of the xyz file.
        Default is c2db-xyz

    Returns
    -------
    array-like
        The unit cell of structure, with vectors as columns.
    or None
        When the lattice vector cannot be extracted 
    """
    if format == "c2db-xyz":
        def xyz_format(line):
            s, x, y, z = [x for x in line.split(" ") if x != ""][:4]
            return (s, np.array([float(x), float(y), float(z)]))

        with open(filepath) as f:
            npoints = f.readline()
            lat = f.readline()
            try:
                atoms = [xyz_format(line) for line in f]
                return atoms
            except ValueError:
                print("Could not properly parse input_file from xyz")
                return None