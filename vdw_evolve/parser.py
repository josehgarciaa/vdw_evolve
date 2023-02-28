"""

"""

import wget
import os.path
import json

import ase.db
import pandas as pd
import numpy as np
from .physics import num2chem

def get_data_as_pd(path, options, props):
    """
    Example:
    options = 'is_magnetic=True, thermodynamic_stability_level=3'
    props = ["formula","spgnum", "spacegroup","uid","asr_id"]
    raw_df_1 = get_data_as_pd(path,options,props)

    :param path:
    :param options:
    :param props:
    :return:
    """
    db = ase.db.connect(path)
    # print("db",db)

    rows = db.select(options)
    print("rows", rows)
    data = [[x.get(p) for p in props] for x in rows]
    raw_df = pd.DataFrame(data, columns=props)
    return raw_df


# get json
def extract_json(uid_list, save_path="JSONcolection"):
    """

    :param uid_list:
    :param save_path:
    :return:
    """
    for uid in uid_list:
        data_url = 'https://cmrdb.fysik.dtu.dk/c2db/row/' + uid + '/all_data'
        file = save_path + "/" + uid + "_data.json"
        if os.path.isfile(file):
            print("file: ", file, "found")
        else:
            print(wget.download(data_url, out="./" + save_path))


def extract_structure(uid_list, save_path="STRUCTUREScolection"):
    for uid in uid_list:

        # https://cmrdb.fysik.dtu.dk/c2db/row/N2O2V3-358facb64a22/data/structure.json
        data_url = 'https://cmrdb.fysik.dtu.dk/c2db/row/' + uid + '/data/structure.json'
        print(data_url)
        file = save_path + "/" + uid + ".json"
        dw_path = save_path + "/" + "structure.json"
        if os.path.isfile(file):
            print("file: ", file, "found")
        else:
            print(wget.download(data_url, out="./" + save_path))
            os.rename(dw_path, file)


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


def get_cell_from_structure_file(f, xyz=False):
    lattice = get_lattice_from_structure_file(f)
    cell = np.array([[lattice[0][0], lattice[1][0]],
                     [lattice[0][1], lattice[1][1]]])
    if xyz == False:
        return cell

    xyz_data= get_atoms_from_structure_file(f)
    return cell, xyz_data


def get_atoms_from_structure_file(f):
    """
    :param f: "path/to/the/json/sith/structure"
    :return:
    """
    with open(f, 'r') as file:
        json_data = json.load(file)
    my_type = "__ndarray__"
    # print(json_data)
    positions = json_data['1']["positions"]['__ndarray__']
    positions = np.reshape(np.array(positions[2]), positions[0])
    type = json_data['1']["numbers"]['__ndarray__']
    type = np.reshape(np.array(type[2]), type[0])
    atoms = []
    for i, atom_n in enumerate(type):
        atom = [atom_n, [positions[i][0], positions[i][1], positions[i][2]]]
        atoms.append(atom)
    return atoms


def extract_molecule(xyz_path):
    molec = open(xyz_path).read()
    molec = molec.split()
    nr_atoms = int(molec[0])
    atoms = {}

    for atom_n in range(0, nr_atoms):
        atom_n = atom_n * 4 + 1
        atom = molec[atom_n]
        atoms[atom] = {}
        atoms[atom]["x"] = float(molec[atom_n + 1])
        atoms[atom]["y"] = float(molec[atom_n + 2])
        atoms[atom]["z"] = float(molec[atom_n + 3])

    return atoms

def t_units_to_x(molec, cel):

    # print("molec:", molec)
    a1 = np.array([cel[0,0],cel[1,0]])
    a2 = np.array([cel[1,0],cel[1,1]])
    new_m ={}
    for atom in molec:
        # print(a1*molec[atom]["x"])
        r = molec[atom]["x"]*a1 + molec[atom]["y"]*a2
        new_m[atom]=r.T
        new_m[atom]= np.append(new_m[atom],molec[atom]["z"])  #z
    # print("new_m:", new_m)
    return new_m