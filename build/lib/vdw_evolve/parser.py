"""

"""

import wget
import os.path
import json

import ase.db
import pandas as pd
import numpy as np


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


def get_lattice_from_structure_file(f):
    """

    :param f: "path/to/the/json/sith/structure"
    :return:
    """
    with open(f, 'r') as file:
        json_data = json.load(file)
    my_type = "__ndarray__"
    data_struct = json_data['1']["cell"]["array"]
    if my_type in data_struct:
        shape = data_struct[my_type][0]
        d_type = data_struct[my_type][1]
        lattice = np.array(data_struct[my_type][2], dtype=d_type).reshape(shape)
    return lattice


def get_cell_from_structure_file(f):
    lattice = get_lattice_from_structure_file(f)
    cell = np.array([[lattice[0][0], lattice[1][0]],
                [lattice[0][1], lattice[1][1]]])
    return cell


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

    print("molec:", molec)
    a1 = np.array([cel[0,0],cel[1,0]])
    a2 = np.array([cel[1,0],cel[1,1]])
    new_m ={}
    for atom in molec:
        print(a1*molec[atom]["x"])
        r = molec[atom]["x"]*a1 + molec[atom]["y"]*a2
        new_m[atom]=r.T
        new_m[atom]= np.append(new_m[atom],molec[atom]["z"])  #z
    print("new_m:", new_m)
    return new_m