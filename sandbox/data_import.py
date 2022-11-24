import wget
import os.path
import json

import ase.db
import pandas as pd
from numpy import array


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
    rows = db.select(options)
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


def get_lattice_from_structure(f):
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
        lattice = array(data_struct[my_type][2], dtype=d_type).reshape(shape)
    return lattice
