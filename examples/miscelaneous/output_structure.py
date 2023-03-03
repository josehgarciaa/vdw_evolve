"""
This is just a temporary solution
Andrei  will build another one when he  will clean the repo
maybe next weekend.
"""

import numpy as np
import json
from .solvers_utils.lat_macher import ChangeBasis


class TemporaryStructure:

    def __init__(self, cell, cell1, cell2, xyz1, xyz2, s1, s2, theta):
#TODO: nwhen xyz coordinates are taken from the json
# they are in carteesian coordinates  plese check to be sure that i do not messed up.
        self.cell = cell
        self.cell1 = cell1
        self.cell2 = cell2
        self.xyz1 = xyz1
        self.xyz2 = xyz2
        self.s1 = s1
        self.s2 = s2
        self.theta = theta

    def json(self, path):
        """
        Constructs the json
        :param path:
        :return:
        """
        # extract the symbols and positions for the parents
        symbol1 = [s[0] for s in self.xyz1]
        symbol2 = [s[0] for s in self.xyz2]
        positions1 = [p[1] for p in self.xyz1]
        positions1 = np.array(positions1).reshape(-1).tolist()
        positions2 = [p[1] for p in self.xyz2]
        positions2 = np.array(positions2).reshape(-1).tolist()

        # TODO: This i s the part that required check:
        s_symbol, s_points = self.xyz()

        # Build the Json
        js = {"1": {
            "cell": {"__ase_objtype__": "cell",
                     "array": {"__ndarray__": [[3, 3], "float64", border(self.cell).T.reshape(-1).tolist()]}},
            "numbers": {"__ndarray__": [[int(len(s_symbol))], "int64", s_symbol]},
            "positions": {"__ndarray__": [[int(len(s_symbol)), 3], "float64", np.array(s_points).reshape(-1).tolist()]},
        },

            "2": {

                "p1": {
                    "cell": {"__ase_objtype__": "cell",
                             "array": {"__ndarray__": [[3, 3], "float64", border(self.cell1).T.reshape(-1).tolist()]}},
                    "numbers": {"__ndarray__": [[int(len(s_symbol))], "int64", symbol1]},
                    "positions": {"__ndarray__": [[int(len(s_symbol)), 3], "float64", positions1]},
                },

                "p2": {
                    "cell": {"__ase_objtype__": "cell",
                             "array": {"__ndarray__": [[3, 3], "float64", border(self.cell2).T.reshape(-1).tolist()]}},
                    "numbers": {"__ndarray__": [[int(len(s_symbol))], "int64", symbol2]},
                    "positions": {"__ndarray__": [[int(len(s_symbol)), 3], "float64", positions2]},
                },

                "transformations": {
                    "s1": self.s1,
                    "s2": self.s2,
                    "theta": self.theta
                }

            }}

        if path != None:
            # json_object = json.dumps(desc, indent=4)
            with open(path, "w") as outfile:
                json.dump(js, outfile, cls=NpEncoder)
        return js

    def xyz(self):
        """
        Builds a grid with atoms from lattice 1 and 2 and
        the select the atoms tha tare in the minimum super cell
        :return:
        """

        # TODO: Check
        s_points, s_symbol = [], []
        l_points, l_symbol = self.stack_2D_xyz()
        xy_l_pints = np.array([ChangeBasis(x[:2], self.cell) for x in l_points])
        # xy_l_pints = ChangeBasis(xy_l_pints.T, self.cell)
        for i, point in enumerate(xy_l_pints):
            # print("point:", point)
            if (0 < point[0] <= 1) & (0 < point[1] <= 1):
                s_points.append(l_points[i])
                s_symbol.append(l_symbol[i])
        # print(s_symbol)
        return s_symbol, s_points

    def stack_2D_xyz(self, zdis=10):
        lat_pointsB, lat_symbolsB = self.xyz2supercell(self.cell1, self.xyz1)
        lat_pointsT, lat_symbolsT = self.xyz2supercell(self.cell2, self.xyz2)

        lat_points = np.array(list(lat_pointsB) + list(lat_pointsT + [0, 0, zdis]))
        lat_symbols = np.array(list(lat_symbolsB) + list(lat_symbolsT))
        return lat_points, lat_symbols

    def xyz2supercell(self, lat_vec, xyz, shape=(1, 1, 1), no_expand=False):
        # Create a grid corresponding to the points in the
        # supercell (sc)

        print("lat_vec", lat_vec)
        lat_vec = np.array([[lat_vec[0][0], lat_vec[0][1], 0],
                            [lat_vec[1][0], lat_vec[1][1], 0],
                            [0, 0, 1]])
        nx, ny, nz = shape
        grid = np.mgrid[-nx:nx, -ny:ny, 0:nz].reshape(3, 2 * nx * 2 * ny * nz).T

        if no_expand is True:
            grid = np.mgrid[0:1, 0:1, 0:1].reshape(3, 1).T

        sc_points = grid.dot(lat_vec)

        # Shift all points in the xyz file to the differnt cells in the supercell (sc)
        symbols, points = list(zip(*xyz))
        lat_points = [list(sc_points + point) for point in points];
        lat_points = np.array(lat_points).reshape(len(xyz) * len(sc_points), 3)
        lat_symbols = np.repeat(symbols, len(sc_points))  # We repeated them following the order
        return lat_points, lat_symbols


class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


def border(cell):
    a = np.array([[cell[0][0], cell[0][1], 0], [cell[1][0], cell[1][1], 0], [0, 0, 1]])
    return a
