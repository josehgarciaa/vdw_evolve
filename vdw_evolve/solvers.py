"""

"""

import numpy as np
import scipy as sp

from .solvers_utils import annealing_sc, mechanic_sc, genetic_sc, LatMatch
from .solvers_utils.fit_function import base_fit_function
from .structure import SuperCell


class AnnealingSolver:

    def __init__(self, nr_epochs=10, bond=30, strain_boundary=[[-0.3, 0.3], [-0.3, 0.3]]):
        self.nr_epochs = nr_epochs
        self.model_par = {
            "bounds": [[-bond, bond] for _ in range(4)],
            "strain_boundary": strain_boundary,
            "start_point": [1, 5, 3, 2],
            "known_min": -999999,

            "initialTemp": 4,
            "finalTemp": 0.0002,
            "beta": 10,

            "nr_neighbours": 1,
            "step_size": 3,
            "gaussian_sd": 3
        }
        self.fit_f = base_fit_function

    def solve(self, cel1, cel2):
        # cel1 = np.dot(cel)
        ps = True
        while ps:
            try:
                ta, strain_tb, t_cel2_no_strain, diagonal_strain, strain = annealing_sc(cel1, cel2, self.nr_epochs,
                                                                                        self.model_par, self.fit_f)
                # print("ps:", ps)
                ps = False
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                pass

        sc = SuperCell(parents=(cel1, cel2), transformation=(ta, strain_tb), strains=(strain, diagonal_strain))
        return sc


class SuperMatcher:
    opt_angle = True
    opt_strain = True
    # theta_min= 15*np.pi/180;
    theta_min = 5 * np.pi / 180
    theta_range = (-np.pi / 2, np.pi / 2)
    smax = [0.05, 0.05]
    bounds = None
    result = None

    def __init__(self, scdim = [5,5], optimize_angle=True, optimize_strain=True):
        self.dim = scdim
        self.updated_target_cell = None
        self.sc_vec = None
        self.opt_angle = optimize_angle
        self.opt_strain = optimize_strain



    def solve(self, cel1, cel2, ch =False):
        # Todo: Please check!!:
        matcher1=LatMatch(scdim=self.dim, reference=cel1, target=cel2, optimize_strain=self.opt_strain,optimize_angle=self.opt_angle);
        sc_vec = matcher1.supercell().T
        print(sc_vec)
        oBlat = sc_vec#matcher1.updatedCell()
        print("\n",oBlat)

        result= matcher1.Result() # S1 S2 Theta
        s1=result[0]
        s2=result[1]
        theta=result[2]

        r = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        strain = np.diag([1 + s1, 1 + s2])

        diagonal_strain = strain

        ta= np.dot(oBlat, np.linalg.inv(cel1))
        strain_tb= np.dot(matcher1.updatedCell(), np.linalg.inv(cel2))

        # strain_tb = strain @ r
        #oBlat = Ta*cel1= TB*cel2

        sc = SuperCell(parents=(cel1, cel2), transformation=(ta, strain_tb), strains=(strain, diagonal_strain))
        if ch:
            return sc, s1,s2,theta
        return sc





class GeneticSolver:

    def __init__(self, nr_epochs=10, bond=30, strain_boundary=[[-0.3, 0.3], [-0.3, 0.3]]):
        self.nr_epochs = nr_epochs
        self.model_par = {
            "bounds": [[-bond, bond] for _ in range(4)],
            "strain_boundary": strain_boundary,
            "start_point": [1, 5, 3, 2],
            'cell_split_number': 2,
            'subjects_in_cell': 2,
            '0_in_cell_enhancement': 2,

            'nr_clones': 10,
            'mutation_gaussian_sd': 4,

            'pins': 7,
            'gene_quality': 1,
            "known_min": -999999,
        }
        self.fit_f = base_fit_function

    def solve(self, cel1, cel2):
        # cel1 = np.dot(cel)
        ps = True
        while ps:
            try:
                ta, strain_tb, t_cel2_no_strain, diagonal_strain, strain = genetic_sc(cel1, cel2, self.nr_epochs,
                                                                                      self.model_par, self.fit_f)
                ps = False
            except:
                pass
        sc = SuperCell(parents=(cel1, cel2), transformation=(ta, strain_tb), strains=(strain, diagonal_strain))
        return sc


class MechanicSolver:

    def __init__(self, exploring_range=10, tolerance=0.1):
        self.exploring_range = exploring_range
        self.tolerance = tolerance
        self.paralel_limit = 0.00001

    def solve(self, cel1, cel2):
        ta, strain_tb, diagonal_strain, strain = mechanic_sc(cel1, cel2, exploring_range=self.exploring_range,
                                                             tolerance=self.tolerance, paralel_limit=self.paralel_limit)
        sc = SuperCell(parents=(cel1, cel2), transformation=(ta, strain_tb), strains=(strain, diagonal_strain))
        return sc
