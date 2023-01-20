"""

"""

from .solvers_utils import annealing_sc, mechanic_sc, genetic_sc
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
        ps =True
        while ps:
            try:
                ta, strain_tb, t_cel2_no_strain, diagonal_strain, strain = annealing_sc(cel1, cel2, self.nr_epochs,
                                                                            self.model_par, self.fit_f)
                ps=False
            except:
                pass

        sc = SuperCell(parents=(cel1, cel2), transformation=(ta, strain_tb), strains=(strain, diagonal_strain))
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
