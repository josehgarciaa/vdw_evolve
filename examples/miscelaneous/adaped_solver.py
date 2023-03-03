
import numpy as np

from vdw_evolve import AnnealingSolver
from vdw_evolve.solvers_utils.fit_function import rectangle_fit_function

from vdw_evolve import allign_along_10


# Example
graphene_cell = np.array([[2.467, 2.467 / 2], [0, 2.467 * np.sqrt(3) / 2]])
cell = allign_along_10([graphene_cell])[0]

print("\n\nAnnealing_r_Solver:")
solver1 = AnnealingSolver()
solver1.fit_f = rectangle_fit_function
solver1.nr_epochs = 9
super_cell_annealing = solver1.solve(graphene_cell, np.array([[1,0],[0,1]]))
print(super_cell_annealing.description_txt())