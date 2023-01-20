from vdw_evolve import get_data_as_pd, extract_structure, get_cell_from_structure_file  # vdw.parser
from vdw_evolve import AnnealingSolver, MechanicSolver, GeneticSolver  # vdw.solvers
from vdw_evolve import allign_along_10  # vdw.solvers_utils

import numpy as np

# Download data:



# extract the cell vectors from structure file
cell1 = np.array([[2.467, 2.467/2],[0, 2.467*np.sqrt(3)/2]])
cell2 = np.array([[5.75,0],[0,5.92]])

# Option step align the cells
cells =  allign_along_10([cell1, cell2])
cell1 = cells[0]
cell2 = cells[1]
print("=== \n cel1:\n {}".format(cell1))
print("\n cel2:\n {}\n===".format(cell2))


# Select a solver and modify the solver parameters

print("\n\nAnnealingSolver:")
solver1 = AnnealingSolver()
strain_boundary = [[-0.3, 0.3], [-0.3, 0.3]]
solver1.nr_epochs = 11

# Calculate super cell
super_cell1 = solver1.solve(cell1, cell2)
print(super_cell1.description_txt())


# # Experiment with other solvers

print("\n\nGeneticSolver:")
solver2 = GeneticSolver()
strain_boundary = [[-0.3, 0.3], [-0.3, 0.3]]
solver2.model_par["subjects_in_cell"] = 2
solver2.model_par["pins"] = 4

# Calculate super cell
super_cell2 = solver2.solve(cell1, cell2)
print(super_cell2.description_txt())

print("\n\nMechanicSolver:")
solver3 = MechanicSolver()
solver3.exploring_range = 50
solver3.tolerance = 0.4

# Calculate super cell
super_cell3 = solver3.solve(cell1, cell2)
print(super_cell3.description_txt())


# Select a solver and modify the solver parameters
