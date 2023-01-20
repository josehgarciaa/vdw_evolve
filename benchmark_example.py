import numpy as np

from vdw_evolve import allign_along_10
from vdw_evolve import AnnealingSolver, MechanicSolver, GeneticSolver
from vdw_evolve import SuperLattice

graphene_xyz_path = "zx_files/graphene.xyz"
pbse_xyz_path = "zx_files/PdSe2.xyz"

graphene_cell = np.array([[2.467, 2.467 / 2], [0, 2.467 * np.sqrt(3) / 2]])
pbse_cell = np.array([[5.75, 0], [0, 5.92]])

aligned_cells = allign_along_10([graphene_cell, pbse_cell])
graphene_cell = aligned_cells[0]
pbse_cell = aligned_cells[1]

print("\n\nAnnealingSolver:")
solver1 = AnnealingSolver()
solver1.model_par["strain_boundary"] = [[-0.05, 0.05], [-0.05, 0.05]]
solver1.nr_epochs = 11
super_cell_annealing = solver1.solve(graphene_cell, pbse_cell)
annealing_lattice = SuperLattice(super_cell_annealing, graphene_xyz_path, pbse_xyz_path)
print(super_cell_annealing.description_txt())

print("\n\nGeneticSolver:")
solver2 = GeneticSolver()
solver2.model_par["strain_boundary"] = [[-0.05, 0.05], [-0.05, 0.05]]
solver2.model_par["subjects_in_cell"] = 2
solver2.model_par["pins"] = 4
# Calculate super cell
super_cell_genetic = solver2.solve(graphene_cell, pbse_cell)
genetic_lattice = SuperLattice(super_cell_genetic, graphene_xyz_path, pbse_xyz_path)
print(super_cell_genetic.description_txt())

print("\n\nMechanicSolver:")
solver3 = MechanicSolver()
solver3.exploring_range = 50
solver3.tolerance = 0.4
# Calculate super cell
super_cell_mechanic = solver3.solve(graphene_cell, pbse_cell)
mechanic_lattice = SuperLattice(super_cell_mechanic, graphene_xyz_path, pbse_xyz_path)
print(super_cell_mechanic.description_txt())


print("\nNr atoms  in annealing_lattice:", annealing_lattice.nr_atoms)
print("Nr atoms  in genetic_lattice  :", genetic_lattice.nr_atoms)
print("Nr atoms  in mechanic_lattice :", mechanic_lattice.nr_atoms)



