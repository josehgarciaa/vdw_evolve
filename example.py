from vdw_evolve import get_data_as_pd, extract_structure, get_cell_from_structure_file  # vdw.parser
from vdw_evolve import AnnealingSolver, MechanicSolver  # vdw.solvers
from vdw_evolve import allign_along_10  # vdw.solvers_utils

# Download data:

# select db file
database_path = "database/c2db-20211702.db"
structure_path = "database/2Dmaterials"

# Choose the selection criteria
options = 'is_magnetic=True, thermodynamic_stability_level=3'
props = ["formula", "spgnum", "spacegroup", "uid", "asr_id"]
raw_df = get_data_as_pd(database_path, options, props)
raw_df = raw_df[(raw_df["spacegroup"] != 'P1') & (raw_df["spacegroup"] != 'Pc') & (raw_df["spacegroup"] != 'P-1')]

# Print
print("The number of elements is:", len(list(raw_df)))
print(raw_df[["formula", "uid"]])

# Download
# extract_structure(raw_df["uid"], save_path=structure_path)


# Choose the two structures that will be overlapped
uid = ["Ru2F8-5b1d25d726e0", "V2F8-6d78fbe605b3"]  # "Pd2Se4-12f02221b8c5", "C2-a6735a4a3797"]
path1 = structure_path + "/" + uid[0] + ".json"
path2 = structure_path + "/" + uid[1] + ".json"

# extract the cell vectors from structure file
cell1 = get_cell_from_structure_file(path1)
cell2 = get_cell_from_structure_file(path2)

# Option step align the cells
cell1, cel2 = allign_along_10(cell1, cell2)
print("=== \n cel1:\n {}".format(cell1))
print("\n cel2:\n {}\n===".format(cell2))

# Select a solver and modify the solver parameters
solver1 = AnnealingSolver()
strain_boundary = [[-0.3, 0.3], [-0.3, 0.3]]
solver1.nr_epochs = 11

# Calculate super cell
super_cell = solver1.solve(cell1, cell2)
print(super_cell.description_txt())

# Experiment with other solvers
solver2 = MechanicSolver()
solver2.exploring_range = 50
solver2.tolerance = 0.1

# Calculate super cell
super_cell2 = solver2.solve(cell1, cell2)
print(super_cell2.description_txt())
