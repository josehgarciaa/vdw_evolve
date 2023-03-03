from vdw_evolve import get_data_as_pd, extract_structure, get_cell_from_structure_file  # vdw.parser
from vdw_evolve import AnnealingSolver, MechanicSolver, GeneticSolver  # vdw.solvers
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

# Download (only first 10 for the sake of example)
extract_structure(raw_df["uid"][:10], save_path=structure_path)

# Choose the two structures that will be overlapped
graphene= raw_df["uid"][0]#"C2-a6735a4a3797"
for u in raw_df["uid"]:
    path = "GRAPHENE_COLECTION/"+graphene+"__"+u
    uid = [graphene, u]
    path1 = structure_path + "/" + uid[0] + ".json"
    path2 = structure_path + "/" + uid[1] + ".json"

    # extract the cell vectors from structure file
    cell1 = get_cell_from_structure_file(path1)
    cell2 = get_cell_from_structure_file(path2)

    # Option step align the cells
    aligned_cells = allign_along_10([cell1, cell2])
    cel1 = aligned_cells[0]
    cel2 = aligned_cells[1]
    print("=== \n cel1:\n {}".format(cell1))
    print("\n cel2:\n {}\n===".format(cell2))


    # Select a solver and modify the solver parameters

    print("\n\nAnnealingSolver:")
    solver1 = AnnealingSolver()
    solver1.model_par["strain_boundary"] = [[-0.2, 0.2], [-0.2, 0.2]]
    solver1.nr_epochs = 11

    # Calculate super cell
    super_cell1 = solver1.solve(cell1, cell2)
    print(super_cell1.description_txt())
    super_cell1.description_txt(path +"_Annealing_"+".json")

    # # Experiment with other solvers

    print("\n\nGeneticSolver:")
    solver2 = GeneticSolver()
    solver2.model_par["strain_boundary"] = [[-0.2, 0.2], [-0.2, 0.2]]
    solver2.model_par["subjects_in_cell"] = 2
    solver2.model_par["pins"] = 4

    # Calculate super cell
    super_cell2 = solver2.solve(cell1, cell2)
    print(super_cell2.description_txt())
    super_cell2.description_txt(path +"_Genetic_"+".json")


