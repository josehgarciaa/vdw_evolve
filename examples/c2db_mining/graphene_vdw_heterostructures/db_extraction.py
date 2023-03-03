from vdw_evolve import get_data_as_pd, extract_structure, get_cell_from_structure_file  # vdw.parser
from vdw_evolve import SuperMatcher
from vdw_evolve.output_structure import TemporaryStructure
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
extract_structure(raw_df["uid"][:100], save_path=structure_path)#TODO: take just first 100

# Choose the two structures that will be overlapped

#TODO: Set the main material to be grapheen
main_material = raw_df["uid"][0] #"C2-a6735a4a3797"

for i, material in enumerate(raw_df["uid"]):
    print("{}/{}".format(i, len(raw_df["uid"])))


    # This part just builds the path to the json files
    uid = [material, main_material ]
    path1 = structure_path + "/" + uid[0] + ".json"
    path2 = structure_path + "/" + uid[1] + ".json"


    # Extract the cell vectors from structure file
    cell1, xyz1 = get_cell_from_structure_file(path1, True)
    cell2, xyz2 = get_cell_from_structure_file(path2, True)

    # This is just your solver should be fine
    solver0 = SuperMatcher()
    super_cell0,s1,s2,t = solver0.solve(cell1, cell2, ch=True)
    cell = super_cell0.cell

    # Puts the output in the right format
    # TODO: Please check  the output_structure.py
    ts = TemporaryStructure(cell,cell1, cell2, xyz1,xyz2,s1,s2,t)
    path = "SUPERCELLS/"+uid[0]+"__"+uid[1]+"scell.json"
    ts.json(path)
    print("\n\n")

