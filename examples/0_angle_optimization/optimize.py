import vdw_evolve as lt
import vdw_evolve.geometry as geo

# Read Structure files
str1 = lt.Structure().read_from("C2_0deg.xyz", format="c2db-xyz")
str2 = lt.Structure().read_from("C2_30deg.xyz", format="c2db-xyz")
vdws = lt.VdWStructure(str1, str2)

supercell_params = vdws.get_minimalcell(dims=(1000 ,1000))

print(supercell_params)
