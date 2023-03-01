import vdw_evolve as lt
import vdw_evolve.geometry as geo

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt

# Read Structure files
str1 = lt.Structure().read_from("c2db-PtSe2.json", format="c2db-json")
str2 = lt.Structure().read_from("c2db-PtSe2_strained_0.021.xyz", format="c2db-xyz")
vdws = lt.VdWStructure(str1, str2)


# Define the type of optimization
optimizer = lt.LatMatch().UniformStrain(0.1).OptAngle(False)

# Get the minimal cell
dims = (10,10)
optVdW = vdws.get_minimalcell(dims=dims, optimizer=optimizer)
print("Strain: ", optVdW.strain, "angle: ", optVdW.angle)


plt.axes().set_aspect('equal')
plt.scatter(*(optVdW.host.supercell_points(dims)[:2]), s=16.0)
plt.scatter(*(optVdW.complement.supercell_points(dims)[:2]), s=16.0)

# Parallelogram
cell = optVdW.cell[:2,:2]
x = [0.0,cell[0][0],cell[0][0]+cell[0][1],cell[0][1],0.0];
y = [0.0,cell[1][0],cell[1][0]+cell[1][1],cell[1][1],0.0];
plt.gca().add_patch(patches.Polygon(xy=list(zip(x, y)), fill=True, alpha=0.4, color="green"))
plt.savefig('matched_lattice.pdf')  

