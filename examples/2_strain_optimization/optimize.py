import vdw_evolve as lt
import vdw_evolve.geometry as geo

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt


# Read Structure files
str1 = lt.Structure().read_from("C2_s0.xyz", format="c2db-xyz")
str2 = lt.Structure().read_from("C2_s0.1_0.3.xyz", format="c2db-xyz")
vdws = lt.VdWStructure(str1, str2)


# Define the type of optimization
optimizer = lt.LatMatch().opt_angle(False)
optimizer.set_max_strain((0.5, 0.5))

# Get the minimal cell
dims = (10,10)
optVdW = vdws.get_minimalcell(dims=dims, optimizer=optimizer)
print("Strain: ", optVdW.complement_strain(),
      "angle: ", optVdW.complement_angle())

plt.axes().set_aspect('equal')
plt.scatter(*(optVdW.host.supercell_points( dims )[:2]), s=16.0)
plt.scatter(*(optVdW.complement.supercell_points( dims )[:2]), s=16.0)

# Parallelogram
cell = optVdW.cell[:2,:2]
x = [0.0,cell[0][0],cell[0][0]+cell[0][1],cell[0][1],0.0];
y = [0.0,cell[1][0],cell[1][0]+cell[1][1],cell[1][1],0.0];
plt.gca().add_patch(patches.Polygon(xy=list(zip(x,y)), fill=True, alpha=0.4,color="green"))
plt.savefig('matched_lattice.pdf')  

