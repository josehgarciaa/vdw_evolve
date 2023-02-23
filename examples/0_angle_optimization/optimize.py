import vdw_evolve as lt


# Read Structure files
str1 = lt.Structure().read_from("C2_0deg.xyz", format="c2db-xyz")
str2 = lt.Structure().read_from("C2_45deg.xyz", format="c2db-xyz")
vdws = lt.VdWStructure(str1, str2)

# Optimize using a single core
vdw_str = lt.VdWStack( vdws, opt_strain=False, opt_angle=True, ncores=1);
#vdw_str.write("vdw_optimal.xyz");
#vdw_str.write("vdw_optimal.json");

#Plot Structures
#vdw_structure.plot(ndim=(2,2));


