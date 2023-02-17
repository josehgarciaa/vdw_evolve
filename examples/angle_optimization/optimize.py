import latt_match as lt


#Read Structure files
str1 = lt.Structure("C2_0deg.json", format="c2db");
str2 = lt.Structure("C2_0deg.xyz" , format="c2db");

#Optimize using a single core
vdw_str = Optimize( str1, str2, opt_strain=False, opt_angle=True, ncores=1);
vdw_str.write("vdw_optimal.xyz");
vdw_str.write("vdw_optimal.json");

#Plot Structures
vdw.plot(ndim=(2,2));


