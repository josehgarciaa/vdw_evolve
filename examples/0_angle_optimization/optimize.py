import vdw_evolve as lt
import vdw_evolve.geometry as geo

# Read Structure files
str1 = lt.Structure().read_from("C2_0deg.xyz", format="c2db-xyz")
str2 = lt.Structure().read_from("C2_0deg.xyz", format="c2db-xyz")
vdws = lt.VdWStructure(str1, str2)


import numpy as np
vector_list = vdws.supercell_points();
tol = 1e-7
alpha=0.5;

V = np.transpose(vector_list).T
basis = []
#cost = np.inf
#for i in range(len(V)):
#        v, ws = V[i], V[i+1:]
#        areas = np.abs(np.cross([v[:2]], ws[:, :2]))
#        print(areas)
#        ws = ws[areas > min_area]
#        areas = areas[areas > min_area]
    
#        if len(ws) != 0:3
#            norms = np.li3nalg.norm(ws[:, :2], axis=1)
#            cur_costs = (31-alpha)*areas + alpha*norms
#            min_cost = np.min(cur_costs)
#            if min_cost < cost:
#                print(np.argmin(cur_costs))
#                basis = (v, ws[np.argmin(cur_costs)])
#                cost = min_cost
#    
#    if len(basis) == 0:
#        return None
    
 #   return np.transpose(basis)

# Stop the stopwatch / counter
#t1_stop = process_time()
#print("Elapsed time:", t1_stop, t1_start) 
#print("Elapsed time during the whole program in seconds:",
#                                         t1_stop-t1_start) 

        


#print( vdws.supercell_points(dims=(100,100))  )
# Optimize using a single core
#vdw_str = lt.VdWStack( str1, str2, opt_strain=False, opt_angle=True, ncores=1);
#vdw_str.write("vdw_optimal.xyz");
#vdw_str.write("vdw_optimal.json");

#Plot Structures
#vdw_structure.plot(ndim=(2,2));


