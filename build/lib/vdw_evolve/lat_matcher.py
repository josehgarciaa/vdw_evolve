
import numpy as np
from scipy import optimize
from .geometry import ChangeBasis, MinimalBasis
import copy

class LatMatch:

    _angle_range = (-np.pi/2, np.pi/2)
    _strain_range = [(-0.05, 0.05), (-0.05, 0.05)]    
    _opt_angle = True
    _opt_strain = True
    _uni_strain = False
    max_dims = (100, 100)
    target = None
    reference = None
    cost = None

    def __init__(self, optimize_angle=True, optimize_strain=True):
        self._opt_angle = optimize_angle
        self._opt_strain = optimize_strain

    def opt_strain(self, x=True):
        self._opt_strain = x
        return self

    def opt_angle(self, x=True):
        self._opt_angle = x
        return self

    def uniform_strain(self, abs_maxstrain, format="perc"):
        self._uni_strain = True
        self.opt_strain()
        smax = abs_maxstrain
        scal = 1.0
        if format == "perc":
            scal = 1/100
        self._strain_range = [(scal*-smax, scal*smax)]
        print("strain range", self._strain_range)
        
        return self
   
    def set_max_strain(self, smax, format="perc"):
        self.opt_strain()
        s1, s2 = smax
        scal = 1.0
        if format == "perc":
            scal = 1/100
        self._strain_range = [(-scal*s1, scal*s1), (-scal*s2, scal*s2)]
        return self

    def set_angle_range(self, angle_range, format="deg"):
        self.opt_angle()
        amin, amax = angle_range
        scal = 1.0
        if format == "deg":
            scal = np.pi/180.
        self._angle_range = (scal*amin, scal*amax)
        return self

    def costFuncion(self, r, eta=0.001):
        dx, dy, dz = (r - np.floor(r))
        cost = 0
        etac = np.sqrt(2) * eta
        n = int(np.ceil(eta))
        for (nx, ny) in np.mgrid[-n:n, -n:n].T.reshape(n ** 2 * 4, 2):
            d2 = (dx + nx) ** 2 + (dy + ny) ** 2
            cost += np.mean(np.exp(- d2 / (2 * etac))) / n ** 2
        return -cost

    def strain_range(self, format="perc"):
        if format == "perc":
            return self._strain_range*100
        return self._strain_range

    def angle_range(self, format="deg"):
        if format == "deg":
            return self._angle_range*180/np.pi
        return self._angle_range
    
    def unpack_params(self, x):
        s1, s2, angle = 0, 0, 0

        if self._opt_strain:
            if self._uni_strain:
                s1 = float(x[0])
                s2 = s1
            else:
                s1, s2 = float(x[0]), float(x[1])
        if self._opt_angle:
            angle = float(x[-1])
        return s1, s2, angle
    
    def fitness(self, x):
        s1, s2, angle = self.unpack_params(x)
        self.opt = copy.copy(self.target).transform2D(strain=(s1, s2),
                                                      strain_format="abs",
                                                      angle=angle,
                                                      angle_format="rad",
                                                      transform_atoms=False)
        SCpoints = self.opt.supercell_points(self.max_dims)
        rBinA = ChangeBasis(SCpoints, self.reference.cell)
        return self.costFuncion(rBinA)

    def minimalCell(self, target, reference,
                    max_dims=(100, 100),  force=False):
        
        self.max_dims = max_dims
        self.target = target
        self.reference = reference

        bounds = []
        if self._opt_strain:
            bounds += self._strain_range
        if self._opt_angle:
            bounds += [self._angle_range]

        init_cost = self.fitness((0, 0, 0))

        res = optimize.differential_evolution(self.fitness, bounds,
                                              maxiter=1000, popsize=1000,
                                              polish=True)
        self.cost = res.fun
        print("For bounds", bounds, "Initial cost", init_cost, "Final cost", res.fun, "improvement", res.fun/init_cost )
        s1, s2, angle = self.unpack_params(res.x)
        self.opt = copy.copy(target).transform2D(strain=(s1, s2),
                                                 strain_format="abs",
                                                 angle=angle,
                                                 angle_format="rad")

        basis = MinimalBasis(self.opt.supercell_points(self.max_dims))
        if basis is None:
            return None
        return (((s1, s2), angle), basis, self.opt)
