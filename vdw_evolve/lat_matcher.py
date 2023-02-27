
import numpy as np
from scipy import optimize
from .geometry import ChangeBasis, MinimalBasis
import copy

class LatMatch:
    opt_angle = True
    opt_strain = True
    uni_strain = False
    angle_range = (-np.pi / 2, np.pi / 2)
    strain_range = [(-0.05, 0.05), (-0.05, 0.05)]
    max_dims = (100,100)
    target = None
    reference = None
    
    def __init__(self, optimize_angle=True, optimize_strain=True):
        self.opt_angle = optimize_angle
        self.opt_strain = optimize_strain

    def OptStrain(self, x=True):
        self.opt_strain = x
        return self

    def OptAngle(self, x=True):
        self.opt_angle = x
        return self

    def UniformStrain(self, smax, smin=None):
        if smin is None:
            smin = -smax
        self.uni_strain = True
        self.OptStrain()
        self.strain_range = [(smin, smax)]
        return self
   
    def Strain(self, smax):
        self.OptStrain()
        self.max = smax
        return self

    def AngleRange(self, angle_range):
        self.OptAngle()
        self.angle_range = angle_range
        return None

    def costFuncion(self, r, eta=0.001):
        dx, dy, dz = (r - np.floor(r))
        cost = 0
        etac = np.sqrt(2) * eta
        n = int(np.ceil(eta))
        for (nx, ny) in np.mgrid[-n:n, -n:n].T.reshape(n ** 2 * 4, 2):
            d2 = (dx + nx) ** 2 + (dy + ny) ** 2
            cost += np.mean(np.exp(- d2 / (2 * etac))) / n ** 2
        return -cost

    def unpack_params(self, x):
        s1, s2, angle = 0, 0, 0

        if self.opt_strain:
            if self.uni_strain:
                s1 = float(x[0])
                s2 = s1
            else:
                s1, s2 = float(x[0]), float(x[1])
        if self.opt_angle:
            angle = float(x[-1])
        return s1, s2, angle
    
    def fitness(self, x):
        s1, s2, angle = self.unpack_params(x)
        self.opt = copy.copy(self.target).transform2D((s1, s2), angle, transform_atoms=False)
        SCpoints = self.opt.supercell_points(self.max_dims)
        rBinA = ChangeBasis(SCpoints, self.reference.cell)
        return self.costFuncion(rBinA)

    def minimalCell(self, target, reference, max_dims=(100, 100),  force=False):
        self.max_dims = max_dims
        self.target = target
        self.reference = reference
        bounds = []
        if self.opt_strain:
            bounds += self.strain_range
        if self.opt_angle:
            bounds += [self.angle_range]
        res = optimize.differential_evolution(self.fitness, polish=True)
        s1, s2, angle = self.unpack_params(res.x)
        self.opt = copy.copy(target).transform2D((s1, s2), angle)
        basis = MinimalBasis(self.opt.supercell_points(self.max_dims))
        if basis is None:
            return None
        return ( ((s1,s2),angle), basis, self.opt)