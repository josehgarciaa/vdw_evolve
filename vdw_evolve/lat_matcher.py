
import numpy as np
from scipy import optimize
from .geometry import ChangeBasis, MinimalBasis
import copy

class LatMatch:
    opt_angle = True
    opt_strain = True
    uni_strain = False
    angle_range = (-np.pi / 2, np.pi / 2)
    smax = [0.05, 0.05]
    bounds = None
    result = None

    def __init__(self, optimize_angle=True, optimize_strain=True):
        self.ref = None
        self.tar = None
        self.opt = None
        self.dim = None
        self.updated_target_cell = None
        self.sc_vec = None
        self.opt_angle = optimize_angle
        self.opt_strain = optimize_strain

    def OptStrain(self, x=True):
        self.opt_strain = x
        return self

    def OptAngle(self, x=True):
        self.opt_angle = x
        return self

    def UniformStrain(self, smax):
        self.uni_strain = True
        self.OptStrain()
        self.smax = smax
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
        dx, dy = (r - np.floor(r))
        cost = 0
        etac = np.sqrt(2) * eta
        n = int(np.ceil(eta))
        for (nx, ny) in np.mgrid[-n:n, -n:n].T.reshape(n ** 2 * 4, 2):
            d2 = (dx + nx) ** 2 + (dy + ny) ** 2
            cost += np.mean(np.exp(- d2 / (2 * etac))) / n ** 2
        return -cost

    def fitness(self, x):
        if self.uni_strain:
            s1 = x[0]
            s2 = x[0]
        else:
            s1, s2 = self.smax
        
        if (len(x) == 3):
            s1, s2, theta = x
        if (len(x) == 2):
            s1, s2 = x
        if (len(x) == 1):
            if not self.opt_angle:
                s1 = s2 = float(x)
            else:
                theta = float(x)
        self.opt = copy.copy(self.tar).transform2D((s1, s2), theta)
        SCpoints = self.opt.supercell_points(self.dim)
        rBinA = ChangeBasis(SCpoints, self.ref.cell)
        return self.costFuncion(rBinA)

    def supercellVectors(self, force=False):
        if self.sc_vec is not None and not force:
            return self.sc_vec

        # Catch the proper variables. This should be properly improved
        if not self.opt_strain:
            self.bounds = [(self.theta_range[0], self.theta_range[1])]
        elif not self.opt_angle:
            smax = self.smax
            if (len(smax) == 1):
                self.bounds = [(-smax[0], smax[0])]
                print(self.bounds)
            else:
                self.bounds = [(-smax[0], smax[0]), (-smax[1], smax[1])]
        else:
            smax = self.smax;
            self.bounds = [(-smax[0], smax[0]), (-smax[1], smax[1]), (self.theta_range[0], self.theta_range[1])]

        print("cost without otimization", self.fitness([0, 0, 0]))

        res = optimize.differential_evolution(self.fitness, self.bounds, maxiter=1000, popsize=1000, polish=True)

        # Catch the proper variables. This should be properly improved
        s1, s2, theta = 0, 0, 0
        if not self.opt_strain:
            theta = float(res.x)
        elif not self.opt_angle:
            if (len(smax) == 1):
                s1, s2 = float(res.x), float(res.x)
            else:
                s1, s2 = res.x
        else:
            s1, s2, theta = res.x
        result = [s1, s2, theta]
        print(result)
        self.result = result
        print("cost after otimization", self.fitness([s1, s2, theta]))

        self.opt = copy.copy(self.tar).transform2D((s1, s2), theta)
        
        basis = MinimalBasis()
        if basis is None:
            return None;
        
        return ((s1,s2),theta), MinimalBasis(self.dim, self.updated_target_cell, self.ref)

    def supercell(self, force=False):
        L = self.supercellVectors( force=force ).T
        N = np.linalg.norm(L,axis=1)
        idx = np.argsort(N)
        L = L[idx][1:]
        N = N[idx][1:]
        thetas = np.diag(1/N)@L
        thetas = thetas@(thetas.T)
        good_thetas = np.abs(thetas) < np.cos(self.theta_min)

        iop,jop, min = 0,0,1e14
        for i,row in enumerate(good_thetas):
            for j,good_angle in enumerate(row):
                pair_norm = N[i]**2+N[j]**2 + np.abs(np.cross(L[i],L[j]))
    #        pair_norm = np.abs(np.cross(L[i],L[j]));
                if (good_angle and pair_norm<min):
                  min =pair_norm
                  iop,jop = i,j
        sc_vec = np.transpose([ L[iop],L[jop] ])
        return sc_vec