import numpy as np
from .parser import extract_molecule, t_units_to_x


class SuperLattice:
    def __init__(self, super_cell, path_xyz1, path_xyz2):

        self.s_cell = super_cell.cell
        self.cells = super_cell.parents
        self.lattices = [t_units_to_x(extract_molecule(path_xyz1), self.cells[0]),
                         t_units_to_x(extract_molecule(path_xyz2), self.cells[1])]
        self.s_lattice = self.s_lattice_()
        self.nr_atoms = len(self.s_lattice)

    def s_lattice_(self):

        # print("\n\ns_cell",self.s_cell)
        # take atoms from 1st lattice:
        cell_atoms = {}
        for i, cell in enumerate(self.cells):
            a1s, a2s = self.get_scales(cell)
            atomic_grid = atom_grid(self.lattices[i], cell, a1s, a2s)


            for atom in atomic_grid:
                if in_cell_check(atomic_grid[atom], self.s_cell):
                    cell_atoms["l" + str(i) + "_" + atom] = atomic_grid[atom]
        # print("cell_atoms:\n",cell_atoms)
        return cell_atoms

    def get_scales(self, cell_x):
        s_cell = self.s_cell
        c0 = np.array([0, 0])
        c1 = np.array([s_cell[0][0], s_cell[1][0]])
        c2 = np.array([s_cell[0][1], s_cell[1][1]])
        c3 = c1 + c2

        c_ = [c0, c1, c2, c3]

        d = [np.sqrt(c[0] ** 2 + c[1] ** 2) for c in c_]

        a1s = int(max(d) / np.sqrt(cell_x.T[0][0] ** 2 + cell_x.T[0][1])) + 1
        a2s = int(max(d) / np.sqrt(cell_x.T[1][0] ** 2 + cell_x.T[1][1])) + 1

        return a1s, a2s


def atom_grid(lattice, cell_x, a1s, a2s):
    atoms = {}

    for e, atom in enumerate(lattice):
        atom_x = lattice[atom][0]
        atom_y = lattice[atom][1]
        counter = 1
        for i in range(-a1s, a1s):
            for j in range(-a2s, a2s):
                x = atom_x + i * cell_x.T[0][0] + j * cell_x.T[1][0]
                y = atom_y + i * cell_x.T[0][1] + j * cell_x.T[1][1]
                atoms[str(counter) + "_" + atom] = [x, y]
                counter += 1
    return atoms


# v = t1C1+t2C2
# v0 = t1c10 + t2c20
# v1 = t1c11 + t2c21
# v0*c11 = t1c10c11 +t2c20c11
# v1*c10 = t1c11c10 +t2c21c10
# v0*c11-v1*c10= t2(c20c11-t2c21c10) => t2= (v0*c11-v1*c10)/(c20c11-c21c10)

def in_cell_check(v, cell):
    c1 = cell.T[0]
    c2 = cell.T[1]
    t2 = (v[0] * c1[1] - v[1] * c1[0]) / (c2[0] * c1[1] - c2[1] * c1[0])
    t1 = (v[0] - t2 * c2[0]) / c1[0]

    if 0 < t1 <= 1 and 0 < t2 <= 1:
        return True
    else:
        return False
