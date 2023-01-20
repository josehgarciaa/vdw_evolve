"""

"""
import matplotlib.pyplot as plt
import numpy as np

def plot_unit_cell_vectors(ax, cell,label=None, alpha=0.3,color='b'):

    # f, ax = plt.subplots()
    if label!=None:
        a1 = ax.arrow(0, 0, cell[0][0], cell[1][0], head_width=0.05, head_length=0.5, color=color, label="a1_"+label, alpha=alpha)
        a2 = ax.arrow(0, 0, cell[0][1], cell[1][1], head_width=0.05, head_length=0.5, color=color, label="a2_"+label, alpha=alpha)
    else :
        a1 = ax.arrow(0, 0, cell[0][0], cell[1][0], head_width=0.05, head_length=0.5, color=color,
                      alpha=alpha)
        a2 = ax.arrow(0, 0, cell[0][1], cell[1][1], head_width=0.05, head_length=0.5, color=color,
                      alpha=alpha)
    return a1,a2


def plot_lattice_grid(ax, s_lattice,parent, label, color, alpha =0.2):
    a1s,a2s = s_lattice.get_scales(s_lattice.cells[parent])
    cell_x = s_lattice.cells[parent]
    # f, ax = plt.subplots()
    plots = []
    plots.append(ax.plot([-a1s * cell_x[0][0], a1s * cell_x[0][0]], [-a1s * cell_x[1][0], a1s * cell_x[1][0]], color=color, label=label,
                         alpha=alpha))
    plots.append(ax.plot([-a2s * cell_x[0][1], a2s * cell_x[0][1]], [-a2s * cell_x[1][1], a2s * cell_x[1][1]], color=color, alpha=alpha))

    for k in range(-a1s, a1s):
        plots.append(ax.plot([-a2s * cell_x[0][1] + k * cell_x[0][0], a2s * cell_x[0][1] + k * cell_x[0][0]],
                             [-a2s * cell_x[1][1] + k * cell_x[1][0], a2s * cell_x[1][1] + k * cell_x[1][0]],
                             color=color,
                             alpha=alpha))
    for k in range(-a2s, a2s):
        plots.append(ax.plot([-a1s * cell_x[0][0] + k * cell_x[0][1], a1s * cell_x[0][0] + k * cell_x[0][1]],
                             [-a1s * cell_x[1][0] + k * cell_x[1][1], a1s * cell_x[1][0] + k * cell_x[1][1]],
                             color=color, alpha=alpha))
    return plots


def plot_atom_grid(ax,s_lattice,parent, color, alpha =0.2):

    lattice = s_lattice.lattices[parent]
    cel0 = s_lattice.cells[parent]
    a1s,a2s = s_lattice.get_scales(s_lattice.cells[parent])

    at =[]

    for e, atom in enumerate(lattice):
        atom_x = lattice[atom][0]
        atom_y = lattice[atom][1]
        x = []
        y = []
        for i in range(-a1s, a1s):
            for j in range(-a2s, a2s):
                x.append(atom_x + i * cel0.T[0][0] + j * cel0.T[1][0])
                y.append(atom_y + i * cel0.T[0][1] + j * cel0.T[1][1])
        at.append(ax.scatter(x, y, color=color, alpha=alpha))
    return at

def plot_cell_atoms(cell_atoms,color="b",marker="*",alpha=0.3,):
    x =[]
    y = []
    for atom in cell_atoms:
        x.append(cell_atoms[atom][0])
        y.append(cell_atoms[atom][1])
        plt.scatter(x,y, marker=marker,color=color,alpha = alpha)

def plot_cell(ax,cell,label="",color='r'):
    c0 = np.array([0,0])
    c1 = np.array([cell[0][0], cell[1][0]])
    c2 = np.array([cell[0][1], cell[1][1]])
    c3 = c1+c2
    c=[]
    if label!="":
        c.append(ax.plot([c0[0], c1[0]], [c0[1], c1[1]], color=color,label=label, alpha =0.2))
    else:
        c.append(ax.plot([c0[0], c1[0]], [c0[1], c1[1]], color=color,  alpha=0.2))
    c.append(ax.plot([c0[0], c2[0]], [c0[1], c2[1]], color=color, alpha =0.2))
    c.append(ax.plot([c1[0], c3[0]], [c1[1], c3[1]], color=color, alpha =0.2))
    c.append(ax.plot([c2[0], c3[0]], [c2[1], c3[1]], color=color, alpha =0.2))
    return c


