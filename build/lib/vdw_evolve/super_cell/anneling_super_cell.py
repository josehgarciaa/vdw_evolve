"""
This file wraps up the routine for calculating the minimum supercell for overlapping the 2D lattices.
Example:

"""

import numpy as np
from numpy.linalg import det, inv

from scipy.optimize import minimize


from ..anneling import Annealing1


def super_cell(cel1, cel2, nr_epochs, model_par): # cost function as parameter
    # extra file with cost functions -> fit function not cost
    """

    :param model_par:
    :param cel1:
    :param cel2:
    :param nr_epochs:
    :return:
    """

    def cost(params):
        """
        Builds the function that needs to be optimized for the given case.
        :param params:[int,init,int,int] parameters  (tA =params.transpose(2,2))
        :return: float
        """
        return fit_function(params, cel1=cel1, cel2=cel2, strain_boundary=model_par["strain_boundary"])

    # Build the experiment setup.
    experiment = Annealing1(cost, model_par["start_point"], model_par)
    # Evolve the experiment to found the solution.
    history_book = experiment.evolve(nr_epochs, prints_p=5)
    solution = experiment.actual_solution

    t_cel1, t_cel2 = t_cel1t_cel2(solution, cel1, cel2)
    strain = brut_strain(t_cel2)
    diagonal_strain = strained_proces(t_cel2, strain_boundary=model_par["strain_boundary"])
    t_cel2_no_strain = np.dot(inv(strain), t_cel2)

    return t_cel1, t_cel2, t_cel2_no_strain, diagonal_strain, strain,


def t_cel1t_cel2(params, cel1, cel2):
    """

    :param params: [int,int,int, int] parameters for the transformation that will apply on cel1
    :param cel1:
    :param cel2:
    :return:
    """
    t_cel1 = np.array([[params[0], params[1]],
                       [params[2], params[3]]])

    t_cel2 = np.dot(np.dot(t_cel1, cel1), inv(cel2))  # tAa=tBb

    return t_cel1, t_cel2


def fit_function(params, cel1, cel2, strain_boundary):
    """
    This creates a cost that the annealing process will minimize.
    The cost is proportional to the area of the supercell
    and the length of the new based vectors.
    :param params: 
    :param cel1: 
    :param cel2: 
    :param strain_boundary: 
    :return:
    """
    t_cel1, t_cel2 = t_cel1t_cel2(params, cel1, cel2)
    super_c_1 = np.dot(t_cel1, cel1)

    # It's a chance that a solution with t_cel2 integer does not exist. So if the t_cel2 has float values,
    # we will decompose it in a new t_cel2 integer and a strain(diagonal matrix), which multiply and reconstruct
    # the initial matrix.
    strain = strained_proces(t_cel2, strain_boundary)
    t_cel2 = np.dot(inv(strain), t_cel2)

    # area of the new cel area = det(cel1)*det(t_cel1)
    super_cel_area = det(super_c_1) * det(super_c_1)  # minimum but bigger than 0

    # t_cel1 integer # cost for t_cel1 not be an integer
    cons = 99999
    round_cost = 0
    for row in t_cel2:
        for e in row:
            z_ero1 = (round(e) - e) * (round(e) - e)
            round_cost += z_ero1
    round_cost = round_cost * cons

    # cost for minimizing the new vectors length
    length_cost = ((t_cel1[0][0] ** 2 + t_cel1[0][1] ** 2) + (t_cel1[1][0] ** 2 + t_cel1[1][1] ** 2)) * 100

    f = partial_cost_shape(super_cel_area) + round_cost ** 3 + length_cost*2

    return f


def partial_cost_shape(x, up=99999999):
    """
    This function will take the area of the new cell as a parameter
    it is built in such that it will take the minimum value when
    the area is close to 1 and will go to inf in 0 and infinity.

    :param x: float  super_cel_area
    :param up: float  modul value of the peek in 1.
    :return: float
    """
    tr_x = (x + 0.5)

    if x < 1:
        res = up * (1 / (1 + np.sqrt((x ** 2)))) * np.sin(tr_x * np.pi) + 1 / (x + 0.000000000000001)
    else:
        res = (x - 1) ** 2 - up * (1 / (1 + x ** 2)) + 1 / (x + 0.00000000000001)
    return res


def strain_tune(x, target_matrix, optimize=True):
    """
    builds a diagonal matrix with the X parameters that will have the  job of strain
    :param x:
    :param target_matrix:
    :param optimize:
    :return:
    """
    strain = np.array([[1 + x[0], 0],
                       [0, 1 + x[1]]])

    # Strain*tB_strained =tB_ => tB_strained=inv(Strain)*tB
    if optimize:
        cell_strain = np.dot(inv(strain), target_matrix)
        cons = 99999
        round_count = 0
        for row in cell_strain:
            for e in row:
                round_count += (round(e) - e) * (round(e) - e)  # e*e
        round_count = round_count * cons + (x[1] ** 2 + x[0] ** 2) * cons
        return round_count

    return strain


def strained_proces(target_matrix, strain_boundary):
    """
    Return a strain matrix with the propriety that multiplying  the inverse with the target matrix
    will get a matrix of integers.

    :param target_matrix:
    :param strain_boundary:
    :return:
    """

    def strain_cost(x):
        return strain_tune(x, target_matrix=target_matrix)

    x0 = np.array([0, 0])
    res = minimize(strain_cost, x0,
                   # constraints =cons,
                   bounds=strain_boundary,
                   method='L-BFGS-B')

    strain = strain_tune(res.x, target_matrix, optimize=False)
    return strain


def brut_strain(initial_matrix):
    """
    Will generate a matrix
    that multiplied with the target matrix round it to the closest integer matrix.
    :param initial_matrix: [[float,float],[float,float]]
    :return: [[float,float],[float,float]]
    """
    round_matrix = initial_matrix.copy()
    for i in range(len(initial_matrix)):
        for j in range(len(initial_matrix[0])):
            round_matrix[i][j] = round(initial_matrix[i][j])
    strain = np.dot(round_matrix, inv(initial_matrix))
    return strain
