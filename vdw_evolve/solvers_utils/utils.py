from scipy.optimize import minimize

import numpy as np
from numpy.linalg import inv


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


def strain_tune(x, strain_tb, optimize=True):
    """
    builds a diagonal matrix with the X parameters that will have the  job of strain
    :param x:
    :param strain_tb:
    :param optimize:
    :return:
    """
    strain = np.array([[1 + x[0], 0],
                       [0, 1 + x[1]]])

    # print("000_st:\n",strain)
    if optimize:
        tb = np.dot(inv(strain), strain_tb)
        cons = 99999
        round_count = 0
        for row in tb:
            for e in row:
                round_count += (round(e) - e) * (round(e) - e)  # e*e
        round_count = round_count * cons + (x[1] ** 2 + x[0] ** 2) * cons
        return round_count

    # print("stt:", strain)
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
        # print("x", x)
        sc = strain_tune(x, strain_tb=target_matrix)
        return sc

    x0 = np.array([0, 0])

    res = minimize(strain_cost, x0,
                   # constraints =cons,
                   bounds=strain_boundary,
                   method='L-BFGS-B')

    strain = strain_tune(res.x, target_matrix, optimize=False)
    return strain


def brut_strain(strain_tb):
    """
    Will generate a matrix
    that multiplied with the target matrix round it to the closest integer matrix.
    :param strain_tb: [[float,float],[float,float]]
    :return: [[float,float],[float,float]]
    """

    # strain * tb = strain_tb
    # strain = strain_tb * (tb^(-1))
    # tb =

    round_matrix = strain_tb.copy()
    for i in range(len(strain_tb)):
        for j in range(len(strain_tb[0])):
            round_matrix[i][j] = round(strain_tb[i][j])
    strain = np.dot( strain_tb, inv(round_matrix))
    return strain


def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def rotation_10(v):
    th = angle_between(v, np.array([1,0]))
    rotation = np.array( [[np.cos(-th), -np.sin(-th)],
                          [np.sin(-th), np.cos(-th)]])
    return rotation

def allign_along_10(celss):

    c = []
    for cel in celss:
        ro = rotation_10(np.array([cel[0][0], cel[1][0]]))
        c.append(np.dot(ro, cel))

    return c



def subject_test(params, cel1,cel2):
    tA, tB= t_cel1t_cel2(params,cel1,cel2)
    if np.linalg.det(tA)==0:
        return False
    if np.linalg.det(tB)==0:
        return False
    return True


