import numpy as np
from numpy.linalg import det, inv


from anneling import Annealing1


def build_annealing(cel1, cel2, nr_epochs, model_par):
    """

    :param model_par:
    :param cel1:
    :param cel2:
    :param nr_epochs:
    :return:
    """
    fit = lambda params: fit_function(params, a=cel1, b=cel2)
    experiment = Annealing1(fit, model_par["start_point"], model_par)
    history_book = experiment.evolve(nr_epochs, prints_p=5)

    solution = experiment.actual_solution
    tA, tB = tAtB(solution, cel1, cel2)
    strain = get_strain(tB)
    return tA, tB, strain


def at_sin(x, up=99999999):
    """

    :param x:
    :param up:
    :return:
    """
    tr_x = (x + 0.5)

    if x < 1:
        res = up * (1 / (1 + np.sqrt((x ** 2)))) * np.sin(tr_x * np.pi) + 1 / (x + 0.00000000000001)
    else:
        res = (x - 1) ** 2 - up * (1 / (1 + x ** 2)) + 1 / (x + 0.00000000000001)
    return res


def tAtB(params, a, b):
    tA = np.array([[params[0], params[1]],
                   [params[2], params[3]]])

    tB = np.dot(np.dot(tA, a), inv(b))  # tAa=tBb

    return tA, tB


def fit_function(params, a, b):
    tA, tB = tAtB(params, a, b)
    tAa = np.dot(tA, a)
    tBb = np.dot(tB, b)

    # main condition
    zero_mat = tAa - tBb
    s = 0
    for row in zero_mat:
        for e in row:
            s += e * e

    # mimimum TA
    detTAa = det(tAa) * det(tAa)  # minimum but biger than 0
    detTBb = det(tBb) * det(tBb)

    # TB integer
    cons = 9999999
    tB_con = 0
    for row in tB:
        for e in row:
            tB_con += ((round(e) - e)) * ((round(e) - e))  # e*e
    tB_con = tB_con * cons

    tA_lenghth = ((tA[0][0] ** 2 + tA[0][1] ** 2) + (tA[1][0] ** 2 + tA[1][1] ** 2)) * 100

    f = at_sin(detTAa) + tB_con ** 2 + tA_lenghth
    # ((1-detTAa)**2)*k_p + tB_con**2

    return f


def get_strain(tB):
    tBr = tB.copy()
    for i in range(len(tB)):
        for j in range(len(tB[0])):
            tBr[i][j] = round(tB[i][j])
    S = np.dot(tBr, inv(tB))
    return (S)
