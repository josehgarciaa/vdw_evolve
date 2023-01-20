"""
This file wraps up the routine for calculating the minimum supercell for overlapping the 2D lattices.
Example:

"""
import numpy as np
from .anneling import Annealing1
from .utils import t_cel1t_cel2, brut_strain, strained_proces


def super_cell(cel1, cel2, nr_epochs, model_par, fit_function):  # cost function as parameter
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

    t_cel1, strain_t_cel2 = t_cel1t_cel2(solution, cel1, cel2)
    strain = brut_strain(strain_t_cel2)
    diagonal_strain = strained_proces(strain_t_cel2, strain_boundary=model_par["strain_boundary"])
    t_cel2 = np.dot(np.linalg.inv(strain), strain_t_cel2)

    return t_cel1, strain_t_cel2, t_cel2, diagonal_strain, strain,
