"""
** file:///C:/Users/tomut/Downloads/neighbour%20simulated%20anneling.pdf
"""

import random
import numpy as np
from scipy.stats import truncnorm

from .annealing_one import Annealing1


class GradientAnnealing(Annealing1):

    def get_neighbours(self):
        delta_p = self.model_par["delta_p"]
        gradient_step_size = self.model_par["g_step_size"]
        step_size = self.model_par["step_size"]
        influence = self.model_par["gradient_influence"]

        neighbours = []

        for _ in range(self.model_par["nr_neighbours"]):
            neighbour = []

            flag = random.uniform(0, 1)
            if flag <= influence:
                grad_flag = True
            else:
                grad_flag = False

            if grad_flag:
                for p_in in range(self.nr_parameters):
                    p_in_m = self.actual_solution.copy()
                    p_in_m[p_in] -= delta_p
                    p_in_p = self.actual_solution.copy()
                    p_in_p[p_in] += delta_p
                    grad = (self.ev_function(p_in_m) - self.ev_function(p_in_p)) / (2 * delta_p)
                    p = self.actual_solution[p_in] + gradient_step_size * grad  # maybe heare i could add the gradient
                    # print("p:",p)
                    # print("bounds:",self.bounds[p_in] )
                    p = self.g_distribution(round(p), self.bounds[p_in])
                    p = round(p)
                    neighbour.append(p)
            else:
                for p_in in range(self.nr_parameters):
                    mi_plu = random.randint(0, 1)
                    if mi_plu == 0:
                        mi_plu = -1
                    p = self.actual_solution[p_in] + mi_plu * step_size
                    p = self.g_distribution(p, self.bounds[p_in])
                    p = round(p)
                    neighbour.append(p)

            neighbours.append(neighbour)

        return neighbours
