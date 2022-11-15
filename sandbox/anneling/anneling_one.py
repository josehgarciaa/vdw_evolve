"""

"""

import random
import numpy as np
from scipy.stats import truncnorm


class Annealing1:

    def __init__(self, fitness_function, start_point, model_par, history=True):

        self.ev_function = fitness_function
        self.model_par = model_par
        self.actual_solution = start_point
        self.actual_value = self.ev_function(self.actual_solution)
        self.actual_temp = self.model_par["initialTemp"]
        self.final_temp = self.model_par["finalTemp"]
        self.history = history
        self.bounds = self.model_par["bounds"]
        self.nr_parameters = len(start_point)

    def temp_reduction(self):
        self.actual_temp = self.actual_temp / (1 + self.model_par["beta"] * self.actual_temp)

    def stop_criteria(self):
        flag = False
        flag = self.actual_temp <= self.final_temp
        if self.actual_value <= self.model_par['known_min']:
            flag = True

        return flag

    def get_neighbours(self):
        neighbours = []

        for _ in range(self.model_par["nr_neighbours"]):
            neighbour = []
            for p_in in range(self.nr_parameters):
                mi_plu = random.randint(0, 1)
                if mi_plu == 0:
                    mi_plu = -1
                p = self.actual_solution[p_in] + mi_plu * self.model_par[
                    "step_size"]  # maybe heare i could add the gradient
                p = self.g_distribution(p, self.bounds[p_in])
                p = round(p)

                neighbour.append(p)

            neighbours.append(neighbour)

        return neighbours

    def evolve(self, epochs, prints_p=5, tr_print=True):
        history_book = {'solutions': [self.actual_solution],'values':[self.actual_value], 'changes': [0], 'temperature': [self.actual_temp]}
        while self.stop_criteria() == False:

            for epoch in range(epochs):

                # generate neighbours
                neighbours = self.get_neighbours()

                # neighbour picking
                choice = np.random.choice([i for i in range(len(neighbours))])
                new_solution = neighbours[choice]
                new_val = self.ev_function(new_solution)
                # change in solution
                change = -self.ev_function(self.actual_solution) + new_val

                if change <= 0:
                    self.actual_solution = new_solution
                    self.actual_value = new_val
                else:
                    if random.uniform(0, 1) < np.exp(-change / self.actual_temp):
                        self.actual_solution = new_solution
                        self.actual_value = new_val

                val = self.ev_function(self.actual_solution)
                if epoch % prints_p == 0 and tr_print:
                    print("temp:{}|epoch:{}|change:{}|value:{}".format(self.actual_temp, epoch, change, val, ))
            # temperature adjustment
            self.temp_reduction()
            if self.history:
                history_book['solutions'].append(self.actual_solution)
                history_book['values'].append(self.actual_value)
                history_book['changes'].append(change)
                history_book['temperature'].append(self.actual_temp)

        return history_book

    def g_distribution(self, c_value, bounds):
        """

        """
        mean = c_value
        low = bounds[0]
        upp = bounds[1]
        sd = self.model_par['gaussian_sd']
        # print("p1:{}|p2:{}|mean:{}|sd:{}".format((low - mean) / sd, (upp - mean) / sd, mean,sd))
        trunc_norm = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()
        return round(trunc_norm)
