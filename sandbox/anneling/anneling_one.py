"""
Base class for a anneling optimisation algorithm that works with fitness functions that get as parameters vectors or integers.

##example:
def function_a(x):
    x=x[0]
    return x*x+3*x-5+np.sqrt(np.sqrt(x**2))+5*np.sin(x-5)*2*x
# Experiment settings

#Hyperparameters of the model
model_par = {

    'initialTemp': 1,
    'finalTemp': 0.0003,

    'beta': 10,
    'bounds': [[-a,a]],

    'nr_neighbours': 5,
    'step_size': 4,
    'gaussian_sd':5,

    'known_min': -100

}

#
input_size = 1 # since our function has one variable this is a trivial scenario.
start_point = [40]

# Experiment
experiment = Annealing1(function_a ,start_point, model_par)

epochs= 5
history_book = experiment.evolve( epochs, prints_p=5)
"""

import random
import numpy as np
from scipy.stats import truncnorm


class Annealing1:

    def __init__(self, fitness_function, start_point, model_par, history=True):
        """
        base class for annealing algorithm,
        :param fitness_function: The function of which the minimum value needs to be found.
        :param start_point: [int , int..] starting parameters for solution search.
        :param model_par: {} A python dictionary that contains the hyper-parameters of the experiment.
        :param history: True/False if true the evolution of the model will return the searching steps as a dictionary.
        """
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
        """
        Update the temperature of the system at each step.
        The speed of temperature decrees is dictated by the model parameter:"beta"
        """
        self.actual_temp = self.actual_temp / (1 + self.model_par["beta"] * self.actual_temp)

    def stop_criteria(self):
        """

        :return: True  if the system satisfy the stopping criteria False otherwise.
        """
        flag = False
        flag = self.actual_temp <= self.final_temp
        if self.actual_value <= self.model_par['known_min']:
            flag = True

        return flag

    def get_neighbours(self):
        """
        Returns a list of neighbouring solutions of the length model_par["nr_neighbours"].
        :return:[[int, int..],[int, int ..]] neighbours of actual solution
        """
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
        """
        This function evolves the experiment until the nr of evolution [nr_epochs]
        is reached or until a specific condition is fulfilled.
        it returns the history os the states.

        :param epochs: int nr. of epochs spent at etch temperature.
        :param prints_p: int how often the state of the system should be printed.
        :param tr_print: if False the system informations are not printed during the evolution.
        :return: {} history of the system states during evolution
        """
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
        Generates an integer from a gaussian distribution around the value c_value
        but in the bounded range.
        :param c_value: float
        :param bounds: [flaot, float]
        :return: float
        """
        mean = c_value
        low = bounds[0]
        upp = bounds[1]
        sd = self.model_par['gaussian_sd']
        # print("p1:{}|p2:{}|mean:{}|sd:{}".format((low - mean) / sd, (upp - mean) / sd, mean,sd))
        trunc_norm = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()
        return round(trunc_norm)
