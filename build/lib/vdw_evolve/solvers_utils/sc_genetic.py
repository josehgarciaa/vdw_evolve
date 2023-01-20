"""

"""
import numpy as np
import random

from .genetic import Gen1
from .utils import subject_test, t_cel1t_cel2,brut_strain,strained_proces


# Centered genetic -> Gen1 adapted for this scenario
class MinimalSupercel_generation(Gen1):

    def __init__(self,cells, fitness_function, input_ranges, model_par,input_size=4):
        self.cel1= cells[0]
        self.cel2=cells[1]
        self.ev_function = fitness_function
        self.input_size = input_size
        self.input_ranges = input_ranges
        self.model_par = model_par
        self.actual_generation = []



    def first_generation(self):
        """
        This function will generate a random generation of possible solutions based on the model hyper parameters
        :return::[[int, int ,...],[int,..]]
        """
        generation_0 = []

        number_of_cells = self._nr_cells()

        for cell_nr in range(number_of_cells):
            cell_constrains = []

            for f in range(self.input_size):
                feature = self.input_size - f - 1
                cell_coordinate = int(cell_nr / (self.model_par["cell_split_number"] ** feature))
                cell_nr = cell_nr - cell_coordinate * int((self.model_par["cell_split_number"] ** feature))

                interval = (-self.input_ranges[feature][0] + self.input_ranges[feature][1]) / self.model_par[
                    'cell_split_number']
                c_0 = self.input_ranges[feature][0] + interval * cell_coordinate
                c_1 = self.input_ranges[feature][0] + interval * (cell_coordinate + 1)
                constrain = [c_0, c_1]
                cell_constrains.append(constrain)

            # 0_in_cell_enhancement
            subjects_in_cell = self.model_par["subjects_in_cell"]
            for feature_bound in cell_constrains:
                if feature_bound[0] <= 0 <= feature_bound[0]:
                    subjects_in_cell = subjects_in_cell * self.model_par["0_in_cell_enhancement"]

            for subject_number in range(subjects_in_cell):
                subject = []
                for feature_bound in cell_constrains:
                    new_subject = random.randint(int(feature_bound[0]), int(feature_bound[1]))
                    subject.append(new_subject)

                while not subject_test(subject, self.cel1, self.cel2):
                    subject = []
                    for feature_bound in cell_constrains:
                        new_subject = random.randint(int(feature_bound[0]), int(feature_bound[1]))
                        subject.append(new_subject)

                generation_0.append(subject)

        return generation_0

    def mating(self, parent_1, parent_2):
        """

        :param parent_1:: [int, int, ..]
        :param parent_2:: [int, int, ..]
        :return:: [int, int, ..]
        """
        children = []
        for char_nr in range(len(parent_1)):
            new_char = round((parent_1[char_nr] + parent_2[char_nr]) / 2)
            children.append(new_char)

        while not subject_test(children,self.cel1, self.cel2):
            new_children = []
            for i in range(len(children)):
                new_children.append(self.mutation(children[i], bounds=self.input_ranges[i]))

            children = new_children
        return children

###

def genetic_cell(cel1, cel2, nr_epochs, model_par, fit_function):  # cost function as parameter
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
    experiment = MinimalSupercel_generation([cel1,cel2],cost , model_par["bounds"], model_par)


    # Evolve the experiment to found the solution.
    last_generation = experiment.evolve(nr_epochs)
    fit = [cost(x) for x in last_generation]
    best = last_generation[fit.index(min(fit))]
    t_cel1, strain_t_cel2 = t_cel1t_cel2(best, cel1, cel2)

    strain = brut_strain(strain_t_cel2)
    diagonal_strain = strained_proces(strain_t_cel2, strain_boundary=model_par["strain_boundary"])
    t_cel2 = np.dot(np.linalg.inv(strain), strain_t_cel2)

    return t_cel1, strain_t_cel2, t_cel2, diagonal_strain, strain,

