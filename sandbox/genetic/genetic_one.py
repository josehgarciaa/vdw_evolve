"""
Base class for a genetic algorithm that works with fitness functions that get as parameters vectors or integers.

##example:

def function_a(x):
    x=x[0]
    return x*x+3*x-5+np.sqrt(np.sqrt(x**2))+5*np.sin(x-5)*2*x

a = 50
X =[i for i in range(-a,a)]
Y =[function_a([x]) for x in X]

model_pars = {'cell_split_number': 5, 'subjects_in_cell': 2,
              'nr_clones': 5, 'mutation_gaussian_sd': 2,
              'pins': 9,'gene_quality': 1,}
input_size = 1 # since our function has one variable this is a trivial scenario.
bounds =[[-a,a]] # searching intervals for the solution

# Experiment
experiment = Gen1(function_a ,input_size, bounds, model_par)
nr_epochs = 20
last_generation = experiment.evolve(nr_epochs)
fit = [function_a(x) for x in last_generation]
plt.scatter(last_generation, fit , c='b', label='generation', alpha=0.3)
plt.plot(X,Y)
plt.legend()

plt.plot(X,Y)

##
"""

import random
import numpy as np
from scipy.stats import truncnorm


def dummy_condition(generation):
    """
    Dummy condition.
    :return:: False
    """
    return False


class Gen1:

    def __init__(self, fitness_function, input_size, input_ranges, model_par):
        """
        Base class for genetic algorithm,
        this  example is constrained to optimize over integers.
        :param fitness_function: The function of which the minimum value needs to be found. .
        :param input_size::int Number of variables of the fitness function.
        :param input_ranges:[[int, int]...] Boundaries for each parameter of the fitness_function
        :param model_par:{} A python dictionary that contains the hyper-parameters of the experiment.

        """
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

            for subject_number in range(self.model_par["subjects_in_cell"]):
                subject = []
                for feature_bound in cell_constrains:
                    subject.append(random.randint(feature_bound[0], feature_bound[1]))
                generation_0.append(subject)

        return generation_0

    def mitosis(self):
        """
        This function will clone the individuals from a generation and add a mutation according to the predefined
        hyper parameter "mutation_gaussian_sd".
        :return::[[[int,int,..],[int,int,..]..]...]
        """
        family_generation = []
        for subject in self.actual_generation:
            family = [subject]
            for clone_nr in range(self.model_par['nr_clones']):

                clone = []
                for i, feature in enumerate(subject):
                    # mutation
                    bound = self.input_ranges[i]
                    mutation = self.mutation(c_value=feature, bounds=bound)
                    clone.append(mutation)

                family.append(clone)

            family_generation.append(family)

        return family_generation

    def kill(self, family_generation):
        """
        Surviving of the fittest, this function kill will keep only the best from each family and will return
        survivors and the performance of them.
        :param family_generation:: [[[int,int,..],[int,int,..]..]...] generated in th mitosis step
        :return:: [[int, int ,...],[int,..]], [float, float...]
        """
        generation = []
        performances = []
        for family in family_generation:
            champ = family[0]
            champ_performance = -self.ev_function(family[0])
            for children_nr in range(len(family)):
                performance = -self.ev_function(family[children_nr])
                if performance > champ_performance:
                    champ_performance = performance
                    champ = family[children_nr]
            generation.append(champ)
            performances.append(champ_performance)

        return generation, performances

    def reproduction(self, generation, performances):
        """
        This function mixed the solutions in order to give birth to new possible solutions.
        It requires as an input an old generation, and the output will be a new generation made by mixing the previous
        solutions. Hyper parameters that affect this function are 'pins' and 'gene_quality'.

        :param generation:: [[int, int ,...],[int,..]]
        :param performances:: [float, float...]
        :return::[[int, int ,...],[int,..]]
        """
        nr_pins = self.model_par["pins"]
        lower_pick = np.sqrt(min(performances) * min(performances))
        performances = [p + lower_pick for p in performances]
        s = 0
        for performance in performances:
            s += performance

        pi_portion_sizes = [performance / s for performance in performances]
        pi_portions = []
        p = 0
        for portion in pi_portion_sizes:
            p = p + portion
            pi_portions.append(p)

        random_displacement = random.uniform(0, 1)  # the spinning of the wheel
        pins = [np.sin(random_displacement + 2 * np.pi / nr_pins * i) for i in range(nr_pins)]

        parents = []
        for pin in pins:
            for i, candidate in enumerate(pi_portions):
                if pin <= candidate:
                    parents.append(generation[i])
                    break

        new_generation = []
        for i, parent in enumerate(parents):
            for j in range(i + 1, len(parents)):
                parent_1 = parent
                parent_2 = parents[j]
                new_generation.append(self.mating(parent_1, parent_2))

        return new_generation

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
        return children

    def evolve(self, nr_epochs, condition=dummy_condition, initial_generation=0, history=False):
        """
        This function evolves the experiment until the nr of evolution [nr_epochs]
        is reached or until a specific condition is fulfilled.
        it returns the last generation.

        :param nr_epochs: int
        :param condition: function with a boolean output
        :param initial_generation:: [[int, int ,...],[int,..]] or 0 if no generation is provided
        :param history: bool if true the experiment history is return.
        :return: [[int, int ,...],[int,..]]
        """
        if initial_generation == 0:
            self.actual_generation = self.first_generation()

        new_generation = self.actual_generation

        for e in range(nr_epochs):
            print("Epochs {}/{}".format(e, nr_epochs))
            if history:
                evolution_history = {}

            # mutate
            gen_0 = self.actual_generation
            clone_family = self.mitosis()
            # kill
            generation, performances = self.kill(clone_family)
            # new born
            new_generation = self.reproduction(generation, performances)
            # update
            self.actual_generation = new_generation
            if history:
                epoch_status = {'gen_0': gen_0,
                                'clone_family': clone_family,
                                'generation': generation,
                                'performances': performances,
                                'new_generation': new_generation}
                evolution_history[e] = epoch_status

            # update
            self.actual_generation = new_generation
            # actualisation


            if condition(new_generation):
                print("Condition fulfilled. Evolution process ended!\n Last generation{}".format(new_generation))
                if history:
                    return new_generation, evolution_history
                return new_generation

            print("Evolution process ended!\nLast generation: {}".format(new_generation))
        if history:
            return new_generation, evolution_history
        return new_generation

    # utility functions
    def _nr_cells(self):
        """
        Returns the number of cells (subspaces) in which the first solution attempts are generated.
        :return:: int
        """
        nc = 1
        for i in range(self.input_size):
            nc = nc * self.model_par["cell_split_number"]
        return nc

    def mutation(self, c_value, bounds):
        """
        Induce a small mutation in the parameter c_value by sampling from a gaussian distribution centered around him.
        Hyper parameter 'mutation_gaussian_sd' control the gaussian standard deviation.
        :param c_value:: int Current value.
        :param bounds:: [int int] The range of values in which this parameter must be.
        :return:: int Mutated value.
        """
        mean = c_value
        low = bounds[0]
        upp = bounds[1]
        sd = self.model_par['mutation_gaussian_sd']
        trunc_norm = truncnorm((low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd).rvs()
        return round(trunc_norm)
