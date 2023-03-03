import numpy as np
from vdw_evolve import AnnealingSolver, MechanicSolver, GeneticSolver


def run_experiment(user_request):
    if user_request["algo"] == "annealing":

        print("\n\nAnnealingSolver:")
        solver1 = AnnealingSolver()
        strain_boundary = [[-user_request["strain"], 0.1], [user_request["strain"], 0.1]]
        solver1.model_par["strain_boundary"] = strain_boundary
        solver1.nr_epochs = user_request["nr_epochs"]

        # Calculate super cell
        super_cell1 = solver1.solve(user_request["cel1"], user_request["cel2"])

        print(super_cell1.description_txt())
        response = super_cell1.description()
        response["algo"]="annealing"
        response["max_strain"]= user_request["strain"]
        return response

    if user_request["algo"] == "genetic":
        print("\n\nGeneticSolver:")
        solver1 = GeneticSolver()
        strain_boundary = [[-user_request["strain"], 0.1], [user_request["strain"], 0.1]]
        solver1.model_par["strain_boundary"] = strain_boundary
        solver1.nr_epochs = user_request["nr_epochs"]

        # Calculate super cell
        super_cell1 = solver1.solve(user_request["cel1"], user_request["cel2"])

        print(super_cell1.description_txt())
        response = super_cell1.description_txt()
        response["algo"] = "genetic"
        response["max_strain"] = user_request["strain"]
        return response

    if user_request["algo"] == "mecanic":
        solver1 = MechanicSolver()
        solver1.exploring_range = user_request["nr_epochs"]
        solver1.tolerance = user_request["strain"]

        # Calculate super cell
        super_cell1 = solver1.solve(user_request["cel1"], user_request["cel2"])
        print(super_cell1.description_txt())

        response = super_cell1.description()
        response["algo"] = "mecanic"
        response["max_strain"] = user_request["strain"]
        return response
