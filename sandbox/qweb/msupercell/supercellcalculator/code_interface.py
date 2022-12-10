import numpy as np
from .super_cell import annealing_super_cell


def get_model_param(algo, strain):
    if algo == "annealing":
        up = 99999999  # max_integer/float. Largest Number in the structure
        bond = 30  # Supercell bounds
        model_par = {

            'initialTemp': 4,  # exp(- 1/(T+beta) ) Check later
            'finalTemp': 0.0002,  #

            'beta': 10,  # How fast the temperature decrease
            'bounds': [[-bond, bond] for _ in range(4)],

            'nr_neighbours': 1,  # Number of randomly choosen neighbors ( Depreceated )
            'step_size': 3,
            # How much do you go in a particular direction. This stepsize will choose the peak of the distribution
            'gaussian_sd': 3,  # Standandard deviation

            'known_min': -up,  # the "minumum" number in the program
            "start_point": [1, 0, 0, 1],  # Initial parameters. IMPORTANT Initial T_A = [ [1 0], [0,1] ];

            "strain_boundary": [[-strain, strain], [-strain, strain]]
            # [[-0.3,0.3],[-0.3,0.3]]#[[-0.3e-14,0.3e-14],[-0.3e-14,0.3e-14]]
        }
        return model_par

    if algo == "genetic":
        return 0

    if algo == "deepnet":
        print("Not implemented yet!")
        return 0


def run_experiment(user_request):
    if user_request["algo"] == "annealing":
        model_parm = get_model_param(user_request["algo"], user_request["strain"])
        tA, tB, t_cel2_no_strain, diagonal_strain, strain = annealing_super_cell(user_request["cel1"],
                                                                                 user_request["cel2"],
                                                                                 user_request["nr_epochs"],
                                                                                 model_parm)
        response = {"cel1": user_request["cel1"],
                    "cel2": user_request["cel2"],
                    "super_cell": np.dot(tA, user_request["cel1"]).tolist(),
                    "tA": tA.tolist(),
                    "det_tA": np.linalg.det(tA),
                    "tB": tB.tolist(),
                    "strain": strain.tolist(),
                    "diagonal_strain": diagonal_strain.tolist(),
                    "tB_op": t_cel2_no_strain.tolist()}

        return response

    if user_request["algo"] == "genetic":
        pass

    if user_request["algo"] == "deepnet":
        # TODO: connect to the  deep net routine
        pass
