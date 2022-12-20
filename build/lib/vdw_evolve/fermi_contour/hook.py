"""
We need to put this in another place.

we have a function and a solution and the task is to fond all connected solutions.
"""


def discrete_hook(obj_function, solution, boundary, step=1, step_range=2,tolerance =0):
    """

    :param obj_function:
    :param solution:
    :param boundary:
    :param step:
    :param step_range:
    :return:
    """

    minimum = obj_function(solution)+tolerance
    # print("minim", minimum)

    solution_points = [solution]
    solution_values = [minimum]
    to_explore = [solution]
    explored_points = []

    while len(to_explore) != 0:
        exploration_point = to_explore[0]

        neighbours = get_neighbours(exploration_point, step, step_range)
        #filter
        filtered=[]
        for n in neighbours:
            for  p_id in range(len(n)):
                if  boundary[p_id][0]<=n[p_id]<boundary[p_id][1]:
                    filtered.append(n)
        neighbours=filtered
        # print("filterd_neighbours:",neighbours )

        for neighbour in neighbours:
            if neighbour not in explored_points:
                explored_points.append(neighbour)
                if neighbour not in to_explore:
                    # print("llll",neighbour)
                    s_val =obj_function(neighbour)
                    if s_val <= minimum:
                        solution_points.append(neighbour)
                        solution_values.append(s_val)
                        to_explore.append(neighbour)


                        # print("good_neighbour:",neighbour)
        # print("len to explore",len(to_explore))
        to_explore.pop(0)

    return solution_points, solution_values


def get_neighbours(solution, step, step_range):
    """
    At the moment it only work with two parameters finctions
    rezoin i am lazy.
    :param solution:
    :param step:
    :param step_range:
    :return:
    """

    s1 = solution.copy()
    s1[0] += 1

    s2 = solution.copy()
    s2[1] += 1

    s3 = solution.copy()
    s3[0] -= 1

    s4 = solution.copy()
    s4[1] -= 1

    s5 = solution.copy()
    s5[0] -= 1
    s5[1] -= 1

    s6 = solution.copy()
    s6[0] += 1
    s6[1] += 1

    s7 = solution.copy()
    s7[0] += 1
    s7[1] -= 1

    s8 = solution.copy()
    s8[0] -= 1
    s8[1] += 1

    s9 = solution.copy()
    s9[0] -= 2
    s9[1] -= 2

    s10 = solution.copy()
    s10[0] += 2
    s10[1] += 2

    s11 = solution.copy()
    s11[0] += 2
    s12 = solution.copy()
    s12[1] += 2

    neighbours = [s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11]
    # print("neighbours:", neighbours)
    return neighbours
