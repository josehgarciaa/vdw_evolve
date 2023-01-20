"""
This is a brut search of the solution, but it is optimized to have lower complexity
than the previous algorithms. Moreover, the complexity of this approach
can be further reduced if, instead of the list of points,
we will keep the optimum points at each step.
"""

import numpy as np
from .utils import brut_strain, strained_proces


def mechanic_super_cell(cel1, cel2, exploring_range=10, tolerance=0.1, paralel_limit=0.00001):
    solution = get_match(cel1, cel2, exploring_range, tolerance, paralel_limit)
    supper_cell = np.array([solution[0]["point2"], solution[1]["point2"]]).transpose()
    t_cel1 = np.dot(supper_cell, np.linalg.inv(cel1))
    strain_t_cel2 = np.dot(supper_cell, np.linalg.inv(cel2))
    strain = brut_strain(strain_t_cel2)

    diagonal_strain = strained_proces(strain_t_cel2, strain_boundary=[[-0.3, 0.3], [-0.3, 0.3]])

    return t_cel1, strain_t_cel2, diagonal_strain, strain


def get_match(cel1, cel2, exploring_range=10, tolerance=0.3, paralel_limit=0.00001):
    c = cel1
    cel1 = cel2
    cel2 = c

    c2 = (cel2[:, 0], cel2[:, 1])

    candidates = []  # TODO: it can be done without list but for the beggining I preffer to keep it
    # in order to get a better view of the data

    for i in range(-exploring_range, exploring_range):
        for j in range(-exploring_range, exploring_range):
            point2 = i * c2[0] + j * c2[1]

            in_c1 = point_in_cell(cel1, point2)

            if in_c1 <= tolerance:
                point_distance = point2[0] ** 2 + point2[1] ** 2
                if point_distance != 0:
                    candidates.append({"point2": point2, "in_c1": in_c1, "point_distance": point_distance})

    candidates = sorted(candidates, key=sort_key_picker(['point_distance', 'in_c1']), )

    #     for c in candidates:
    #         print(c)

    solution = [candidates[0], ]
    k = 1
    while len(solution) < 2 and len(candidates) >= k + 2:
        k += 1
        close = candidates[k]
        min_v = close["point2"]
        if parallel_check(solution[0]["point2"], min_v, trash_hold=paralel_limit) == False:
            solution.append(close)

    return solution


# Helper functions:

def point_in_cell(cell, point):
    """
    Check hpw close is a point to a point of lattice constructed by cell.
    :param cell:
    :param point:
    :return:
    """
    s = decompose_to_basis(point, a=cell[:, 0], b=cell[:, 1])
    s1 = s[0]
    s2 = s[1]

    c_p = round(s1) * cell[:, 0] + round(s2) * cell[:, 1]
    distance_to_close = np.sqrt((c_p[0] - point[0]) ** 2 + (c_p[1] - point[1]) ** 2)

    #     print("s1:{}={}".format(round(s1),s1))
    #     print("s2:{}={}".format(round(s2),s2))
    #     print("distance_to_close:", distance_to_close)
    return distance_to_close


def parallel_check(v1, v2, trash_hold=0.00001):
    t_prod = np.cross(v1, v2)
    ar_2 = np.linalg.norm(t_prod) ** 2

    if ar_2 < trash_hold:
        return True
    else:
        return False


# v1_1= ca1+db1 => c= (v1_1-db1)/a1
# v1_2= ca2+db2 => v1_2= v1_1*a2/a1-db1*a2/a1+db2
# v1_2-v1_1*a2/a1=d(-b1*a2/a1+b2)
# b = (v1_2-v1_1*a2/a1)/(-b1*a2/a1+b2)
def decompose_to_basis(v, a, b):
    """

    :param v:  vector
    :param a:  base vector 1
    :param b:  base vector 2
    :return:
    """
    d = (v[1] - v[0] * a[1] / a[0]) / (b[1] - b[0] * a[1] / a[0])
    c = (v[0] - d * b[0]) / a[0]

    #     print(v-c*a-d*b)
    return [c, d]


def sort_key_picker(keynames):
    negate = set()
    for i, k in enumerate(keynames):
        if k[:1] == '-':
            keynames[i] = k[1:]
            negate.add(k[1:])

    def getit(adict):
        composite = [adict[k] for k in keynames]
        for i, (k, v) in enumerate(zip(keynames, composite)):
            if k in negate:
                composite[i] = -v
        return composite

    return getit
