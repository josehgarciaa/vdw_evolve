"""
Each solver it's a routine that will take two cells (2x2 np. arrays) and will return a supper cell object
"""

from .sc_annealing import super_cell as annealing_sc
from .sc_mechanic import mechanic_super_cell as mechanic_sc
from .utils import allign_along_10
