.. vdw_evolve documentation master file, created by
   sphinx-quickstart on Sat Oct 15 13:03:55 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Examples 
========


Minimal Supercell
___________________

A Bravais lattice is a discrete collection of infinite points all related by a set of discrete translations follow a pattern. For instance, in two dimensions 
this means that a lattices is the collections of points :math:`R(n_1,n_2)` that satisfies the following equation

.. math::
   R(n_1,n_2) = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2, 

where :math:`{\mathbf{a}_1}` and :math:`{\mathbf{a}_2}` two noncollinear two-dimensional vectors known as the lattice vectors. An example of a lattice is presented
below:

.. image:: images/example0_honeycomb_lattice.svg
   :width: 300

Unit cell
*********

For convenience, let us group the lattice vectors as columns in a :math:`2\times2` matrix 

.. math::
   A = \left( \begin{array}{cc}
               \mathbf{a}_{1,1} & \mathbf{a}_{2,1} \\
               \mathbf{a}_{1,2} & \mathbf{a}_{2,1}
      \end{array} \right).

The determinant of matrix :math:`A` describes the area of the plane formed by the lattice vectors, and within this area there is a collection of points. 
This collection of points is known as unit cell, and represent a possible subset of points that can be translate with these unit vectors that fully represents 
the lattice. We will call the matrix :math:`A` unit cell matrix. 


Van der Waal cells
******************

Let us consider know, another lattice B, represent by the unit cell matrix

.. math::
   B = \left( \begin{array}{cc}
               \mathbf{b}_{1,1} & \mathbf{b}_{2,1} \\
               \mathbf{b}_{1,2} & \mathbf{b}_{2,1}
      \end{array} \right).

for instance, this could be a lattice with a different symmetry such a the square Lattice.

.. image:: images/example0_square_lattice.png
   :width: 300

If we stack these two lattices together, a new lattice may form, but its not guarantee since the lattice can be incomesurate, meaning that
the ratio of their periods along any given direction in space is not an integer and therefore the entire collection of points won't repite periodically. 

Let us work under the assumption the are commesurate. Therefore, there exist two transformation :math:`T_A` and :math:`T_B` that will generate the same set of points
of the super cell :math:`R_{sc}`, i.e

.. math::
   R_{sc} = T_A A = T_B B,

where 

.. math::
   T_\alpha = \left( \begin{array}{cc}
   i_\alpha & j_\alpha \\
   k_\alpha & l_\alpha
   \end{array} \right),\quad \alpha=A,B
   :label: super_cell_eq

just discrete linear transformations. Its clear that if :math:`A` describe the points of a bravai lattice, a 
linear combination of the lattice vectors :math:`T_A A` will describe another allowed unit cell. 
The role of equation :eq:`super_cell_eq` is two transform the two initial unit cell represents into a new common cell.  

The problem
************

For two arbitrary sets lattice vector sets :math:`A` and :math:`B`, find the set of 8 integers :math:`I:=(i_A,i_A,j_A,J_B,k_B,k_B,l_A,l_B)`  
that satisfy

.. math::
   T_A A - T_B B

with minium  :math:`|det(T_A A)|>0`.   


Optimal Fermi Contour
_____________________

Another interesting approach for the evolutionary strategy is to efficiently sample a Fermi Contour. The Brillouin zone is the momentum-space 
equivalent of a Bravais lattice, where the allowed momentums satisfies:

.. math::
   k(n_1,n_2) = \frac{n_1}{N_1} \mathbf{b}_1 + \frac{n_2}{N_2} \mathbf{b}_2, 

where :math:`N_1` and :math:`N_2` the number of unit cells contained in your crystal , and the reciprocal lattice vector
are defined as the solution of the following equation

.. math:: 
   B = \left( \begin{array}{cc}
               \mathbf{b}_{1,1} & \mathbf{b}_{2,1} \\
               \mathbf{b}_{1,2} & \mathbf{b}_{2,1}
      \end{array} \right) A = 2\pi \left( \begin{array}{cc}
               1 & 0 \\
               0 & 1
      \end{array} \right).

Let us now consider that the Hamiltonian matrix in the momentum space is given by the following expression

.. math:: 
   H(\textbf{k}) = \left( \begin{array}{cc}
               0             & f(\textbf{k}) \\
               f^*(\textbf{k}) & 0
      \end{array} \right),

where :math:`f(k) = {\rm e}^{k\textbf{k}\cdot \textbf{a}_1}+ {\rm e}^{k\textbf{k}\cdot \textbf{a}_2}+ 1` and :math:`f^*(\textbf{k})` represents its 
complex conjugate

An important quantity for many application is to find the set of :math:`k_F{\rm s}`  that satisfy the following constrain

.. math:: 
   {\rm eigval}[ H(\textbf{k}_F) -  E_F \mathbb{I} ]= 0

where :math:`E_F` the Fermi level of the system, where :math:`\mathbb{I}` the identity matrix, 
eigval a function that compute the eigenvalues of the system and the set :math:`\{k_F\}` the so-called Fermi momentum. 

The problem
***********

For a given system size :math:`(N_1,N_2)` (typicall above 1000 ), look for the set of :math:`(n_1,n_2)` with :math:`n_1=0,\dots, N_1` and :math:`n_2=0,\dots, N_2` that satisfies:

.. math::
   {\rm eigval}\left[ H((\textbf{k})(n_1,n_2)) - E_F \mathbb{I} \right] =0


Please use :math:`a_1=(1/2,\sqrt{3}/2)` and  :math:`a_2=(1/2,-\sqrt{3}/2)`


   
