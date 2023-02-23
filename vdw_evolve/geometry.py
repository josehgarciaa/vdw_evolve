import numpy as np
from scipy import optimize


def ChangeBasis(r, B, A=np.eye(3)):
    """ Changes r from their canonical basis A to a target basis B.

        Parameters
        ----------
        r : array_like
            A two-dimensional vector expressed in the A basis.
        B : matrix_like.
            The target basis in which r will be expressed. 
            If B is the identity, r 
            will be expressed in the coordinates defining A
        A : matrix_like, optional
            The basis of r.
        
        Return
        ----------
        rt : array_like
            the r vector expressed in the B basis. 

        Notes
        -----
        A and B are matrices containing the basis vectors expressed 
        in a common basis as columns.
        
        This function takes the vector r=(\alpha_1,alpha_2) 
        in the basis A basis, i.e:
        $\vec{r}_i = \alpha_1 \vec{a}_1 + \alpha_2 \vec{a}_2$
        
        and compute the coefficients $(\beta_1,\beta_2)$, such that
        $\vec{r}_i = \beta_1 \vec{b}_1 + \beta_2 \vec{b}_2$.
        
        This is achieved by taking the dot product
        of $\vec{r}_i$ agains $\vec{b}_j$
        $\vec{b}_i\cdot\vec{b}_j \beta_i = \vec{b}_i\cdot\vec{a}_j \alpha_j$.
        which defines change of basis matrix $U_{ij} = \vec{b}_i\cdot\vec{a}_j$
    """
    r = np.array(r)
    U = np.linalg.inv(B.T @ B).dot(B.T @ A)
    return U @ r


# Strain and Rotation
def R(theta):
    return np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])


def S(s1, s2):
    return np.diag([1 + s1, 1 + s2])


def SR(s1, s2, theta):
    return S(s1, s2) @ R(theta)


def get_supercell_vectors(dims, ref, target, tol=1e-2):
    """
    Received a set of points x in the reference lattice basis and returns those
    that are close to the reference lattice points given a tolerance tol.
    """
    B, A = ref, target;
    Bpoints = supercell_points(dims, B)  # This points are in cartesians
    rBinA = ChangeBasis(Bpoints, A)  # Map from cartesian to Alat coordinates
    return Bpoints[:, np.linalg.norm(np.round(rBinA) - rBinA, axis=0) < tol]


def supercell_points(dims, lat_vec, fractional=False):
    nx, ny = dims
    grid = np.mgrid[-nx:nx, -ny:ny].T.reshape(nx * 2 * ny * 2, 2).T
    if fractional:
        return grid
    return lat_vec @ grid