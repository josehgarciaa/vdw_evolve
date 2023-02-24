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


def MinimalBasis(self, vector_list, alpha=0.5, min_area=1e-7):
    """ Extract a minimal basis out of a list of vectors.

        Parameters
        ----------
        vector_list : array_like
            The list of potential vectors to construct the basis.
        alpha : float, default.
            Decided whereas to construct the basis based on its area (alpha=0)
            or on the norm of its lattices vectors (alpha=1). 
            The default is 0.5 which balance both
        min_area : float, optional
            Define the minimal area you'll accept for your supercell.
        
        Return
        ----------
        rt : array_like
            a list of two vectors that forms a basis. Or none if nothing.
    """

    V = vector_list

    basis = []
    cost = np.inf
    for i in range(len(V)):
        v, ws = V[i], V[i+1:]
        areas = np.abs(np.cross([v], ws))
        ws = ws[areas > min_area]
        areas = areas[areas > min_area]
    
        if len(ws) != 0:
            norms = np.linalg.norm(ws, axis=1)
            costs = (1-alpha)*areas + alpha*norms
            min_cost = np.min(costs)
            if min_cost < cost:
                basis = (v, ws[np.argmin(costs)])
                cost = min_cost
    
    if len(basis) == 0:
        return None
    
    return np.array(basis)
