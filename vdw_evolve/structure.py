import numpy as np
from numpy.linalg import norm
import json
from .geometry import ChangeBasis
from .lat_matcher import LatMatch
from .parser import cell_fromXYZ, atoms_fromXYZ
from .parser import cell_fromJSON, atoms_fromJSON
from .physics import num2chem

class Structure():
    """
    Structural information of a material
    Methods
    -------
    """
    def __init__(self, cell=None, atoms=None, periodic=(True, True)):
        """
        Attributes
        ----------
        periodic : bool or tuple of bool, optional
            Whereas it is periodic or not along its lattice vectors.
        cell: an 3x3 array of float, optional
        Defining the lattice vectors as columns
        atoms: a list of tuples, optional
        Containing as first element the
        atom type and their posittion in 3D space.
        """
        self.cell = None
        self.atoms = None

    def supercell_points(self, dims):
        """The set of allowed cell points in a supercell defined by dims
            Attributes
            ----------
            dims: 2-tuple of ints
                The number of repeated unit cells along the cell vectors
                that defines the supercell.
        """       
        nx, ny = dims
        grid = np.mgrid[-nx:nx+1, -ny:ny+1, 0:1].T.reshape((nx*2+1)*(ny*2+1), 3).T
        return (self.cell) @ grid
    
    def transform2D(self, strain, angle, transform_atoms=True,
                    angle_format="deg", strain_format="perc"):
        """
            Attributes
            ----------
            strain: float or 3-tuple of float, optional
                The strain in percentile applied to the system.
                When scalar uniform strain will be applied to all axis
                When tuple, each value will be apply to an axis
            angle: float, optional
                The rotation angle along the z direction
            transform_atoms: boolean, optinal
                Apply the transformation to the atoms.
        """
        try:
            one = np.eye(2).astype(float)
            if isinstance(strain, (int, float)):
                StrMat = one*strain
            else:
                StrMat = np.diag(strain).astype(float)
                
            if strain_format == "perc":
                StrMat = StrMat/100.0
            StrMat += one
        except ValueError:
            print("Could not properly parse the strain parameters to a  matrix")
            StrMat = None

        if angle_format == "deg":
            angle = angle*np.pi/180.0
        RotMat = [[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle), np.cos(angle)]]
        try:
            RotMat = np.array(RotMat, dtype=float)
        except ValueError:
            print("Could not properly parse the angle into a Rot matrix")
            RotMat = None

        TrMat = np.eye(3)
        TrMat[:2, :2] = StrMat@RotMat
        self.cell = TrMat@(self.cell)
        if transform_atoms:
            self.atoms = [(s, TrMat@p) for s, p in self.atoms]
        return self

    def rotate2D(self, angle, format="deg"):
        """
            Attributes
            ----------
            angle: float
                The rotation angle along the z direction
            format: str
                The angle format, can be either deg (degrees) or rad (radians)
        """
        return self.transform2D(angle=angle, angle_format=format,
                                strain=0.0)

    def strain2D(self, strain, format="perc"):
        """
            Attributes
            ----------
            strain: float or 3-tuple of float, optional
                The strain in percentile applied to the system.
                When scalar uniform strain will be applied to all axis
                When tuple, each value will be apply to an axis
        """
        return self.transform2D(strain=strain, strain_format=format,
                                angle=0.0)

    def ChangeUnitCell(self, new_cell):
        """Change the unit cell and modify the atoms
        
            Attributes
            ----------
            cell : _type_
                _description_
        """
        sc_points = self.supercell_points(dims=(1000, 1000))
        XX, YY, ZZ = ChangeBasis(sc_points, new_cell)
        inSC = (XX >= 0.0)*(XX < 0.999)*(YY >= 0.0)*(YY < 0.999)
        lat_vecs = np.transpose([XX[inSC], YY[inSC], ZZ[inSC]])

        new_atoms = []
        for lat_vec in lat_vecs:
            for s, p in self.atoms:
                new_atoms += [(s, p+lat_vec)]
        self.atoms = new_atoms
        self.cell = new_cell
        return self

    def read_from(self, filepath, format):
        """
            Attributes
            ----------
            filepath: string  
                The location of the input file to be read.
            format: string
                The format of the input file
        """
        if "-xyz" in format:
            self.cell = cell_fromXYZ(filepath, format)
            self.atoms = atoms_fromXYZ(filepath, format)
        if "-json" in format:
            self.cell = cell_fromJSON(filepath, format)
            self.atoms = atoms_fromJSON(filepath, format)

        return self

    def write_to(self, output_filename, format):
        """
            Attributes
            ----------
            output_file: string
                The location of the output file to be written.
            format: string
                The format of the output file
        """
        if format == "c2db-xyz":
            try:
                npoints = len(self.atoms)
                with open(output_filename, "w") as f:
                    xyz_out = str(npoints)+"\n"
                    xyz_out += "Lattice=\""
                    for lat in (self.cell).T:
                        xyz_out += "{} {} {} ".format(*lat)
                    xyz_out = xyz_out[:-1] + "\""  # The-1 remove the last space and replace
                    xyz_out += "Properties=Generated from vdw_evolve. No other property here\n"  
                    for sym, pos in self.atoms:
                        xyz_out += "{} {} {} {}\n".format(sym, *pos)
                    f.write(xyz_out)
            except ValueError:
                print("Could not properly parse data to the output format")

        if format == "c2db-json":
            json_cell =self.to_dict()
            with open(output_filename, "w") as outfile:
                json.dump(json_cell, outfile, cls=NumpyEncoder)



        return self

    def to_dict(self):
        """
            Attributes
            ----------

            format: dict
                a dictionary with the proprieties of the class
        """
        json_data = {}
        arr_type = "__ndarray__"
        obj_type = "array"

        json_data['1'] = {}
        json_data['1']["cell"] = {}
        json_data['1']["cell"][obj_type] = {arr_type: []}

        lattice = np.transpose(self.cell)
        shape = lattice.shape
        d_type = "float"

        json_data['1']["cell"][obj_type][arr_type].append(shape)  # shape = data_struct[0]
        json_data['1']["cell"][obj_type][arr_type].append(d_type)  # d_type = data_struct[1]
        json_data['1']["cell"][obj_type][arr_type].append(lattice.reshape(-1))
        # lattice = np.array(data_struct[2], dtype=float).reshape(shape)
        # np.transpose(lattice.reshape(shape))


        json_data['1']["positions"] = {arr_type: []}

        shape = (len(self.atoms), 3)
        json_data['1']["positions"][arr_type].append(shape)
        json_data['1']["positions"][arr_type].append(d_type)
        positions = []
        for atom in self.atoms:
            positions.append(atom[1][0])
            positions.append(atom[1][1])
            positions.append(atom[1][2])
        json_data['1']["positions"][arr_type].append(positions)
        # positions = np.reshape(np.array(positions[2]), positions[0])

        json_data['1']["numbers"] = {arr_type: []}
        shape = (len(self.atoms),)
        json_data['1']["numbers"][arr_type].append(shape)
        d_type = "int"
        json_data['1']["numbers"][arr_type].append(d_type)

        chem2num = {}
        for key in num2chem:
            chem2num[num2chem[key]] = key

        atom_n = []
        for atom in self.atoms:
            atom_n.append(chem2num[atom[0]])
        # type = np.reshape(np.array(type[2]), type[0])
        # atom = [num2chem[atom_n], [positions[i][0], positions[i][1], positions[i][2]]]
        json_data['1']["numbers"][arr_type].append(atom_n)
        return json_data
    
    def get_cell(self):
        """
            Return
            -------
            The lattice vectors as the columns of a square matrix.
        """
        return self.cell


class VdWStructure(Structure):

    _strain = (0, 0)
    _angle = 0

    def __init__(self, host, complement, supercell=False):
        """
        Attributes
        ----------
        host : Structure
            The substrate or inmutable lattice that wont be modified
            throught stacking.
        complement : Structure
            The lattice which will be placed on the substrate and adapt to it,
            through by strain and twisting.
        supercell: bool
            If true, the optimal periodic supercell between host and complement will be computed. 
        """
        self.host = host
        self.complement = complement

    def get_cell(self, tol=1e-2):
        """ The lattice vectors as the columns of a square matrix.
        
        Note
        -----
        When the van der Waal supercell is not yet constructed returns none 
        """
        return None

    def supercell_points(self, dims=(300, 300), tol=1e-4):
        """ Lattice points compatible with the host and the compoment
        
        Attributes
        ----------
        dims: 2-tuple of ints
            The number of repeated unit cells along the cell vectors
            that defines the supercell.
        tol : float
            The allowed tolerance to consider a point beloning to both lattices
        Note
        -------
        When no points exist, returns none
        
        TODO
        -----
        The tolerance should be better addressed
        
        """
        
        host, comp = self.host, self.complement
        comp_points = self.complement.supercell_points(dims)
        CinH = ChangeBasis(comp_points, host.cell)
        return comp_points[:, norm(np.round(CinH) - CinH, axis=0) < tol]

    def get_minimalcell(self, dims, optimizer):

        (strain, angle), cell, opt_str = optimizer.minimalCell(self.complement,
                                                               self.host, max_dims=dims)
        opt_vdw = VdWStructure(self.host, opt_str)
        opt_vdw._strain = strain
        opt_vdw._angle = angle
        opt_vdw.host.ChangeUnitCell(cell)
        opt_vdw.complement.ChangeUnitCell(cell)
        opt_vdw.cell = cell
        opt_vdw.atoms = [*opt_vdw.host.atoms, *opt_vdw.complement.atoms]
        return opt_vdw

    def complement_strain(self, format="perc"):
        """Return the strain applied to complement to make the vdw cell 

        Attributes
        ----------
        format : str, optional
            The format which can be perc for percentile and abs for absolute, by default "perc"
        """
        s1, s2 = self._strain
        scal = 1.0
        if format == "perc":
            scal = 100
        return (scal*s1, scal*s2)
    
    def complement_angle(self, format="deg"):
        """Return the angle applied to complement to make the vdw cell 

        Attributes
        ----------
        format : str, optional
            The format which can be deg for degrees and rad for radians, by default "deg"
        """
        if format == "deg":
            return self._angle*180/np.pi
        return self._angle

    def write_to_json(self, output_filename):

        host_dict = self.host.to_dict()
        complement_dict = self.complement.to_dict()
        super_cell_dict = self.to_dict()
        super_cell_dict['host']=host_dict
        super_cell_dict['complement']=complement_dict
        super_cell_dict['complement_strain']=self._strain
        super_cell_dict['complement_angle']=self._angle
        with open(output_filename, "w") as outfile:
            json.dump(super_cell_dict, outfile, cls=NumpyEncoder)










class SuperCell():

    def __init__(self, parents, transformation, strains):
        """
        :param parents: (np.array[],np.array[]) (cel1,cel2)
        :param transformation:(np.array[],np.array[]) (ta,tb)
        :param strains:(np.array[],np.array[]) (strain,diag_strain)
        TA*cel1-TB*cel2 = 0
        super_cell = TA*cel1
        TB=strain*TB_unstrain
        """
        self.parents = parents
        self.ta = transformation[0]
        self.tb = transformation[1]
        self.strain = strains[0]  # exact strain
        self.diagonal_strain = strains[1]  # not so accurate
        self.cell = np.dot(self.ta, parents[0])

    def description(self):
        desc = {
            "cell": self.cell,
            "cell1": self.parents[0],
            "cell2": self.parents[1],
            "ta": self.ta,
            "strain_tb": self.tb,
            "det_ta": np.linalg.det(self.ta),
            "strain": self.strain,
            "tb": np.dot(np.linalg.inv(self.strain), self.tb),
            "diagonal_strain": self.diagonal_strain,
            "diagonal_strain_tb:": np.dot(np.linalg.inv(self.diagonal_strain), self.tb)
        }

        return desc

    def description_txt(self, path=None):
        desc = self.description()
        tx = "\n=== \n "
        for k in desc:
            tx += k + " :\n" + "{}".format(desc[k]) + "\n\n"
        if path!=None:
            #json_object = json.dumps(desc, indent=4)
            with open(path, "w") as outfile:
                json.dump(desc, outfile, cls=NumpyEncoder)
        return tx


# maybe you can suggest a better place for this
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)
