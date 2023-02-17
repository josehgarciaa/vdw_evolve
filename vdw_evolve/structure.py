import numpy as np
import json


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
    
    def transform2D(self, strain=1, angle=0):
        """
            Attributes
            ----------
            strain: float or 3-tuple of float, optional
                The strain in percentile applied to the system.
                When scalar uniform strain will be applied to all axis
                When tuple, each value will be apply to an axis
            angle: float, optional
                The rotation angle along the z direction
        """
        try:
            if isinstance(strain, (int, float)):
                StrMat = np.eye(2, cdtype=float)*strain
            StrMat = np.diag(1.0 + np.array(strain)/100.0)
        except ValueError:
            print("Could not properly parse the strain parameters to a  matrix")
            StrMat = None

        RotMat = [[np.cos(angle), -np.sin(angle)],
                  [np.sin(angle), np.cos(angle)]]
        try:
            RotMat = np.array(RotMat, dtype=float)
        except ValueError:
            print("Could not properly parse the angle into a Rot matrix")
            RotMat = None
        TrMat = StrMat@RotMat
        self.cell = list(TrMat@self.cell)
        self.atoms = [(s, list(TrMat@p)) for s, p in self.atoms];
        return self

    def rotate2D(self, angle):
        """
            Attributes
            ----------
            angle: float
                The rotation angle along the z direction
        """
        return self.transform(angle)

    def strain2D(self, strain):
        """
            Attributes
            ----------
            strain: float or 3-tuple of float, optional
                The strain in percentile applied to the system.
                When scalar uniform strain will be applied to all axis
                When tuple, each value will be apply to an axis
        """
        return self.transform(strain)
    
    
    def read_from(self, input_file, format):
        """
            Attributes
            ----------
            input_file: string  
                The location of the input file to be read.
            format: string
                The format of the input file
        """
        if format == "c2db-xyz":
            try:
                # A temporal function to format a line
                def _format(line):
                    s, x, y, z = [x for x in line.split(" ") if x != ""][:4]
                    return (s, [float(x), float(y), float(z)])

                with open(input_file) as f:
                    npoints = int(f.readline())
                    lat = f.readline()
                    lat = lat[lat.find("\"")+1:lat.find("Properties")-2]
                    lat = lat.split(" ")
                    lat = np.array(list(map(float, lat))).reshape(3, 3)
                    xyz = [_format(line) for line in f]
                    self.cell = np.transpose(lat)
                    self.atoms = xyz
            except ValueError:
                print("Could not properly parse input_file")
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
        return self


        


    





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