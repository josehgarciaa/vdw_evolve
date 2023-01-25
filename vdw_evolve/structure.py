import numpy as np
import json

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