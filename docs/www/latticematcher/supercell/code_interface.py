import numpy as np
import vdw_evolve as lt
from .models import Experiment
from django.core.files.storage import default_storage

def demo_solver(cel1, cel2, max_angle, max_strain):
    str1 = lt.Structure(cell=np.array(cel1))
    str2 = lt.Structure(cell=np.array(cel2))
    print(str1)
    vdws = lt.VdWStructure(str1, str2)
    print(str2)
    if max_angle == 0:
        optimizer = lt.LatMatch().opt_angle(False).uniform_strain(max_strain, "perc")
    if max_strain == 0:
        optimizer = lt.LatMatch().opt_strain(False)
    elif max_angle != 0:
        optimizer = lt.LatMatch().uniform_strain(max_strain, "perc")

    # Get the minimal cell
    dims = (max_strain, max_strain)
    optVdW = vdws.get_minimalcell(dims=dims, optimizer=optimizer)
    strain = optVdW.complement_strain(),
    angle = optVdW.complement_angle()

    return angle, strain


def file_solver(file1, file2, max_angle, max_strain):

    if file_json(file1):
        str1 = lt.Structure().read_from(file1, format="c2db-json")
    else:
        str1 = lt.Structure().read_from(file1, format="c2db-xyz")

    if file_json(file2):
        str2 = lt.Structure().read_from(file2, format="c2db-json")
    else:
        str2 = lt.Structure().read_from(file2, format="c2db-xyz")

    vdws = lt.VdWStructure(str1, str2)

    if max_angle == 0:
        optimizer = lt.LatMatch().opt_angle(False).uniform_strain(max_strain, "perc")
    if max_strain == 0:
        optimizer = lt.LatMatch().opt_strain(False)
    elif max_angle != 0:
        optimizer = lt.LatMatch().set_angle_range((-max_angle, max_angle), "deg").uniform_strain(max_strain, "perc")

    # Get the minimal cell
    dims = (50, 50)
    optVdW = vdws.get_minimalcell(dims=dims, optimizer=optimizer)
    strain = optVdW.complement_strain(),
    angle = optVdW.complement_angle()

    f1 =get_file_name(file1)
    f2 =get_file_name(file2)
    name = f1+"-a-"+f2
    # optVdW.write_to_json("/SuperCells/"+f1+"-"+f2+".json")
    return optVdW, name

def file_json(file_path):
    file= file_path.split("/")[-1]
    ext = file.split(".")[-1]
    if ext == "json":
        return True
    else:
        return  False
def get_file_name(file_path):
    file = file_path.split("/")[-1]
    file_name = file.split(".")[-2]
    return file_name
def update_db(ip, cel1=None, cel2=None,angle=None,strain=None,super_cell=None):

    if super_cell is not None:
        pass

    else:
        pass
    cel1_=""
    for r in cel1:
        for l in r:
            cel1_+=str(l)+","
        cel1_ += ";"

    cel2_ = ""
    for r in cel2:
        for l in r:
            cel2_ += str(l) + ","
        cel2_ += ";"

    exp = Experiment(ip_name = ip,
        parent_lattice=cel1_,
        parent_name=0,
        host_name=0,
        angle_range=angle,
        strain_max=strain,

        )
    exp.save()
