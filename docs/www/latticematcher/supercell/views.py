import os
import numpy as np 
from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.core.files.storage import default_storage
from django.conf import settings

from .models import Experiment, Document
from .forms import FileForm, SimpleForm
from .code_interface import demo_solver, file_solver, update_db
from datetime import  datetime
import shutil

def index(response):
    return HttpResponse("random text!")


def basic_calcualtor(response):
    form = SimpleForm()
    return render(response, 'supercell/supercellcalcualtor.html', {"form": form})


def file_calcualtor(response):
    form = FileForm()
    return render(response, 'supercell/supercellcalculator_files.html', {"form": form})


def calculate_basic(response):
    # Form check
    if response.method == "POST":
        form = SimpleForm(response.POST)

        if form.is_valid():
            ip = get_client_ip(response)
            cel1 = [[form.cleaned_data["c1_00"], form.cleaned_data["c1_01"], 0],
                    [form.cleaned_data["c1_10"], form.cleaned_data["c1_11"], 0],
                    [0, 0, 1]]
            cel2 = [[form.cleaned_data["c2_00"], form.cleaned_data["c2_01"], 0],
                    [form.cleaned_data["c2_10"], form.cleaned_data["c2_11"], 0],
                    [0, 0, 1]]
            max_angle = form.cleaned_data["angle"]
            max_strain = float(form.cleaned_data["strain"])
            media_path = store_directory(ip, max_angle, max_strain)

            # angle, strain = demo_solver(cel1, cel2, max_angle, max_strain)
            file1_path= cell_to_xyx(cel1, media_path, "file1.xyz")
            file2_path = cell_to_xyx(cel2, media_path, "file2.xyz")
            super_cell, name = file_solver(file1_path, file2_path, max_angle, max_strain)
            best_angle = super_cell.complement_angle()
            best_strain =super_cell.complement_strain()

            # update the database
            # update_db(ip, cel1, cel2, angle, strain)

            return render(response, 'supercell/supercellcalcualtor.html', {"form": form,"angle":round(best_angle,3),"strain":(round(best_strain[0]),round(best_strain[1]))})
            # return basic_calcualtor(response)
            #return render(response, 'supercell/supercellcalculator.html', {"form": form})
#render(response, 'supercell/supercellcalculator.html', {"form": form,"angle":best_angle,"strain":best_strain})

    return basic_calcualtor(response)
def calculate_file(response):
    if response.method == "POST":
        import pprint
        pprint.pprint(response.POST)
        pprint.pprint(response.FILES)

        form = FileForm(data=response.POST, files=response.FILES)
        print(form.errors)

        if form.is_valid():
            ip = get_client_ip(response)

            file1 = form.cleaned_data["file1"]
            file2 = form.cleaned_data["file2"]
            file_name1 = default_storage.save(file1.name, file1)
            file_name2 = default_storage.save(file2.name, file2)

            max_angle = form.cleaned_data["angle"]
            max_strain = float(form.cleaned_data["strain"])
            media_path = store_directory(ip, max_angle, max_strain)
            shutil.move(str(settings.MEDIA_ROOT)+"/" + file_name1, media_path + "/" + file_name1 )
            shutil.move(str(settings.MEDIA_ROOT) + "/" + file_name2, media_path + "/" + file_name2)

            file1_path = media_path + "/" + file_name1
            file2_path = media_path + "/" + file_name2



            super_cell, name = file_solver(file1_path, file2_path, max_angle, max_strain)
            # json_cell = media_path + "/SuperCells/" + name + ".json"
            # super_cell.write_to_json(json_cell)
            xyz_cell = media_path + "/SuperCells/" + name + ".xyz"
            super_cell.write_to(xyz_cell, format="c2db-xyz")

            d_file = FileResponse(open(xyz_cell, 'rb'),as_attachment=True)



            #update_db(ip,cel1=super_cell.host.cell, cel2=super_cell.complement.cell,angle=super_cell.complement_angle(),strain=super_cell.complement_angle(),super_cell=super_cell.cell)

            return d_file

    return render(response, 'supercell/supercellcalculator_files.html', {"form": form})


def get_client_ip(request):
    x_forwarded_for = request.META.get('HTTP_X_FORWARDED_FOR')
    if x_forwarded_for:
        ip = x_forwarded_for.split(',')[0]
    else:
        ip = request.META.get('REMOTE_ADDR')
    return ip

def store_directory(ip,max_angle, max_strain):
    now = datetime.now()
    date_time_str = now.strftime("%Y-%m-%d-%H-%M-%S")
    path=str(settings.MEDIA_ROOT)+"/{}_{}_{}_{}_{}".format(ip, date_time_str, max_angle, max_strain, np.random.randint(0,9))
    os.mkdir(path)
    os.mkdir(path+"/SuperCells")
    return path

def cell_to_xyx(cell, partial_path,name):
    file_pth= partial_path+"/"+name
    f = open(file_pth, "a")
    f.write("1\n")
    # Lattice = "Lattice=\"3.5126690426238314 2.0280404173329334 0.0  -3.512681993375633 2.028062848693055 0.0  3.9535022123585904e-21 2.4764297047182165e-15  20.22163904\""
    Lattice="Lattice=\"{} {} {} {} {} {} {} {} {}\"".format(float(cell[0][0]),float(cell[0][1]),float(cell[0][2]),
                                                  float(cell[1][0]),float(cell[1][1]),float(cell[1][2]),
                                                  float(cell[2][0]),float(cell[2][1]),float(cell[2][2]),)
   # Lattice = "3.5126690426238314 2.0280404173329334 0.0" \
    #          " -3.512681993375633 2.028062848693055 0.0 " \
    #         "3.9535022123585904e-21 2.4764297047182165e-15 20.22163904"

    f.write(Lattice)
    f.write("Properties = Generated for vdw_evolve.No other property here\n")
    f.write("V 1.1706582931978193 2.7041709592936956 12.686315530583629\n")
    f.close()
    return file_pth

