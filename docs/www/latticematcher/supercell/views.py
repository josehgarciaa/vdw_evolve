import os

from django.shortcuts import render
from django.http import HttpResponse, FileResponse
from django.core.files.storage import default_storage
from django.conf import settings

from .models import Experiment, Document
from .forms import FileForm, SimpleForm
from .code_interface import demo_solver, file_solver, update_db


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
            cel1 = [[form.cleaned_data["c1_00"], form.cleaned_data["c1_01"], 0],
                    [form.cleaned_data["c1_10"], form.cleaned_data["c1_11"], 0],
                    [0, 0, 1]]
            cel2 = [[form.cleaned_data["c2_00"], form.cleaned_data["c2_01"], 0],
                    [form.cleaned_data["c2_10"], form.cleaned_data["c2_11"], 0],
                    [0, 0, 1]]

            max_angle = form.cleaned_data["angle"]

            max_strain = form.cleaned_data["strain"]

            angle, strain = demo_solver(cel1, cel2, max_angle, max_strain)

            ip = "x"
            # update the database
            update_db(ip, cel1, cel2, angle, strain)

            print(angle)


def calculate_file(response):
    if response.method == "POST":
        import pprint
        pprint.pprint(response.POST)
        pprint.pprint(response.FILES)

        form = FileForm(data=response.POST, files=response.FILES)
        print(form.errors)

        if form.is_valid():
            file1 = form.cleaned_data["file1"]
            file2 = form.cleaned_data["file2"]
            file_name1 = default_storage.save(file1.name, file1)
            file_name2 = default_storage.save(file2.name, file2)

            file1_path = str(settings.MEDIA_ROOT) + "/" + file_name1
            file2_path = str(settings.MEDIA_ROOT) + "/" + file_name2

            max_angle = form.cleaned_data["angle"]

            max_strain = float(form.cleaned_data["strain"])

            super_cell, name = file_solver(file1_path, file2_path, max_angle, max_strain)
            json_cell = str(settings.MEDIA_ROOT) + "/SuperCells/" + name + ".json"
            super_cell.write_to_json(json_cell)
            xyz_cell = str(settings.MEDIA_ROOT) + "/SuperCells/" + name + ".xyz"
            super_cell.write_to(xyz_cell, format="c2db-xyz")

            d_file = FileResponse(open(json_cell, 'rb'))

            ip = get_client_ip(response)

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
