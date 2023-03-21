"""

"""
import json

from django.shortcuts import render
from django.http import HttpResponse

from .forms import ExperimentForm

from .code_interface import run_experiment
from .models import Experimentt

# Create your views here.


# render the main page:
def cell_calculator(response):
    form = ExperimentForm()
    return render(response, "supercellcalculator/supercellcalculator_files.html", {"form": form})


# compute the form
def calculate(response):
    if response.method == "POST":
        form = ExperimentForm(response.POST)
        if form.is_valid():
            # run experiment:
            cel1 = [[form.cleaned_data["c1_00"], form.cleaned_data["c1_01"]],
                    [form.cleaned_data["c1_10"], form.cleaned_data["c1_11"]]]
            cel2 = [[form.cleaned_data["c2_00"], form.cleaned_data["c2_01"]],
                    [form.cleaned_data["c2_10"], form.cleaned_data["c2_11"]]]
            algo = form.cleaned_data["algo"]
            strain = form.cleaned_data["strain"]
            nr_epochs = form.cleaned_data["nr_epochs"]
            user_request = {
                "cel1": cel1, "cel2": cel2, "algo": algo, "strain": strain, "nr_epochs": nr_epochs,
            }
            exp_results = run_experiment(user_request)

            # save the experiment in to data base:
            exp = Experimentt(**exp_results)
            exp.save()
            # return the experiment results in a file
            return HttpResponse(json.dumps(exp_results), content_type="application/json")
    form = ExperimentForm()
    return render(response, "supercellcalculator/supercellcalculator_files.html", {"form": form})

# Examples
# def index(response):
#     return HttpResponse("random text!")
# def home(response):
#     return render(response, "supercellcalculator/extaension_example.html", {})
