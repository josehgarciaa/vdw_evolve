from django.urls import path

from . import views

urlpatterns = [
    path("cell_calculator/", views.cell_calculator, name="cell_calculator"),
    path("calculate/", views.calculate, name="calculate")
]

# http://127.0.0.1:8000/cell_calculator/

# http://127.0.0.1:8000/admin/