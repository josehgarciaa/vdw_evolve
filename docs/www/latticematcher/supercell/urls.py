from django.urls import path

from . import views

urlpatterns = [
    path("basic_matcher", views.basic_calcualtor, name="basic"),
    path("file_matcher", views.file_calcualtor, name="file"),
    path("calculate_basic", views.calculate_basic, name="basic_c"),
    path("calculate_file", views.calculate_file, name="file_c"),
]