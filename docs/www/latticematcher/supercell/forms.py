from django import forms
from django.db import models

class SimpleForm(forms.Form):
    # Cel1
    c1_00 = forms.FloatField()
    c1_01 = forms.FloatField()
    c1_10 = forms.FloatField()
    c1_11 = forms.FloatField()

    # Cel2
    c2_00 = forms.FloatField()
    c2_01 = forms.FloatField()
    c2_10 = forms.FloatField()
    c2_11 = forms.FloatField()

    # Angle
    angle = forms.FloatField(min_value=0, max_value=180, label="angle (*)")

    # Strain
    strain = forms.FloatField(min_value=0, max_value=1, label="strain (*)",step_size=0.01)

class FileForm(forms.Form):
    # Cel1
    file1 = forms.FileField(label="cell1 file (*)", )

    # Cel2
    file2 = forms.FileField(label="cell2 file (*)", )

    # Angle
    angle = forms.FloatField(min_value=0, max_value=180, label="angle (*)")

    # Strain
    strain = forms.FloatField(min_value=0, max_value=1, label="strain (*)",step_size=0.01)

