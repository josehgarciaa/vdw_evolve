from django import forms


class ExperimentForm(forms.Form):
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

    # Algorithm
    algo = forms.ChoiceField(choices=(("annealing", "annealing"), ("genetic", "genetic"), ("mecanic", "mecanic")))

    # Strain
    strain = forms.FloatField(min_value=0, max_value=1, label="strain (*)")
    # Nr. epochs
    nr_epochs = forms.IntegerField(min_value=1, max_value=100, label="nr. epochs (*)")
