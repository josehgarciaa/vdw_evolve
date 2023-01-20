from django.db import models





class Experimentt(models.Model):

    cell = models.CharField(max_length=500)
    cell1 = models.CharField(max_length=500)
    cell2 = models.CharField(max_length=500)
    ta = models.CharField(max_length=500)
    strain_tb = models.CharField(max_length=500)
    det_ta = models.FloatField()
    strain = models.CharField(max_length=500)
    tb = models.CharField(max_length=500)
    diagonal_strain = models.CharField(max_length=500)
    diagonal_strain_tb = models.CharField(max_length=500)

    algo = models.CharField(max_length=100)
    max_strain = models.FloatField()
    nr_epochs = models.IntegerField()

    def __str__(self):
        return self.super_cell
