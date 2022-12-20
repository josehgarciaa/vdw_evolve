from django.db import models


# Create your models here.
# class Experiments(models.Model):
#     cel1 = models.CharField(max_length=100)
#     cel2 = models.CharField(max_length=100)
#     max_strain = models.FloatField()
#     nr_epochs = models.IntegerField()
#
#     # ip_str = models.CharField(max_length=100)
#
#     def __str__(self):
#         return self.cel1


class Experiment(models.Model):
    cel1 = models.CharField(max_length=500)
    cel2 = models.CharField(max_length=500)
    super_cell = models.CharField(max_length=500)
    tA = models.CharField(max_length=500)
    det_tA = models.FloatField()
    tB = models.CharField(max_length=500)
    strain = models.CharField(max_length=500)
    diagonal_strain = models.CharField(max_length=500)
    tB_op = models.CharField(max_length=500)

    algo = models.CharField(max_length=100)
    max_strain = models.FloatField()
    nr_epochs = models.IntegerField()

    def __str__(self):
        return self.super_cell
