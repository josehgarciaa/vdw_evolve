from django.db import models

class Experiment(models.Model):

    ip_name = models.CharField(max_length=300)

    parent_lattice = models.CharField(max_length=100)
    parent_name = models.CharField(max_length=50)

    c = models.CharField(max_length=100)
    host_name = models.CharField(max_length=50)

    angle_range = models.CharField(max_length=50)
    strain_max = models.CharField(max_length=50)
    angle = models.FloatField()
    c= models.FloatField()

    def __str__(self):
        return self.ip_name

class Document(models.Model):
    docfile = models.FileField(upload_to='documents/')


