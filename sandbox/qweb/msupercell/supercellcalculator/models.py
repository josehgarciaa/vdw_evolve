from django.db import models


# Create your models here.
class Experiments(models.Model):
    cel1 = models.CharField(max_length=100)
    cel2 = models.CharField(max_length=100)
    max_strain = models.FloatField()
    nr_epochs = models.IntegerField()
    # ip_str = models.CharField(max_length=100)

    def __str__(self):
        return self.cel1
