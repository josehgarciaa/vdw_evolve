# Generated by Django 4.1.4 on 2022-12-10 18:23

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('supercellcalculator', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='Experiment',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cel1', models.CharField(max_length=500)),
                ('cel2', models.CharField(max_length=500)),
                ('super_cell', models.CharField(max_length=500)),
                ('ta', models.CharField(max_length=500)),
                ('det_ta', models.FloatField()),
                ('tb', models.CharField(max_length=500)),
                ('strain', models.CharField(max_length=500)),
                ('diagonal_strain', models.CharField(max_length=500)),
                ('tb_op', models.CharField(max_length=500)),
                ('max_strain', models.FloatField()),
                ('nr_epochs', models.IntegerField()),
                ('algo', models.CharField(max_length=100)),
            ],
        ),
    ]
