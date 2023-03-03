# Generated by Django 4.1.5 on 2023-01-20 19:07

from django.db import migrations, models


class Migration(migrations.Migration):

    initial = True

    dependencies = [
    ]

    operations = [
        migrations.CreateModel(
            name='Experimentt',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('cell', models.CharField(max_length=500)),
                ('cell1', models.CharField(max_length=500)),
                ('cell2', models.CharField(max_length=500)),
                ('ta', models.CharField(max_length=500)),
                ('strain_tb', models.CharField(max_length=500)),
                ('det_ta', models.FloatField()),
                ('strain', models.CharField(max_length=500)),
                ('tb', models.CharField(max_length=500)),
                ('diagonal_strain', models.CharField(max_length=500)),
                ('diagonal_strain_tb', models.CharField(max_length=500)),
                ('algo', models.CharField(max_length=100)),
                ('max_strain', models.FloatField()),
                ('nr_epochs', models.IntegerField()),
            ],
        ),
    ]
