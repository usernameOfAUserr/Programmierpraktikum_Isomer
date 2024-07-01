# Generated by Django 5.0.2 on 2024-06-30 20:37

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('scraper_app', '0018_substances_id_alter_substances_source'),
    ]

    operations = [
        migrations.AlterField(
            model_name='substances',
            name='cas_num',
            field=models.CharField(blank=True, default='0', max_length=1000, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='categories',
            field=models.JSONField(blank=True, default=list, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='formula',
            field=models.CharField(blank=True, max_length=1000, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='inchi',
            field=models.CharField(blank=True, max_length=1000, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='inchi_key',
            field=models.CharField(blank=True, max_length=1000, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='iupac_names',
            field=models.JSONField(blank=True, default=list, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='molecular_mass',
            field=models.FloatField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='names',
            field=models.JSONField(blank=True, default=list, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='smiles',
            field=models.CharField(blank=True, max_length=1000, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='source',
            field=models.JSONField(blank=True, null=True),
        ),
        migrations.AlterField(
            model_name='substances',
            name='validated',
            field=models.BooleanField(blank=True, default=False, null=True),
        ),
    ]