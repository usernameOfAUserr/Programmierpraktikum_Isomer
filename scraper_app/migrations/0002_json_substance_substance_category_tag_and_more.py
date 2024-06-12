# Generated by Django 5.0.2 on 2024-03-22 16:14

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('scraper_app', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='json_substance',
            fields=[
                ('id', models.BigAutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('all', models.JSONField()),
            ],
        ),
        migrations.AddField(
            model_name='substance',
            name='category_tag',
            field=models.CharField(default='not_known', max_length=50),
        ),
        migrations.AddField(
            model_name='substance',
            name='formular',
            field=models.CharField(default='not_known', max_length=100),
        ),
        migrations.AddField(
            model_name='substance',
            name='inchl',
            field=models.CharField(default='not_known', max_length=100),
        ),
        migrations.AddField(
            model_name='substance',
            name='inchl_key',
            field=models.CharField(default='not_known', max_length=100),
        ),
        migrations.AddField(
            model_name='substance',
            name='iupac_names',
            field=models.CharField(default='not_known', max_length=500),
        ),
        migrations.AddField(
            model_name='substance',
            name='molecular_weight',
            field=models.IntegerField(default=0),
            preserve_default=False,
        ),
        migrations.AddField(
            model_name='substance',
            name='smiles',
            field=models.CharField(default='not_known', max_length=100),
        ),
        migrations.AlterField(
            model_name='substance',
            name='id',
            field=models.IntegerField(primary_key=True, serialize=False),
        ),
        migrations.AlterField(
            model_name='substance',
            name='names',
            field=models.CharField(default='not_known', max_length=500),
        ),
    ]
