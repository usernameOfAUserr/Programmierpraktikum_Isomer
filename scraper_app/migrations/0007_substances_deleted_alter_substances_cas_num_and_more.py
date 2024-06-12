# Generated by Django 5.0.2 on 2024-04-24 11:29

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('scraper_app', '0006_substances_delete_json_substance_delete_substance'),
    ]

    operations = [
        migrations.AddField(
            model_name='substances',
            name='deleted',
            field=models.BooleanField(default=False),
        ),
        migrations.AlterField(
            model_name='substances',
            name='cas_num',
            field=models.IntegerField(blank=True, max_length=1000),
        ),
        migrations.AlterField(
            model_name='substances',
            name='iupac_name',
            field=models.JSONField(blank=True, default=list),
        ),
    ]
