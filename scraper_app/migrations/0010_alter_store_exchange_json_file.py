# Generated by Django 5.0.2 on 2024-04-26 19:00

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('scraper_app', '0009_alter_store_exchange_json_file'),
    ]

    operations = [
        migrations.AlterField(
            model_name='store_exchange_json',
            name='file',
            field=models.FileField(upload_to='./files_from_others/'),
        ),
    ]
