# Generated by Django 5.0.2 on 2024-06-07 10:04

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('scraper_app', '0014_alter_substances_last_modified'),
    ]

    operations = [
        migrations.DeleteModel(
            name='Store_exchange_json',
        ),
    ]