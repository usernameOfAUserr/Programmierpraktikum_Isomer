from django.db import models

# Create your models here.

class Substances(models.Model):
    ID = models.AutoField(primary_key=True, editable=False)
    names = models.JSONField(default=list, blank=True,null=True)
    formula = models.CharField(max_length=1000, blank=True,null=True)
    iupac_names = models.JSONField(default=list, blank=True,null=True)
    smiles = models.CharField(max_length=1000, blank=True,null=True)
    inchi = models.CharField(max_length=1000, blank=True,null=True)
    inchi_key = models.CharField(max_length=1000, blank=True,null=True)
    molecular_mass = models.FloatField(blank=True,null=True)
    cas_num = models.CharField(max_length=1000, default="0", blank=True, null=True)
    categories = models.JSONField(default=list, blank=True,null=True)
    source = models.JSONField(blank=True,null=True)
    validated = models.BooleanField(default=False, blank=True,null=True)
    deleted = models.BooleanField(default=False)
    last_modified = models.CharField(max_length=200, default="", blank=True,null=True)
    version = models.CharField(max_length=200, default="0", blank=True,null=True)
    details = models.JSONField(blank=True,null=True)

