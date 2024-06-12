import asyncio
import aiohttp
import json
import os
import time
import re
from periodictable import formula
from .getCategory import getCategorys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .models import Substances
from django.apps import apps
from .Store import Store


class Exchanger:
    Storer = Store()
    all_databases = apps.get_models()

    def generate(self, which_database):
        which_database = "Substances"
        target_db = None
        for db in self.all_databases:
            if db.__name__ == which_database:
                target_db = db
                print(target_db)
                break
        if target_db is None:
            return False
        
        substances_in_db = target_db.objects.all()
        first_sub = target_db.objects.first()
        keys = first_sub.__dict__.keys()
        keys_to_remove = ['_state', 'id']
        usefull_keys = [key for key in keys if key not in keys_to_remove]
        i = 0
        print(usefull_keys)
        data_for_file = []
        for sub in substances_in_db:
            sub_data = {}
            for key in usefull_keys:
                if key != "last_modified":
                    sub_data[key] = getattr(sub, key)
                else:
                     sub_data[key] = str(getattr(sub, key))
            data_for_file.append(sub_data)
            i+=1
         
         
        print(str(i)+" Substances written to File")
        self.store_as_file(data_for_file=data_for_file)
        absolute_path = os.path.abspath("PIHKAL.json")
        print("Absoluter Pfad:", absolute_path)
        return absolute_path

    def store_as_file(self, data_for_file):
        with open("PIHKAL.json", "w") as f:  
            json.dump(data_for_file, f)

    def process(self,file):
        self.Storer.Substances(file)

        return True