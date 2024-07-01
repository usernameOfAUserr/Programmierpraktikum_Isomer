import asyncio
import aiohttp
import json
import os
import re
from periodictable import formula
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .models import Substances  # Importing the Substances model
from django.apps import apps  # Importing apps to get all models in the project
from .Store import Store  # Importing the Store class

class Exchanger:
    Storer = Store()  # Instance of Store class for storing data
    all_databases = apps.get_models()  # Fetch all models in the project

    def generate(self, which_database):
        """Generate JSON file from the specified database."""
        which_database = "Substances"  # Hardcoded to "Substances" for this example
        target_db = None

        # Find the target database model
        for db in self.all_databases:
            if db.__name__ == which_database:
                target_db = db
                print(target_db)
                break

        if target_db is None:
            return False

        # Fetch all records and keys from the target database
        substances_in_db = target_db.objects.all()
        first_sub = target_db.objects.first()
        keys = first_sub.__dict__.keys()

        # Define keys to remove from the data
        keys_to_remove = ['_state', 'id', 'ID']
        usefull_keys = [key for key in keys if key not in keys_to_remove]

        i = 0
        print(usefull_keys)
        data_for_file = []

        # Extract useful data for each substance
        for sub in substances_in_db:
            sub_data = {}
            for key in usefull_keys:
                if key != "last_modified":
                    sub_data[key] = getattr(sub, key)
                else:
                    sub_data[key] = str(getattr(sub, key))  # Convert datetime to string
            data_for_file.append(sub_data)
            i += 1

        print(f"{i} Substances written to File")
        self.store_as_file(data_for_file=data_for_file)
        absolute_path = os.path.abspath("PIHKAL.json")
        print("Absoluter Pfad:", absolute_path)
        return absolute_path

    def store_as_file(self, data_for_file):
        """Store the extracted data in a JSON file."""
        with open("PIHKAL.json", "w") as f:
            json.dump(data_for_file, f)

    def process(self, file):
        """Process the given file and store the data."""
        self.Storer.Substances(file)
        return True
