import datetime
from .models import Substances
import time
class Store():
     def Substances(self, list_of_substances):
        for substance in list_of_substances:
            try:
                data_dict = substance
                data_source = data_dict["source"]
                existing_objects = Substances.objects.filter(source=data_source).all()
                if existing_objects:
                    existing_objects.delete()
                    print("existing object deleted")
              
                Substances.objects.create(
                    smiles=data_dict['smiles'],
                    names=data_dict['names'],
                    iupac_names=data_dict['iupac_names'],
                    formula=data_dict['formula'],
                    inchi=data_dict['inchi'],
                    inchi_key=data_dict['inchi_key'],
                    molecular_mass=data_dict['molecular_mass'],
                    cas_num=data_dict['cas_num'],
                    categories=data_dict['categories'],
                    source=data_dict["source"],
                    validated=data_dict['validated'],
                    deleted=data_dict['deleted'],
                    details=data_dict['details'],
                    last_modified=data_dict['last_modified']
                )

            except:
                print(str(substance['names'])+ " couldnt be stored")

    