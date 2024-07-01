from .models import Substances

class Store:
    def Substances(self, list_of_substances):
        """
        Process a list of substance dictionaries and store them in the database.
        If a substance with the same source, names, and formula exists, it is deleted before creating a new record.
        """
        for substance in list_of_substances:
            try:
                data_dict = substance  # Get the dictionary representing the substance
                
                # Extract the necessary fields from the dictionary
                data_source = data_dict["source"]
                data_names = data_dict["names"]
                data_formula = data_dict["formula"]
                
                # Check if any existing objects match the source, names, and formula
                existing_objects = Substances.objects.filter(source=data_source, names=data_names, formula=data_formula).all()
                print(existing_objects)
                
                if existing_objects:
                    # Delete existing objects if found
                    existing_objects.delete()
                    print("existing object deleted")
                
                try:
                    # Create a dictionary with validated data, using defaults for missing values
                    validated_data_dict = {
                        'smiles': data_dict.get('smiles', ''),
                        'names': data_dict.get('names', []),
                        'iupac_names': data_dict.get('iupac_names', []),
                        'formula': data_dict.get('formula', ''),
                        'inchi': data_dict.get('inchi', ''),
                        'inchi_key': data_dict.get('inchi_key', ''),
                        'molecular_mass': data_dict.get('molecular_mass', 0.0),
                        'cas_num': data_dict.get('cas_num', 'unknown'),
                        'categories': data_dict.get('categories', []),
                        'source': data_dict.get('source', {}),
                        'validated': data_dict.get('validated', False),
                        'deleted': data_dict.get('deleted', False),
                        'details': data_dict.get('details', {}),
                        'last_modified': data_dict.get('last_modified', ''),
                        'version': data_dict.get('version', '0'),
                    }
                    
                    # Create a new Substances object with the validated data
                    Substances.objects.create(
                        smiles=validated_data_dict['smiles'],
                        names=validated_data_dict['names'],
                        iupac_names=validated_data_dict['iupac_names'],
                        formula=validated_data_dict['formula'],
                        inchi=validated_data_dict['inchi'],
                        inchi_key=validated_data_dict['inchi_key'],
                        molecular_mass=validated_data_dict['molecular_mass'],
                        cas_num=validated_data_dict['cas_num'],
                        categories=validated_data_dict['categories'],
                        source=validated_data_dict["source"],
                        validated=validated_data_dict['validated'],
                        deleted=validated_data_dict['deleted'],
                        details=validated_data_dict['details'],
                        last_modified=validated_data_dict['last_modified']
                    )
                    print(validated_data_dict["smiles"])
                except Exception as e:
                    # Log any errors that occur during the creation of the Substances object
                    print(f"An error occurred: {e}")
                                
            except Exception as e:
                # Log any errors that occur during the processing of the substance
                print(f"An error occurred: {e}")
