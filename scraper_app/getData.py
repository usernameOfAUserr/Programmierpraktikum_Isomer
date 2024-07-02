import asyncio
import datetime
import aiohttp
import json
import os
import time
from bs4 import BeautifulSoup
import re
from periodictable import formula
from .getCategory import getCategorys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .models import Substances
from .Store import Store
from datetime import datetime, timezone

class getData:
    html_errors = 0
    Storer = Store()
    progress = 0
    gatherd_substances = []
    Substances = Substances.objects.all()
    find_id_regex = re.compile(r'&id=\d{1,5}')
    just_id_regex = re.compile(r'\d{1,5}')
    categorys = {}
    def __init__(self) -> None:
        pass

   
    def getProgress(self):
        return self.progress
    

    async def fetch_url(self, session, url, index):

        try:  # for the case that the uls isnt accessible
            async with session.get(url) as response:
                if response.status != 200:
                    return True 
        
                if index %100 == 0:
                    self.progress = index / 150
                    print(f"done for {url}")
                html = await response.text()  # wait until the server responses
                # parse to process data better
                data = BeautifulSoup(html, "html.parser")

                # extrac informations from data
                
                all_clippable = data.findAll(class_='clippable')  #so i can calculate exactly
                number_of_clippable = len(all_clippable)

                # names and iupac, can be tricky because sometimes, there are no names as well as no iupac, so it has to be calculated what exactly is there
                try:
                    name_list = data.findAll('ul', class_='name-list')
                    names = []
                    iupac = []
                    # names
                    if len(name_list) == 2:
                        for name in name_list[0]:
                            if name.name == "li":
                                names.append(name.text.strip())
                        # iupac
                        for name in name_list[1]:
                            if name.name == "li":
                                iupac.append(name.text.strip())
                    elif len(name_list) == 1:
                        # now, I calculate the position of iupac and names, in addition regex to pin down where the iupac begins
                        parent_of_name_list = name_list[0].find_parent().find_parent()
                        child_elements = parent_of_name_list.children
                        just_iupac_exists = False
                        for child in child_elements:
                            if child.text == "IUPAC names" or child.text == "IUPAC name":
                                just_iupac_exists = True
                                names.append("Unknown")
                                for name in name_list[0]:
                                    if name.name == "li":
                                        iupac.append(name.text.strip())
                        if not just_iupac_exists:
                            iupac.append("Unknown")
                            for name in name_list[0]:
                                if name.name == "li":
                                    names.append(name.text.strip())
                    else:
                        iupac.append("Unknown")
                        names.append("Unknown")
         
                except:
                    print(f"error by url {url}, data_quantity: {number_of_clippable}")
                    self.html_errors = self.html_errors + 1

                # index & category
                category = []
                # substance_index = index

                # category
                # look, if the index is known in categorys
                id = data.find(id='substance-id').text
                if id in self.categorys.keys():
                    category =self.categorys[id]
       
                if len(category) == 0:
                    category.append("unknown")
                # formular
                number_of_clippable = len(data.findAll(class_='clippable'))
                formular = all_clippable[number_of_clippable-5].text
                # molecular_weight
                molecular_weight = all_clippable[number_of_clippable-4].text
                molecular_weight = float(molecular_weight)
                # inChl
                InChl_regex = re.compile(r"InChI=(.*)")
                match = InChl_regex.search(all_clippable[number_of_clippable-3].text)
                #to remove the inchi= 
                inchi = match.group(1)
                # inchl_key
                inchi_key = all_clippable[number_of_clippable-2].text
                #smiles
                smiles = data.find(id='smiles').text
                #datetime
                current_datetime = datetime.now(timezone.utc)

                # Formatieren des datetime-Objekts in das gewünschte Format
                formatted_string = current_datetime.strftime('%Y-%m-%dT%H:%M:%SZ')
            
                #validation
                is_valid = True
                molecule = Chem.MolFromSmiles(smiles)
                canonical_smiles = Chem.MolToSmiles(molecule)
                weight_chem = Descriptors.MolWt(molecule) 
                formula_chem = rdMolDescriptors.CalcMolFormula(molecule)
                if molecular_weight < round(weight_chem,2) - 0.1 or molecular_weight > round(weight_chem, 2)+ 0.1:
                    is_valid = False
                if formula_chem != formular:
                    is_valid = False
                data_dict = {
                    "smiles": canonical_smiles,
                    "names": names,
                    "iupac_names": iupac,
                    "formula": formular,
                    "inchi": inchi,
                    "inchi_key": inchi_key,
                    "molecular_mass": molecular_weight,
                    "cas_num": 0,
                    "last_modified":formatted_string,
                    "categories": category,
                    "source": {"name": "PIHKAL", "url": url},
                    "validated": is_valid,
                    "deleted": False,
                    "version": "0.0",
                    "details": {},
                }
                self.gatherd_substances.append(data_dict)
                return True
        except Exception as e:
            print(f"Fehler beim Abrufen der URL {url}: {e}")
            return True



    async def get_responses(self, urls):

        try:
            timeout = aiohttp.ClientTimeout(total=None, connect=None, sock_connect=None, #total timelimit for the whole reques /

                                            sock_read=60)  # Timeout von 60 Sekunden für Socken lesen
            async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=30),timeout=timeout) as session: #to limit the requests in order to respekt the server """connector=aiohttp.TCPConnector(limit=10)"""
                    # creates ClientSession-Object to manage the asynchron Communication; within the async with-block, get-requests are done; clears after himself when done
                    tasks = [self.fetch_url(session, url, index) for index, url in
                            enumerate(urls, 1)]  # create a list of tasks, where every item is a call of the fetch_url function
                    responses = await asyncio.gather(
                        *tasks)  # wait until all tasks(calls of the fetch_url-function) came to a result
        except aiohttp.ClientConnectionError:
            print("Conntection was dropped before finishing")
        except aiohttp.ServerTimeoutError:
            print("Timeouterror")

    def start(self, urls, testing=False):
        self.progress = 0
        self.gatherd_substances = []
        if urls is None:
            urls = [f"https://isomerdesign.com/pihkal/explore/{i}" for i in range(1, 17000)]
        reuse_old_categorys = True  #when this is True, then categorys from the last time they were scraped are used
        if reuse_old_categorys:
            with open("category_id_json_file.json", 'r') as f:
                self.categorys = json.load(f)
        else:
            try:
                self.categorys = getCategorys()
            except:
                with open("category_id_json_file.json", 'r') as f:
                    self.categorys = json.load(f)
        
        # Load the categories from an external file to improve speed
        asyncio.run(self.get_responses(urls))
        if testing is False:
            self.Storer.Substances(self.gatherd_substances)
            print("Errors: " + str(self.html_errors))
        else:
            return self.gatherd_substances


    def test_getData(self):  #is used in views.py
        try:
            # Create instances of getData and Store
            getData_Object = getData()
            Store_Object = Store()
            
            # Gather data from the URL
            gathered_substance_list = getData_Object.start(["https://isomerdesign.com/pihkal/explore/1"], testing=True)
            gathered_substance = gathered_substance_list[0]
            # Check if gathered_substance is a dictionary
            if not isinstance(gathered_substance, dict):
                print("getData: Error -> gathered_substance is not a dictionary")
                return
            # Process the gathered substance
            substance_query_set = Store_Object.Substances([gathered_substance], testing=True)
            if len(substance_query_set) == 1:
                substance = substance_query_set.first()
            # Extract values
            names = substance.names                # List of names (JSONField)
            formula = substance.formula            # Formula (CharField)
            iupac_names = substance.iupac_names    # List of IUPAC names (JSONField)
            smiles = substance.smiles              # SMILES string (CharField)
            inchi = substance.inchi                # InChI string (CharField)

            # Perform checks and print results
            if names != gathered_substance.get('names', []):
                print(f"Store: Failed -> Names do not match: {names} != {gathered_substance.get('names', [])}")
            elif formula != gathered_substance.get('formula', ''):
                print(f"Store: Failed -> Formula does not match: {formula} != {gathered_substance.get('formula', '')}")
            elif iupac_names != gathered_substance.get('iupac_names', []):
                print(f"Store: Failed -> IUPAC Names do not match: {iupac_names} != {gathered_substance.get('iupac_names', [])}")
            elif smiles != gathered_substance.get('smiles', ''):
                print(f"Store: Failed -> SMILES do not match: {smiles} != {gathered_substance.get('smiles', '')}")
            elif inchi != gathered_substance.get('inchi', ''):
                print(f"Store: Failed -> InChI does not match: {inchi} != {gathered_substance.get('inchi', '')}")
            else:
                print("getData: Successful\nStore: Successful")
        
        except Exception as e:
            print(f"getData: Failed -> {str(e)}")