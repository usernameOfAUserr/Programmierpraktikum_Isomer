import asyncio
import datetime
import aiohttp
import json
import os
import time
from bs4 import BeautifulSoup
import re
from .new_getCategory import getCategorys
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .models import Substances
from .Store import Store
from datetime import datetime, timezone

class getData:
    html_errors = 0  # Counter for HTML errors
    Storer = Store()  # Instance of Store class for storing data
    progress = 0  # Progress indicator
    gatherd_substances = []  # List to collect gathered substances
    Substances = Substances.objects.all()  # Fetch all substances from the database
    find_id_regex = re.compile(r'&id=\d{1,5}')  # Regex to find IDs in URLs
    just_id_regex = re.compile(r'\d{1,5}')  # Regex to extract just the ID
    categorys = {}  # Dictionary to store categories

    def __init__(self) -> None:
        pass

    def getProgress(self):
        """Return the current progress."""
        return self.progress

    async def fetch_url(self, session, url, index):
        """Fetch data from a URL and process it."""
        try:
            async with session.get(url) as response:
                if response.status != 200:
                    return True

                if index % 100 == 0:
                    self.progress = index / 150
                    print(f"Done for {url}")

                html = await response.text()
                data = BeautifulSoup(html, "html.parser")

                # Extract information from the data
                all_clippable = data.findAll(class_='clippable')
                number_of_clippable = len(all_clippable)

                # Extract names and IUPAC
                try:
                    name_list = data.findAll('ul', class_='name-list')
                    names = []
                    iupac = []

                    if len(name_list) == 2:
                        for name in name_list[0]:
                            if name.name == "li":
                                names.append(name.text.strip())
                        for name in name_list[1]:
                            if name.name == "li":
                                iupac.append(name.text.strip())
                    elif len(name_list) == 1:
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
                    print(f"Error by URL {url}, data quantity: {number_of_clippable}")
                    self.html_errors += 1

                # Extract index and category
                category = []
                id = data.find(id='substance-id').text
                if id in self.categorys.keys():
                    category = self.categorys[id]

                if len(category) == 0:
                    category.append("unknown")

                # Extract formula, molecular weight, InChI, InChI key, and SMILES
                formular = all_clippable[number_of_clippable-5].text
                molecular_weight = float(all_clippable[number_of_clippable-4].text)
                InChl_regex = re.compile(r"InChI=(.*)")
                match = InChl_regex.search(all_clippable[number_of_clippable-3].text)
                inchi = match.group(1)
                inchi_key = all_clippable[number_of_clippable-2].text
                smiles = data.find(id='smiles').text

                # Get current datetime in UTC
                current_datetime = datetime.now(timezone.utc)
                formatted_string = current_datetime.strftime('%Y-%m-%dT%H:%M:%SZ')

                # Validate data
                is_valid = True
                molecule = Chem.MolFromSmiles(smiles)
                canonical_smiles = Chem.MolToSmiles(molecule)
                weight_chem = Descriptors.MolWt(molecule)
                formula_chem = rdMolDescriptors.CalcMolFormula(molecule)
                if molecular_weight < round(weight_chem, 2) - 0.1 or molecular_weight > round(weight_chem, 2) + 0.1:
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
                    "last_modified": formatted_string,
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
            print(f"Error fetching URL {url}: {e}")
            return True

    async def get_responses(self, urls):
        """Fetch data from multiple URLs asynchronously."""
        try:
            timeout = aiohttp.ClientTimeout(total=None, connect=None, sock_connect=None, sock_read=60)
            async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=30), timeout=timeout) as session:
                tasks = [self.fetch_url(session, url, index) for index, url in enumerate(urls, 1)]
                await asyncio.gather(*tasks)
        except aiohttp.ClientConnectionError:
            print("Connection was dropped before finishing")
        except aiohttp.ServerTimeoutError:
            print("Timeout error")

    def start(self, urls):
        """Start the data fetching process."""
        self.progress = 0
        if urls is None:
            urls = [f"https://isomerdesign.com/pihkal/explore/{i}" for i in range(1, 17000)]
        reuse_old_categorys = False
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
        self.Storer.Substances(self.gatherd_substances)
        print("Errors: " + str(self.html_errors))
