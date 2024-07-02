import asyncio
import json
import aiohttp
from bs4 import BeautifulSoup
import re
import time
import requests

# Dictionary to store id and category associations
id_category_dictionary = {}

# Parameters defining the URLs and main categories to scrape
parameter_for_category_scrape = [
    ["https://isomerdesign.com/pihkal/browse/table/si", "Shulgin Index"],
    ["https://isomerdesign.com/pihkal/browse/table/pea", "Phenethylamine"],
    ["https://isomerdesign.com/pihkal/browse/collection/arylcycloalkylamine", "Arylcycloalkylamine"],
    ["https://isomerdesign.com/pihkal/browse/collection/aryldiazepine", "Aryldiazepine"],
    ["https://isomerdesign.com/pihkal/browse/collection/cannabinoid", "Cannabinoid"],
    ["https://isomerdesign.com/pihkal/browse/collection/fentanyl", "Fentanyl"]
]

def extractInfoFromLink(text, MainCategory):
    """
    Extracts information from the HTML to map IDs to categories and subcategories.
    """
    id_regex = re.compile(r'\d+')
    option_element = text.find('option', selected=True)  # Find the selected option element for subcategory
    if option_element:  # If a subcategory exists
        text_of_option_element = option_element.text
        subcategory = text_of_option_element.split()[1]
        allLinks = text.find_all('a')
        for link in allLinks:
            if link.text.strip() == "Explore":
                id_match = id_regex.search(link.get("href"))
                id = id_match.group()
                id_category_dictionary[id] = [MainCategory, subcategory]
    else:  # No subcategories
        allLinks = text.find_all('a')
        for link in allLinks:
            if link.text.strip() == "Explore":
                id_match = id_regex.search(link.get("href"))
                id = id_match.group()
                id_category_dictionary[id] = [MainCategory]
    return True

async def getAllUrlsOfCategorys(url, MainCategory, session):
    """
    Recursively scrapes all URLs of a category to extract substance information.
    """
    print(f"{url}: {MainCategory}")
    async with session.get(url) as response:
        html = await response.text()
        text = BeautifulSoup(html, "html.parser")
        extractInfoFromLink(text, MainCategory)
        all_links = text.find_all("a")
        for link in all_links:
            if link.text.strip() == "Next" and link.get("href") != '#':
                if link.get("href") != "https://isomerdesign.com/pihkal/browse/table/si/5/352":
                    await getAllUrlsOfCategorys(link.get("href"), MainCategory, session)
                else:
                    await getAllUrlsOfCategorys("https://isomerdesign.com/pihkal/browse/table/si/5/353", MainCategory, session)
        return True

async def start_category_scrape():
    """
    Initiates the asynchronous scraping of all defined categories.
    """
    
    timeout = aiohttp.ClientTimeout(total=None, connect=None, sock_connect=None, sock_read=60)
    async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=30), timeout=timeout) as session:
        tasks = []
        for parameter in parameter_for_category_scrape:
            tasks.append(getAllUrlsOfCategorys(parameter[0], parameter[1], session))
        await asyncio.gather(*tasks)

def getCategorys():
    """
    Entry point for the category scraping. Orchestrates the scraping process and saves results to a file.
    """
    start = time.time()
    asyncio.run(start_category_scrape())
    print(f"Runtime:  {int(time.time() - start)} seconds")
    print(f"{len(id_category_dictionary.keys())} ids found")
    with open("category_id_json_file.json", "w") as file:
        json.dump(id_category_dictionary, file)
    return id_category_dictionary

class test_getCategory_Class:

    def test_getCategorys(self):
        try:
            global id_category_dictionary
            id_category_dictionary.clear()
            
            # Load category IDs from file
            with open("category_id_json_file.json", "r") as file:
                categorys_id_from_file = json.load(file)
            
            # Fetch and parse HTML content
            fentanyl_html = requests.get("https://isomerdesign.com/pihkal/browse/collection/fentanyl")
            text = BeautifulSoup(fentanyl_html.text, "html.parser")
            
            # Extract information from link
            extractInfoFromLink(text, "Fentanyl")
            
            # Count previous categories
            count_how_many_previous = sum(1 for category in categorys_id_from_file.values() if category == ["Fentanyl"])
            
            # Count current categories
            count_how_many_current = sum(1 for category in id_category_dictionary.values() if category == ["Fentanyl"])
            
            # Check if counts are within the acceptable range
            if abs(count_how_many_current - count_how_many_previous) > 20:
                print(f"getCategorys: suspicious --> current_found: {count_how_many_current} vs. previous_found: {count_how_many_previous}")
            else:
                print("getCategorys: working")
        
        except Exception as e:
            print(f"getCategorys: failed ({str(e)})")
