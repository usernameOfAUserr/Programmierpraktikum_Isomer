import asyncio
import json
import aiohttp
from bs4 import BeautifulSoup
import re
import time
import requests

# Dictionary to store id and category associations
id_category_dictionary = {}

# List of URLs and their corresponding main categories to scrape
parameter_for_category_scrape = [
    ["https://isomerdesign.com/pihkal/browse/table/si", "Shulgin Index"],
    ["https://isomerdesign.com/pihkal/browse/table/pea", "Phenethylamine"],
    ["https://isomerdesign.com/pihkal/browse/collection/arylcycloalkylamine", "Arylcycloalkylamine"],
    ["https://isomerdesign.com/pihkal/browse/collection/aryldiazepine", "Aryldiazepine"],
    ["https://isomerdesign.com/pihkal/browse/collection/cannabinoid", "Cannabinoid"],
    ["https://isomerdesign.com/pihkal/browse/collection/fentanyl", "Fentanyl"]
]

# Function to extract information form the HTML of a category pade and update the id_category_dictionary
def extractInfoFromLink(text, MainCategory):
    """
    Extracts information from the HTML to map IDs to categories and subcategories.
    """
    id_regex = re.compile(r'\d+') # Regular expression to match numerical IDs
    option_element = text.find('option', selected=True)  # Get the currently selected subcategory
    if option_element:  # If a subcategory is present (e.g. Pent... or Shulgin Index)
        text_of_option_element = option_element.text # Name of the subcategory
        subcategory = text_of_option_element.split()[1]
        allLinks = text.find_all('a')
        for link in allLinks:
            if link.text.strip() == "Explore": # If the link text is "Explore"
                id_match = id_regex.search(link.get("href")) # Find the ID in the href attribute
                id = id_match.group() # Extract to match ID
                id_category_dictionary[id] = [MainCategory, subcategory] # Update dictionary with ID, main category, and subcategory
    else:  # If no subcategories are present
        allLinks = text.find_all('a')
        for link in allLinks:
            if link.text.strip() == "Explore": # If the link text is "Explore"
                id_match = id_regex.search(link.get("href")) # Find the ID in the href attribute
                id = id_match.group() # Extract to match ID
                id_category_dictionary[id] = [MainCategory] # Update dictionary with ID and main category
    return True

# Asynchronus function to scrape all URLs in a category
async def getAllUrlsOfCategorys(url, MainCategory, session): # session is needed for the async request
    """
    Recursively scrapes all URLs of a category to extract substance information.
    """
    print(f"{url}: {MainCategory}")
    async with session.get(url) as response: # Perform an async GET request
        html = await response.text() # Get the response text
        text = BeautifulSoup(html, "html.parser") # Parse the HTML with BeautifulSoup
        extractInfoFromLink(text, MainCategory) # Exract relevant data from the HTML
        all_links = text.find_all("a")
        for link in all_links:
            if link.text.strip() == "Next" and link.get("href") != '#': # If the link text is "Next" and href is not "#"
                if link.get("href") != "https://isomerdesign.com/pihkal/browse/table/si/5/352": # Specific URL check
                    await getAllUrlsOfCategorys(link.get("href"), MainCategory, session)
                else:
                    await getAllUrlsOfCategorys("https://isomerdesign.com/pihkal/browse/table/si/5/353", MainCategory, session)
        return True

#Asynchronous function to start scraping all defined categories 
async def start_category_scrape():

    timeout = aiohttp.ClientTimeout(total=None, connect=None, sock_connect=None, sock_read=60)
    async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=30), timeout=timeout) as session:
        tasks = []
        for parameter in parameter_for_category_scrape: # Loop trough all category parameters
            tasks.append(getAllUrlsOfCategorys(parameter[0], parameter[1], session))
        await asyncio.gather(*tasks) # Gather all tasks and run them concurrently

# Entry point function to start the scraping process
def getCategorys():
    """
    Orchestrates the scraping process and saves results to a file.
    """
    start = time.time()
    asyncio.run(start_category_scrape()) # Start the program that scrapes all the categories and looks which substances are listed there
    print(f"Runtime:  {int(time.time() - start)} seconds")
    print(f"{len(id_category_dictionary.keys())} ids found")
    with open("category_id_json_file.json", "w") as file:
        json.dump(id_category_dictionary, file) # Save the results to a JSON file
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
