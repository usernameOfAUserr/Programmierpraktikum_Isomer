import asyncio
import aiohttp
import json
import re
import os
from bs4 import BeautifulSoup
from periodictable import formula
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from .models import Substances
from .getData import getData

# Global variables to store various states and results
changes_accoure = []  # List to keep track of URLs where changes occurred
not_reachable = []  # List to keep track of URLs that are not reachable
counter = 0  # Counter to track the number of processed URLs

# Regular expressions to find IDs in URLs
find_id_regex = re.compile(r'&id=\d{1,5}')
just_id_regex = re.compile(r'\d{1,5}')

# Objects and lists to manage data and categories
previous_categorys = {}
gathered_substances = []
getDataObject = getData()

async def find_not_reachable_urls(session, url):
    """Check if a URL is reachable and update the respective lists."""
    global counter
    counter += 1
    try:
        # Print progress every 100 URLs
        if counter % 100 == 0:
            print(counter)
        # Make a HEAD request to check if the URL is reachable
        async with session.head(url) as response:
            status_code = response.status
            if status_code >= 400:
                not_reachable.append(url)
            else:
                changes_accoure.append(url)
    except:
        print("Couldn't make connection")
        not_reachable.append(url)

async def manage(urls):
    """Manage the checking of URL accessibility."""
    print(f"Length of urls: {len(urls)}")
    timeout = aiohttp.ClientTimeout(total=None, connect=None, sock_connect=None, sock_read=60)
    # Create a session with a TCP connector and a timeout
    async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=30), timeout=timeout) as session:
        tasks = [find_not_reachable_urls(session, url) for url in urls]
        await asyncio.gather(*tasks)

def getReachableUrls():
    """Determine which URLs are reachable and update the not_reachable_urls.json file."""
    global not_reachable, changes_accoure, counter
    counter = 0  # Reset counter
    try:
        # Load previously not reachable URLs from a JSON file
        with open("not_reachable_urls.json", "r") as json_file:
            file = json.load(json_file)
        previous_not_reachable = file["not_reachable"]
        
        print("Not reachable URLs loaded")
    
        # Check if previously not reachable URLs are still not reachable
        asyncio.run(manage(previous_not_reachable))
 
    except:
        print("Previous not reachable URLs couldn't be loaded")
        # Generate a list of URLs to check if JSON file doesn't exist
        urls = [f"https://isomerdesign.com/pihkal/explore/{i}" for i in range(1, 17000)]
        asyncio.run(manage(urls))
        print("URL accessibility determined")

    # Save the updated list of not reachable URLs
    data = {
        "not_reachable": not_reachable
    }
    with open("not_reachable_urls.json", "w") as file:
        json.dump(data, file)
   
    print(str(changes_accoure))
    if changes_accoure:
        print("Changes occurred")
        getDataObject.start(changes_accoure)
    
    return changes_accoure
