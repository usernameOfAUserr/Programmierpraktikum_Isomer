import asyncio
import json
import aiohttp
from bs4 import BeautifulSoup
import re
import time

id_category_dictionary = {}

parameter_for_category_scrape = [
    ["https://isomerdesign.com/pihkal/browse/table/si", "Shulgin Index"],
    ["https://isomerdesign.com/pihkal/browse/table/pea", "Phenethylamine"],
    ["https://isomerdesign.com/pihkal/browse/collection/arylcycloalkylamine", "Arylcycloalkylamine"],
    ["https://isomerdesign.com/pihkal/browse/collection/aryldiazepine", "Aryldiazepine"],
    ["https://isomerdesign.com/pihkal/browse/collection/cannabinoid", "Cannabinoid"],
    ["https://isomerdesign.com/pihkal/browse/collection/fentanyl", "Fentanyl"]
]


def extractInfoFromLink(text, MainCategory):
    id_regex = re.compile(r'\d+')
    option_element = text.find('option', selected=True)  # gets the current sub_category
    if (option_element):  # is there when Pent.. or Shulgin Index
        text_of_option_element = option_element.text  # the name of the sub_category
        subcategory = text_of_option_element.split()[1]
        allLinks = text.find_all('a')
        for link in allLinks:
            if link.text.strip() == "Explore":
                id_uebereinstimmung = id_regex.search(link.get("href"))
                id = id_uebereinstimmung.group()
                id_category_dictionary[id] = [MainCategory, subcategory]
    else:  # there are no subcategorys
        allLinks = text.find_all('a')
        for link in allLinks:
            if link.text.strip() == "Explore":
                id_uebereinstimmung = id_regex.search(link.get("href"))
                id = id_uebereinstimmung.group()
                id_category_dictionary[id] = [MainCategory]
    return True


async def getAllUrlsOfCategorys(url, MainCategory, session,already_scraped_urls):  # session is needed for the async request

    print(url+": "+MainCategory)
    async with session.get(url) as response:
        html = await response.text()
        text = BeautifulSoup(html, "html.parser")
        extractInfoFromLink(text, MainCategory)  # to get all the relevent data from the html text
        already_scraped_urls.append(url)
        all_links = text.find_all("a")
        for link in all_links:
            if link.text.strip() == "Next" and link.get("href") != '#':
                if link.get("href") != "https://isomerdesign.com/pihkal/browse/table/si/5/352":
                    await getAllUrlsOfCategorys(link.get("href"), MainCategory, session,already_scraped_urls)
                else:
                    await getAllUrlsOfCategorys("https://isomerdesign.com/pihkal/browse/table/si/5/353", MainCategory, session,already_scraped_urls)

        return True


async def start_category_scrape():
    already_scraped_urls = []
    timeout = aiohttp.ClientTimeout(total=None, connect=None, sock_connect=None,sock_read=60)
    async with aiohttp.ClientSession(connector=aiohttp.TCPConnector(limit=30),timeout=timeout) as session:
        tasks = []
        for parameter in parameter_for_category_scrape:
            tasks.append(getAllUrlsOfCategorys(parameter[0], parameter[1], session,already_scraped_urls))
        await asyncio.gather(*tasks)


def getCategorys():  # thats the entry point
    start = time.time()
    already_scraped_urls = []
    asyncio.run(start_category_scrape())  # start the programm that scrapes all the categorys and looks which substances are listed there
    print(f"Runtime:  {int(time.time() - start)}")
    print(str(len(id_category_dictionary.keys())) + " ids found")
    with open("cateogry_id_json_file.json", "w") as file:
        json.dump(id_category_dictionary, file)
    return id_category_dictionary

