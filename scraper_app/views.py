import json
import os
import re
import aiohttp
import asyncio
from django.shortcuts import render
from django.http import HttpResponse, JsonResponse
from scraper_app.models import Substances
from .getData import getData
from .getReachableUrls import getReachableUrls
from .Exchanger import Exchanger

# Initial Setup
latest_search_results = []
getDataObject = getData()
JS = Substances.objects.all()
keys = [field.name for field in Substances._meta.fields]
keys.remove("ID")
model_fields_that_are_lists = ["names", "details", "iupac_names", "categories"]

def scraper(request):
    """Handles search requests and displays results."""
    if request.method == "POST":
        requested = []
        selected = request.POST.get('select')
        searched = request.POST.get('search_field')

        if selected == "molecular_mass":
            toleranz = request.POST.get("abweichung")
            start = float(searched) - float(toleranz)
            end = float(searched) + float(toleranz)
            requested = Substances.objects.filter(molecular_mass__range=[start, end])
        elif selected == "smiles":
            requested = Substances.objects.filter(smiles=searched)
        elif selected == "formula":
            format = ".*" + searched + ".*"
            compare_regex = re.compile(format)
            requested = [formula for formula in JS if compare_regex.search(formula.formula)]
        elif selected == "names":
            requested = [names for names in JS if searched in names.names]
        elif selected == "iupac_names":
            requested = [names for names in JS if searched in names.iupac_names]
        elif selected == "cas_num":
            requested = Substances.objects.filter(cas_num=searched)
        elif selected == "categories":
            requested = [
                {
                    "names": res.names,
                    "formula": res.formula,
                    "iupac_names": res.iupac_names,
                    "smiles": res.smiles,
                    "inchi": res.inchi,
                    "inchi_key": res.inchi_key,
                    "molecular_mass": res.molecular_mass,
                    "cas_num": res.cas_num,
                    "categories": res.categories,
                    "source": res.source,
                    "validated": res.validated,
                    "deleted": res.deleted,
                    "last_modified": res.last_modified,
                    "version": res.version,
                    "details": res.details,
                }
                for res in Substances.objects.filter(categories__icontains=searched)
            ]
        elif selected == "source":
            requested = [req for req in Substances.objects.filter(source__icontains=searched)]
        elif selected == "validated":
            requested = Substances.objects.filter(valid=searched)
        elif selected == "deleted":
            requested = Substances.objects.filter(deleted=searched)
        elif selected == "version":
            requested = Substances.objects.filter(version=searched)
        elif selected == "details":
            requested = Substances.objects.filter(details__icontains=searched)
        else:
            requested = None

        global latest_search_results
        latest_search_results = requested

        return render(request, "scraper.html", {
            "answer": requested,
            "selected": selected,
            "keys": keys,
        })
    else:
        db = Substances.objects.all()
        return render(request, "scraper.html", {
            "db": db,
            "answer": None,
            "keys": keys,
        })

def reset_database(request):
    """Resets the database."""
    getDataObject.start(None)
    return HttpResponse("Database was fully restored")

def request_how_many_json_file(request):
    """Returns the progress of JSON file processing."""
    progress = getDataObject.getProgress()
    return JsonResponse({"progress": progress})

def delete_search_results(request):
    """Deletes search results from the database."""
    global latest_search_results
    print("Delete request for " + str(len(latest_search_results)) + " substances")
    i = 0
    for sub in latest_search_results:
        print(sub)
        Substances.objects.filter(smiles=sub.smiles).delete()
        print(str(i) + " deleted")
        i += 1
    latest_search_results = []
    return HttpResponse(200)

def search_for_newcomers(request):
    """Finds and returns new substances."""
    new_ones = getReachableUrls()
    print("stored in db")
    return JsonResponse({"newSubstances": new_ones})

async def get_witz_asynchron(url):
    """Asynchronously fetches a joke."""
    async with aiohttp.ClientSession() as session:
        try:
            async with session.get(url) as response:
                data = await response.json()
                witz = data[0]["text"]
                return witz
        except Exception as e:
            return {"error": str(e)}

def get_witz(request):
    """Fetches a joke and returns it as JSON."""
    url = "https://witzapi.de/api/joke/"
    witz = asyncio.run(get_witz_asynchron(url))
    return JsonResponse({"witz": witz})

def generateJsonFile(request):
    """Generates a JSON file and returns its path."""
    exchanger = Exchanger()
    path_to_file = exchanger.generate("example_name")
    return HttpResponse(path_to_file)

def processJsonInput(request):
    """Processes JSON files from a folder."""
    how_many_files = 0
    folder_path = "scraper_app/files_to_insert"

    for filename in os.listdir(folder_path):
        how_many_files += 1
        path_to_file = os.path.join(folder_path, filename)
        with open(path_to_file, "r") as f:
            data_to_insert = json.load(f)
            exchanger = Exchanger()
            exchanger.process(data_to_insert)
    
    if how_many_files > 0:
        return HttpResponse(f"{how_many_files} Files successfully processed and stored")
    
    abs_path = os.path.abspath(folder_path)
    return HttpResponse(f"There are no files in \n{abs_path}\n add one in the directory")

def search_suggestion(request):
    """Provides search suggestions based on the category and searched term."""
    category = request.GET.get('category')
    searched = request.GET.get('searched')

    print(f"Category: {category}, Searched: {searched}")

    requested = []
    most_relevant = []

    if category == "smiles":
        requested = Substances.objects.filter(smiles__icontains=searched)
        for req in requested.values():
            if category in model_fields_that_are_lists:
                for li in req["smiles"]:
                    if li.startswith(searched):
                        most_relevant.append(li)
            else:
                if req[category].startswith(searched):
                    most_relevant.append(req[category])
    elif category == "formula":
        requested = Substances.objects.filter(formula__icontains=searched)
        for req in requested.values():
            if category in model_fields_that_are_lists:
                for li in req["formula"]:
                    if li.startswith(searched):
                        most_relevant.append(li)
            else:
                if req["formula"].startswith(searched):
                    most_relevant.append(req[category])
    elif category == "names":
        requested = Substances.objects.filter(names__icontains=searched)
        for req in requested.values():
            if category in model_fields_that_are_lists:
                for li in req["names"]:
                    if li.startswith(searched):
                        most_relevant.append(li)
            else:
                if req["names"].startswith(searched):
                    most_relevant.append(req[category])
    elif category == "iupac_names":
        requested = Substances.objects.filter(iupac_names__icontains=searched)
        for req in requested.values():
            if category in model_fields_that_are_lists:
                for li in req["iupac_names"]:
                    if li.startswith(searched):
                        most_relevant.append(li)
            else:
                if req["iupac_names"].startswith(searched):
                    most_relevant.append(req[category])
    elif category == "cas_num":
        requested = Substances.objects.filter(cas_num__icontains=searched)
        for req in requested.values():
            if category in model_fields_that_are_lists:
                for li in req["cas_num"]:
                    if li.startswith(searched):
                        most_relevant.append(li)
            else:
                if req["cas_num"].startswith(searched):
                    most_relevant.append(req[category])
    elif category == "categories":
        requested = Substances.objects.filter(categories__icontains=searched)
        for req in requested.values():
            if req["categories"] not in most_relevant:
                most_relevant.append(req["categories"])
    elif category == "source":
        requested = Substances.objects.filter(source__icontains=searched)
        for req in requested:
            if req.source["name"] not in most_relevant:
                most_relevant.append(req.source["name"])
    elif category == "validated":
        requested = Substances.objects.filter(validated__icontains=searched)
    elif category == "deleted":
        requested = Substances.objects.filter(deleted__icontains=searched)
    elif category == "version":
        requested = Substances.objects.filter(version__icontains=searched)
    elif category == "details":
        requested = Substances.objects.filter(details__icontains=searched)

    return JsonResponse(most_relevant, safe=False)
