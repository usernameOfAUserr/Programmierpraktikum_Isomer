import json
from .getData import getData
from django.shortcuts import render
import os
import logging
import requests
from scraper_app.models import Substances
from django.http import HttpResponse, JsonResponse
import time
import shutil
import re
from bs4 import BeautifulSoup
from .getReachableUrls import getReachableUrls
from django.template.defaultfilters import stringfilter
import aiohttp, asyncio
from .Exchanger import Exchanger

# Create your views here.
latest_search_results = []
getDataObject = getData()
JS = Substances.objects.all()
keys = [field.name for field in Substances._meta.fields]
keys.remove("ID")
print(keys)
model_fields_that_are_lists = ["names", "details", "iupac_names", "categories"]

def scraper(request):
    if request.method == "POST":
        requested = []
        selected = request.POST.get('select')
        searched = request.POST.get('search_field')
        if selected == "molecular_mass":
            toleranz = request.POST.get("abweichung")
            start = float(searched)-float(toleranz)
            end = float(searched)+float(toleranz)
            requested = Substances.objects.filter(molecular_mass__range=[start,end])
        elif selected == "smiles":
            if Substances.objects.filter(smiles=searched) is not None:
                requested = Substances.objects.filter(smiles=searched)
            else:
                request = None
        elif selected == "formula":
            requested = []
            format = ".*" + searched + ".*"
            compare_regex = re.compile(format)
            for formula in JS:
                if compare_regex.search(formula.formula) is not None:
                    requested.append(formula)
        elif selected == "names":
            requested = []
            for names in JS:
                if searched in names.names:
                    requested.append(names)
        elif selected == "iupac_nams":
            requested = []
            for names in JS:
                if searched in names.iupac_names:
                    requested.append(names)        
        elif selected == "cas_num":
            cas_num_results = Substances.objects.filter(cas_num=searched)

            requested.append(cas_num_results)

        elif selected == "categories":
            category_results = Substances.objects.filter(categories__icontains=searched)
            for res in category_results:
                requested.append({
                   
                    "names": res.names,
                    "formula": res.formula,
                    "iupac_names": res.iupac_names,
                    "smiles": res.smiles,
                    "inchi": res.inchi,
                    "inchi_key": res.inchi_key,
                    "molecular_mass": res.molecular_mass,
                    "cas_num": res.cas_num,
                    "categories": res.categories,
                    "source" :res.source,
                    "validated": res.validated,
                    "deleted": res.deleted,
                    "last_modified": res.last_modified,
                    "version": res.version,
                    "details": res.details,
                
                })

        elif selected == "source":
            requested = []
            print("searched: "+ searched+ " in categroy: " + selected)
            requested_substances = Substances.objects.filter(source__icontains=searched)
            for req in requested_substances:
                requested.append(req)
            print(len(requested))

        elif selected == "validated":
            requested = []
            valid_results = Substances.objects.filter(valid=searched)
            requested.append(valid_results)

        elif selected == "deleted":
            requested = []
            deleted_results = Substances.objects.filter(deleted=searched)
            requested.append(deleted_results)

        elif selected == "version":
            requested = []
            version_results = Substances.objects.filter(version=searched)
            requested.append(version_results)

        elif selected == "details":
            requested = []
            details_results = Substances.objects.filter(details__icontains=searched)
            requested.append(details_results)
    
        else:
            requested = None
        global latest_search_results
        latest_search_results = requested
        return render(request, "scraper.html",{
            "answer": requested, "selected":selected,
                        "keys":keys,

        })
    else:
        db = Substances.objects.all()
        return render(request,"scraper.html",{
            "db": db, "answer": None,
            "keys":keys,
        })


def reset_database(request):
    getDataObject.start(None)
    return HttpResponse("Database was fully restored")
    
def request_how_many_json_file(request):
    progress = getDataObject.getProgress()
    return JsonResponse({"progress":progress})
    
def delete_search_results(request):
    global latest_search_results
    print("Delete request for "+ str(len(latest_search_results)) + " substances")
    i = 0
    for sub in latest_search_results:
        print(sub)
        
        Substances.objects.filter(smiles=sub.smiles).delete()
        print(str(i) + " deleted")
        i+=1
    latest_search_results = []
    return HttpResponse(200)

def search_for_newcomers(request):
    new_ones = getReachableUrls()
    print("stroed in db")
    return JsonResponse({"newSubstances": new_ones})

################################################################
#Get joke

async def get_witz_asynchron(url):
    async with aiohttp.ClientSession() as session:  
        try:
            async with session.get(url) as response:
                data = await response.json()
                witz = data[0]["text"]
                return witz
        except Exception as e:  
            return {"error": str(e)}

async def get_witz_task(session, url):
    async with session.get(url) as response:
        data = await response.json()
        witz = data[0]["text"]
        return witz

def get_witz(request):
    url = "https://witzapi.de/api/joke/"
    witz = asyncio.run(get_witz_asynchron(url))  
    return JsonResponse({"witz": witz})
#####################################################################

def generateJsonFile(request):
    exchanger = Exchanger()
    path_to_file = exchanger.generate("example_name")
    return HttpResponse(path_to_file)

def processJsonInput(request):
    how_many_files = 0
    folder_path = "scraper_app\\files_to_insert"

    for filename in os.listdir(folder_path):
        how_many_files = how_many_files+1
        print(filename)
        path_to_file = "scraper_app\\files_to_insert\\"+filename
        with open(path_to_file, "r") as f:
            data_to_insert = json.load(f)
            exchanger = Exchanger()
            exchanger.process(data_to_insert)
    if how_many_files >0:
        return HttpResponse(str(how_many_files)+ " Files successfuly processed and stored")
    abs_path = os.path.abspath(folder_path)
    return HttpResponse("There are no files in \n"+abs_path+"\n add one in the dictionary")

def search_suggestion(request):  #for the search suggestion
    
        category = request.GET.get('category')
        searched = request.GET.get('searched')
     
        print("Category: " +category+", Searched: "+ searched)
     
        requested = []
        most_relevant = []
        if category == "smiles":
            requested = Substances.objects.filter(smiles__icontains=searched)
            for req in requested.values():
            #print(str(req))
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
            #print(str(req))
                if category in model_fields_that_are_lists:
                    for li in req["formula"]:
                        if li.startswith(searched):
                            most_relevant.append(li)
                else:
                    if req["formular"].startswith(searched):
                        most_relevant.append(req[category])
        elif category == "names":
            requested = Substances.objects.filter(names__icontains=searched)
            for req in requested.values():
                if category in model_fields_that_are_lists:
                    for li in req["names"]:
                        if li.startswith(searched):
                            most_relevant.append(li)
        elif category == "iupac_names":
            requested = Substances.objects.filter(iupac_names__icontains=searched)
            for req in requested.values():
            #print(str(req))
                if category in model_fields_that_are_lists:
                    for li in req["iupac_names"]:
                        if li.startswith(searched):
                            most_relevant.append(li)
        elif category == "cas_num":
            requested = Substances.objects.filter(cas_num__icontains=searched)
            for req in requested.values():
            #print(str(req))
                if category in model_fields_that_are_lists:
                    for li in req["cas_num"]:
                        if li.startswith(searched):
                            most_relevant.append(li)
                else:
        
                        if req[category].startswith(searched):
                            most_relevant.append(req[category])
        elif category == "categories":
            requested = Substances.objects.filter(categories__icontains=searched)
            for req in requested.values():
                    if req["categories"] not in most_relevant:
                        most_relevant.append(req["categories"])
                    
               
        elif category == "source":
            requested = Substances.objects.filter(source__icontains=searched)
            print("found"+str(len(requested)))
            for req in requested :
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

        
       # print("Response data: " +str(most_relevant))
        return JsonResponse(most_relevant, safe=False)
    