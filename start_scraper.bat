@echo off

REM Pfad des Webscrapers
set path_to_webscraper=C:\Users\User\Desktop\Webscraper_PIHKAL


set original_string=Webscraper_PIHKAL
set logging_directory=logging_directory
set path_to_logging_directory=%path_to_webscraper%\%logging_directory%

cd %path_to_webscraper%

REM Installiere benötigte Python-Pakete
pip install aiohttp
pip install asyncio
pip install bs4
pip install django
pip install rdkit
pip install periodictable
pip install json

REM Generiere einen Zeitstempel für die Log-Datei
set timestamp=%date%_%time%
set timestamp=%timestamp::=-%
set timestamp=%timestamp:/=-%
set timestamp=%timestamp: =_%
set timestamp=%timestamp:.=-%
set hour_minute=%timestamp:~0,16%


REM Erstellen der Log-Datei mit Zeitstempel
set logfile=%path_to_logging_directory%\logging_%hour_minute%.txt

echo Starting Django development server...

REM Startet den Django-Server in einem neuen Fenster und leitet die Ausgabe in die Log-Datei um
start cmd /k "python manage.py runserver > %logfile% 2>&1"
echo Waiting for the server to start...
ping 127.0.0.1 -n 5 > nul


echo Opening browser...
start http://127.0.0.1:8000/webscraper


