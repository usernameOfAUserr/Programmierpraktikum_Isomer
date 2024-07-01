USAGE:
Download the Repository

Under 
-Windows: Use start_scraper.bat
-Linux: Use start_scraper.py
Whatever you use, change the path to the directory, and then you can drag the script to wherever you want it to be stored. After that, you can use it to start the web scraper and use it as you please.

INSERT DATA / IMPORT:
If you want to insert data into the database, please make sure that it is in the correct format (schema.json).

NOTE: If you are done with inserting, remove the JSON file. If you don't, the next insertion will have an unnecessarily long duration.

GENERATE DATA / EXPORT:
If you want to generate a JSON file with all the currently stored substances, click "Generate JSON" in the GUI. The file will be generated with the name PIHKAL.json.

FUNCTIONS IN GENERAL:
-Delete
-Scrape: https://isomerdesign.com/pihkal (just new substances or reset the whole DB)
Note: When you reset the DB, that DOES NOT affect stored data from other websites.
-Search
-Insert substances
-Generate

STRUCTURE:
OOP: The code follows Object-Oriented Programming principles.
DB: Uses Django models in models.py. You can change it by writing whatever you want and then run python manage.py makemigrations and then python manage.py migrate.
Note: When there are problems with changing the DB structure, generate a JSON file with all substances of the DB, then delete db.sqlite3 and create a new DB with python manage.py migrate, and then insert the substances.
Coordination and Management: The functions are managed in views.py.
Scraper Application: The names of the .py files in scraper_app are self-explanatory.

LOGGING:
the logging files are in logging_directory