import os
import sys
import subprocess
import time
from datetime import datetime
import webbrowser

# Pfad des Webscrapers
path_to_webscraper = "C:/Users/User/Desktop/Webscraper_PIHKAL"  # Ändere dies für Windows
# path_to_webscraper = "/home/User/Desktop/Webscraper_PIHKAL"  # Ändere dies für Linux

# Logging-Verzeichnis
logging_directory = "logging_directory"
path_to_logging_directory = os.path.join(path_to_webscraper, logging_directory)

# Generiere einen Zeitstempel für die Log-Datei
timestamp = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
logfile = os.path.join(path_to_logging_directory, f"log_{timestamp}.txt")

# Ändere das Verzeichnis
os.chdir(path_to_webscraper)

# Installiere benötigte Python-Pakete
def install_packages():
    packages = ["aiohttp", "asyncio", "beautifulsoup4", "django", "rdkit", "periodictable"]
    for package in packages:
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])

install_packages()

# Umleitung der print-Ausgaben in die Log-Datei
class LogFileWriter:
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.log = open(filename, "a")

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def flush(self):
        self.terminal.flush()
        self.log.flush()

sys.stdout = LogFileWriter(logfile)

print("Starting Django development server...")

# Startet den Django-Server und leitet die Ausgabe in die Log-Datei um
server_process = subprocess.Popen(
    ["python", "manage.py", "runserver"],
    stdout=sys.stdout.log,
    stderr=subprocess.STDOUT
)

print("Waiting for the server to start...")
time.sleep(5)

print("Opening browser...")
# Browser-URL
url = "http://127.0.0.1:8000/webscraper"

# Starten des Browsers
webbrowser.open(url)

# Beispielweise: Hier könntest du eine bestimmte Zeit warten oder eine andere Methode verwenden, um den Browser-Status zu überprüfen

try:
    # Warten auf eine Benutzereingabe oder andere Bedingungen, um den Server zu beenden
    input("Press Enter to stop the server...")
finally:
    print("Stopping Django server...")
    # Beende den Django-Server-Prozess
    server_process.terminate()  # Oder server_process.kill() falls terminate() nicht funktioniert
    server_process.wait()       # Warten bis der Prozess beendet ist

print("Django server stopped.")
