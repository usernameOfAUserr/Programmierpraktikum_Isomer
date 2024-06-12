from datetime import datetime, timezone

# Aktuelles Datum und Uhrzeit erhalten
current_datetime = datetime.now(timezone.utc)

# Formatieren des datetime-Objekts in das gewünschte Format
formatted_string = current_datetime.strftime('%Y-%m-%dT%H:%M:%SZ')

print(formatted_string)