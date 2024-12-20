from flask import Flask, jsonify, request
from flask_cors import CORS
from datetime import datetime
import json

app = Flask(__name__)
CORS(app)  # Allow all origins to access this API

with open('./out_files/data.json', 'r') as json_file:
   events_data = json.load(json_file)
            
#events_data = {
#    "2024-12-01": {
#        "launches": [
#            {"id": "2020-001", "site": "Site1", "rocket": "Falcon 9", "lat": 34.0522, "lon": -118.2437, "emissions": {"CO2": 100, "BC": 10}, "smc": "True"},
#            {"id": "2020-002", "site": "Site2", "rocket": "Delta IV", "lat": 51.5074, "lon": -0.1278, "emissions": {"CO2": 120, "BC": 12}, "smc": "False"}
#        ],
#        "reentries": [
#            {"id": "2020-001", "site": "Site1", "reusability": "Reusable", "lat": 24.0522, "lon": 30.2437,"emissions": {"Al2O3": 100, "NOx": 10, "Mass": 40},
#             "smc": "True", "name": "Starlink 0001", "category": "Payload"},
#             {"id": "2020-003", "site": "Site2", "reusability": "Discarded (known)", "lat": -50.0522, "lon": 60.2437,"emissions": {"Al2O3": 40, "NOx": 20, "Mass": 60},
#             "smc": "False", "name": "Soyuz S1", "category": "S1"},
#        ]
#    },
#    "2024-12-02": {
#        "launches": [
#            {"id": "2020-003", "site": "Site3", "rocket": "Soyuz", "lat": 64.0522, "lon": -15.2437, "emissions": {"CO2": 100, "BC": 10}, "smc": "False"},
#            {"id": "2020-004", "site": "Site4", "rocket": "Electron", "lat": 11.5074, "lon": 30.1278, "emissions": {"CO2": 120, "BC": 12}, "smc": "True"}
#        ],
#        "reentries": []
#    },
#    "2024-12-03": {
#        "launches": [
#            {"id": "2020-005", "site": "Site5", "rocket": "Kuaizhou-1", "lat": 34.0522, "lon": -40.2437, "emissions": {"CO2": 100, "BC": 10}, "smc": "False"},
#        ],
#        "reentries": [
#            {"id": "2020-007", "site": "Site3", "reusability": "Discarded (approx.)", "lat": 10.0522, "lon": 20.2437,"emissions": {"Al2O3": 50, "NOx": 80, "Mass": 10},
#             "smc": "True", "name": "Starlink Deployment Rail", "category": "Component"},
#        ]
#    }
#}

@app.route('/api/launches', methods=['GET'])
def get_launches():
    # Get optional query parameters for filtering by date range
    start_date = request.args.get('start_date')
    end_date = request.args.get('end_date')

    try:
        if start_date:
            start_date = datetime.strptime(start_date, '%Y-%m-%d')
        if end_date:
            end_date = datetime.strptime(end_date, '%Y-%m-%d')
    except ValueError:
        return jsonify({"error": "Invalid date format. Please use YYYY-MM-DD."}), 400

    # Filter the launches data based on the date range
    filtered_data = []
    for date_str, events in events_data.items():
        event_date = datetime.strptime(date_str, '%Y-%m-%d')
        if (not start_date or event_date >= start_date) and (not end_date or event_date <= end_date):
            filtered_data.append({"date": date_str, "launches": events["launches"]})

    return jsonify(filtered_data)

@app.route('/api/reentries', methods=['GET'])
def get_reentries():
    # Get optional query parameters for filtering by date range
    start_date = request.args.get('start_date')
    end_date = request.args.get('end_date')

    try:
        if start_date:
            start_date = datetime.strptime(start_date, '%Y-%m-%d')
        if end_date:
            end_date = datetime.strptime(end_date, '%Y-%m-%d')
    except ValueError:
        return jsonify({"error": "Invalid date format. Please use YYYY-MM-DD."}), 400

    # Filter the reentries data based on the date range
    filtered_data = []
    for date_str, events in events_data.items():
        event_date = datetime.strptime(date_str, '%Y-%m-%d')
        if (not start_date or event_date >= start_date) and (not end_date or event_date <= end_date):
            filtered_data.append({"date": date_str, "reentries": events["reentries"]})

    return jsonify(filtered_data)

if __name__ == "__main__":
    app.run()

