import json

with open('./out_files/data_2020.json', 'r') as json_file:
   events_data = json.load(json_file)
with open('./out_files/data_2021.json', 'r') as json_file:
   events_data = events_data | json.load(json_file)
with open('./out_files/data_2022.json', 'r') as json_file:
   events_data = events_data | json.load(json_file)

print(events_data.items())

from flask import Flask, jsonify, request
from flask_cors import CORS
from datetime import datetime

app = Flask(__name__)
CORS(app)  # Allow all origins to access this API



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

