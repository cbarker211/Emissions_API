import json
from flask import Flask, jsonify, request
from flask_cors import CORS
from datetime import datetime
import os

app = Flask(__name__)
CORS(app)  # Allow all origins to access this API

data_dir = './out_files'
events_data = []

# Load all yearly files from 1957–2025 if they exist
for year in range(1957, 2026):
    file_path = os.path.join(data_dir, f'{(year // 10) * 10}/data_{year}.json')
    if os.path.exists(file_path):
        with open(file_path, 'r') as json_file:
            yearly_data = json.load(json_file)

        for date_str, events in yearly_data.items():
            events_data.append(
                (datetime.strptime(date_str, "%Y-%m-%d"),
                date_str,
                events))

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
    for event_date, date_str, events in events_data:
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
    for event_date, date_str, events in events_data:
        if (not start_date or event_date >= start_date) and (not end_date or event_date <= end_date):
            filtered_data.append({"date": date_str, "reentries": events["reentries"]})

    return jsonify(filtered_data)

@app.route('/api/test')
def test():
    return {"hello": "world"}

if __name__ == "__main__":
    app.run()

