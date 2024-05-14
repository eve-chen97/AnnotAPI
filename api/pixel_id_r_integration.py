import requests
import json

def call_r_service(input_data):
    r_api_url = 'http://localhost:8000/pixel_id'
    headers = {'Content-Type': 'application/json'}

    payload = json.dumps(input_data)

    try:
        response = requests.post(r_api_url, headers=headers, data=payload)
        # Check if the request was successful
        if response.status_code == 200:
            return response.json()  # Return the JSON response
        else:
            # Handle responses other than 200
            response.raise_for_status()
    except requests.exceptions.RequestException as e:
        # Handle connection errors, timeouts, and other requests exceptions
        raise Exception("Error from R service: " + str(e)) from e

