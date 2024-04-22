import requests
import json

def call_r_service(input_data):
    r_api_url = 'http://localhost:8000/pixel_id'
    headers = {'Content-Type': 'application/json'}
    payload = json.dumps(input_data)

    response = requests.post(r_api_url, headers=headers, data=payload)
    if response.status_code == 200:
        return response.json()  # Return the JSON response
    else:
        raise Exception("Error from R service: " + response.text)
