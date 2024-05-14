import requests
import json
import subprocess
import os
import re

## Use plumber
# def call_r_service(input_data):
#     r_api_url = 'http://localhost:8000/pixel_id'
#     headers = {'Content-Type': 'application/json'}

#     payload = json.dumps(input_data)

#     try:
#         response = requests.post(r_api_url, headers=headers, data=payload)
#         # Check if the request was successful
#         if response.status_code == 200:
#             return response.json()  # Return the JSON response
#         else:
#             # Handle responses other than 200
#             response.raise_for_status()
#     except requests.exceptions.RequestException as e:
#         # Handle connection errors, timeouts, and other requests exceptions
#         raise Exception("Error from R service: " + str(e)) from e

# Use Subprocess
def call_r_service(input_data):
    # Convert the Python dictionary to a JSON string
    json_data = json.dumps(input_data)
    input_file_path = os.path.join(os.getcwd(), "input.json")
    # Save this JSON data to a file because R script expects a file path
    try:
        with open("input.json", "w") as f:
            f.write(json_data)
            print("Data written to:", input_file_path)
    except IOError as e:
        print("Error writing to file:", e)
        return None

    # Define the path to the R script
    script_path = os.path.join(os.getcwd(), 'pixel_id_r_scripts', 'pixel_id_service.R')
    if not os.path.exists(script_path):
        raise Exception(f"Script not found: {script_path}")
    
    # Debug the subprocess call
    command = ["Rscript", script_path, input_file_path]
    print("Running command:", " ".join(command))

    # Call the R script
    try:
        result = subprocess.run(
            command, capture_output=True, text=True, env=os.environ # Include environment variables
        )
        print("Raw output from Rscript, STDOUT:", result.stdout)  

        # Extracting the JSON part from the output
        json_output = re.search(r'\{.*\}', result.stdout)
        if json_output:
            output = json.loads(json_output.group())
            return output
        else:
            raise ValueError("No valid JSON found in R output")

        # print("STDERR:", result.stderr)
        # result.check_returncode()  # This will raise an exception if the call failed
        # # Parse the output from R
        # output = json.loads(result.stdout)
        # return output
    except subprocess.CalledProcessError as e:
        print("Error in R script:", e.stderr)
        raise Exception("R service failed failed from subprocess.")
    
    except json.JSONDecodeError as e:
        print("JSON decoding failed:", e)
        raise

# Example use
if __name__ == "__main__":
    payload = {
    "input_new_file": "test_population_data.csv", 
    "polygons_new_file": "new_input_polygons_all.json", 
    "alignment_data": "alignment_results.csv", 
    "input_master_map": "master_map_5_CD4_CD8.csv", 
    "polygons_master_map": "master_map_polygons_5_CD4_CD8.json"
}
    result = call_r_service(payload)
    print("R returned:", result)