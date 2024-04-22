import requests
import json
import os

r_api_url = 'http://localhost:8000/pixel_id'
headers = {'Content-Type': 'application/json'}

directory_path = os.path.join(os.getcwd(), 'test_data')

input_new_file = 'test_population_data.csv'
input_file_path = os.path.join(directory_path, input_new_file)
polygons_new_file = 'new_input_polygons.json'
polygons_new_path = os.path.join(directory_path, polygons_new_file)
alignment_data = 'alignment_data.csv'
alignment_path = os.path.join(directory_path, alignment_data)
master_map_file = 'master_map_5_CD4_CD8.csv'
mastermap_path = os.path.join(directory_path, master_map_file)
master_polygon_file = 'master_map_polygons_5_CD4_CD8.json'
master_polygon_path = os.path.join(directory_path, master_polygon_file)


payload = {
    "input_new_file": input_file_path,
    "polygons_new_file": polygons_new_path,
    "alignment_data": alignment_path,
    "input_master_map": mastermap_path,
    "polygons_master_map": master_polygon_path,
}

response = requests.post(r_api_url, headers=headers, data=json.dumps(payload))
print(response.json())