import os
from pixel_id_r_integration import call_r_service

directory_path = os.path.join(os.getcwd(), 'test_data')

input_new_file = 'test_population_data.csv'
input_file_path = os.path.join(directory_path, input_new_file)
polygons_new_file = 'new_input_polygons_all.json'
polygons_new_path = os.path.join(directory_path, polygons_new_file)
alignment_data = 'alignment_data.csv'
alignment_path = os.path.join(directory_path, alignment_data)
master_map_file = 'master_map_5_CD4_CD8.csv'
mastermap_path = os.path.join(directory_path, master_map_file)
master_polygon_file = 'master_map_polygons_5_CD4_CD8.json'
master_polygon_path = os.path.join(directory_path, master_polygon_file)

if not os.path.exists(polygons_new_path):
    raise Exception(f"File not found: {polygons_new_path}")

# Define whether plots should be generated and the division matrix
do_plots = False
matrix_division = 32

payload = {
    "input_new_file": input_file_path,
    "polygons_new_file": polygons_new_path,
    "alignment_data": alignment_path,
    "input_master_map": mastermap_path,
    "polygons_master_map": master_polygon_path,
    "Do_Plots": do_plots,
    "matrix_division": matrix_division
}

try:
    result = call_r_service(payload)
    print("R returned:", result)
except Exception as ex:
    print(ex)