from flask import Flask, request, jsonify, send_file
import pandas as pd
import os
import tempfile
from werkzeug.utils import secure_filename
#from population_data import PopulationData
#from polygon_data import PolygonData
from Alignment_WebService_func import alignment_2D_flowdata  # Import alignment function
from pixel_id_r_integration import call_r_service

app = Flask(__name__)

# Configure the maximum upload size to 16MB (adjust if necessary)
app.config['MAX_CONTENT_LENGTH'] = 16 * 1024 * 1024

# Create data storage

# @app.route('/align/polygon', methods=['POST'])
# def process_polygon_data():
#     # This assumes the incoming request has a JSON payload with polygon vertex data.
#     json_data = request.get_json()
#     polygon_data = PolygonData.from_json(json_data)
    
#     # Do something with the polygon data, such as passing it to the alignment script
#     # For example:
#     # aligned_data = some_alignment_function(polygon_data.get_all_vertices())
    
#     # Return some result or response
#     return jsonify({'status': 'success', 'data': polygon_data.get_all_vertices()})

# @app.route('/upload', methods=['POST'])
# def upload_population_data():
#     if 'file' not in request.files:
#         return jsonify({"error": "No file part"}), 400
#     file = request.files['file']
    
#     # Read the file to string and parse it as JSON
#     # Assuming the file is a CSV file
#     csv_content = file.stream.read().decode('utf-8')
#     reader = csv.DictReader(csv_content.splitlines())
    
#     # Convert the CSV to a list of dictionaries
#     data_dicts = [row for row in reader]
    
#     # Create an instance of PopulationData with the loaded data
#     population_data = PopulationData(data_dicts)
    
#     # Do something with the population data, like passing it to your alignment function
#     # ...

#     # Return a response or some results
#     return jsonify({"status": "success", "data": population_data.get_marker_data('CD8')})

# @app.route('/align', methods=['POST'])
# def align():
#     # Check if the post request has the file part
#     if 'file' not in request.files:
#         return jsonify({"error": "No file part"}), 400
#     file = request.files['file']

#     # If the user does not select a file, the browser submits an
#     # empty file without a filename.
#     if file.filename == '':
#         return jsonify({"error": "No selected file"}), 400

#     if file:
#         # Save the uploaded file to a temporary file
#         temp_upload_dir = tempfile.mkdtemp()
#         file1_path = os.path.join(temp_upload_dir, file.filename)
#         file.save(file1_path)

#         # Set parameters
#         project = "ANNOTAPI"
#         Fast = True
#         bin_size = 64
#         verbose = True
#         save = False
#         sample_size = 50000

#         # if(save):
#         #     temp_results_dir = tempfile.mkdtemp()
#         #     plots_dir = os.path.join(temp_results_dir, 'plots_pairwise_comparison')
#         #     shift_dir = os.path.join(temp_results_dir, 'pairwise_shift')
#         #     os.makedirs(plots_dir, exist_ok=True)
#         #     os.makedirs(shift_dir, exist_ok=True)

#         # Set the path to your master map file
#         directory_path = os.path.join(os.getcwd(), 'test_data')
#         master_map_file = 'master_map_5_CD4_CD8.csv'
#         file2_path = os.path.join(directory_path, master_map_file)

#         # Run your alignment function
#         results_df = alignment_2D_flowdata(
#             file1=file1_path,
#             file2=file2_path,
#             j=1,
#             k=0,
#             project=project,
#             Fast=Fast,
#             bin_size=bin_size,
#             verbose=verbose,
#             save=save,
#             sample_size=sample_size
#         )

#         # Convert DataFrame to CSV and send it back as a file response
#         result_csv = results_df.to_csv(index=False)

#         # Clean up the temporary directories after use
#         os.remove(file1_path)
#         os.rmdir(temp_upload_dir)
#         # Potentially clean up the temp_results_dir if needed

#         return result_csv


@app.route('/pixel_id', methods=['POST'])
def pixel_id():
    if 'input_new_file' not in request.files or 'polygon_file' not in request.files:
        return jsonify({"error": "Missing file parts (expected 'input_new_file' and 'polygon_file')"}), 400

    input_new_file = request.files['input_new_file']
    polygon_file = request.files['polygon_file']

    # If the user does not select a file, the browser submits an empty file without a filename.
    if input_new_file.filename == '' or polygon_file.filename == '':
        return jsonify({"error": "Missing selected file (expected 'input_new_file' and 'polygon_file')"}), 400

    temp_upload_dir = tempfile.mkdtemp()
    input_file_path = os.path.join(temp_upload_dir, secure_filename(input_new_file.filename))
    polygon_path = os.path.join(temp_upload_dir, secure_filename(polygon_file.filename))
       
    input_new_file.save(input_file_path)
    polygon_file.save(polygon_path)

    # Files of master map
    directory_path = os.path.join(os.getcwd(), 'test_data')
    master_map_file = 'master_map_5_CD4_CD8.csv'
    mastermap_path = os.path.join(directory_path, master_map_file)
    master_polygon_file = 'master_map_polygons_5_CD4_CD8.json'
    master_polygon_path = os.path.join(directory_path, master_polygon_file)

    # Run alignment function first       
    # Get alignment results
    alignment_results_df = alignment_2D_flowdata(
        file1=input_file_path,
        file2=mastermap_path,
        j=1,
        k=0,
        project="ANNOTAPI",
        Fast=True,
        bin_size=64,
        verbose=False,
        save=False,
        sample_size=50000
    )

    # Save the alignment results and get directory
    alignment_path = os.path.join(temp_upload_dir, 'alignment_results.csv')
    # Write the CSV data to a file using DataFrame's to_csv method
    # alignment_result_csv = alignment_results_df.to_csv(index=False)
    alignment_results_df.to_csv(alignment_path, index=False)

    input_data = {
        "input_new_file": input_file_path,
        "polygons_new_file": polygon_path,
        "alignment_data": alignment_path,
        "input_master_map": mastermap_path,
        "polygons_master_map": master_polygon_path
    }

    # Try to call the R service
    try:
        results = call_r_service(input_data)
        return jsonify(results)
    except Exception as e:
        return jsonify({"error": str(e)}), 500
    finally:
        # Clean up files
        os.remove(input_file_path)
        os.remove(polygon_path)
        os.remove(alignment_path)
        os.rmdir(temp_upload_dir)

if __name__ == '__main__':
    app.run(debug=True)