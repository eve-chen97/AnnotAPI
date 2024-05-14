library(jsonlite)
# Source Pixel ID R functions
source('pixel_id_r_scripts/functions/pixel_ID_func.R')
source('pixel_id_r_scripts/functions/pixel_ID_WebService.R')
source('pixel_id_r_scripts/functions/label_to_alphahull.R')

# Read data from a JSON file passed as a command-line argument
args <- commandArgs(trailingOnly = TRUE)
print(args[1])
if (!file.exists(args[1])) {
  stop("File does not exist: ", args[1])
}
json_input <- fromJSON(args[1])

# Extract data from JSON
input_new_file <- json_input$input_new_file
polygons_new_file <- json_input$polygons_new_file
alignment_data <- json_input$alignment_data
input_master_map <- json_input$input_master_map
polygons_master_map <- json_input$polygons_master_map
Do_Plots <- json_input$Do_Plots
matrix_division <- json_input$matrix_division

# Call the R function
result <- tryCatch({
  pixel_ID_WebService(
    input_new_file = input_new_file,
    polygons_new_file = polygons_new_file,
    alignment_data = alignment_data,
    input_master_map = input_master_map,
    polygons_master_map = polygons_master_map,
    Do_Plots = Do_Plots,
    matrix_division = matrix_division
  )
}, error = function(e) {
  list(error = "Error processing data", message = e$message)
})

# Output the result as JSON
cat(toJSON(result))