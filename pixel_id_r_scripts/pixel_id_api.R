# Turning the Pixel ID R functions into a RESTful API using plumber.
library(plumber)
library(jsonlite)

# Load Pixel ID R scripts
# source('label_to_alphahull.R')
source('pixel_id_r_scripts/pixel_ID_func.R')
source('pixel_id_r_scripts/pixel_ID_WebService.R')

#* @post /pixel_id
function(req){
  # Parse JSON body
  data <- fromJSON(req$body)
  
#   # Defaults for the optional parameters
#   Do_Plots <- FALSE # Dont't generate plots
#   matrix_division <- 32

  # Here you run your pixel_ID_WebService function and other operations needed
  pixel_ID_scores <- pixel_ID_WebService(input_new_file = data$input_new_file,  # original population data (csv)
                                polygons_new_file = data$polygons_new_file,   # polygon data (json)
                                alignment_data = data$alignment_data,  # alighment data (csv)
                                input_master_map = data$input_master_map, # test_data master map (csv)
                                polygons_master_map = data$polygons_master_map, # test_data master map polygon (json)
                                Do_Plots = data$Do_Plots, 
                                matrix_division = data$matrix_division)
  
  # Return the result as JSON
  return(toJSON(pixel_ID_scores))
}
# 

# Run the API
plumber::plumb("pixel_id_api.R")$run(port=8000)
