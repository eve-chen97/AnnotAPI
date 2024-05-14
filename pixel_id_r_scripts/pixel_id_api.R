# Turning the Pixel ID R functions into a RESTful API using plumber.
library(plumber)
library(jsonlite)
#library(alphahull) 

# Load Pixel ID R scripts
# source('label_to_alphahull.R')
source('pixel_id_r_scripts/functions/pixel_ID_func.R')
source('pixel_id_r_scripts/functions/pixel_ID_WebService.R')

#* @post /pixel_id
function(req, res) {
    data <- tryCatch({
        jsonlite::fromJSON(req$postBody)
    }, error = function(e) {
        res$status <- 400
        return(list(error = "Invalid JSON input", message = e$message))
    })
    
    # Check if data is correctly parsed before proceeding
    if (!is.list(data)) {
        res$status <- 400
        return(list(error = "Processed data is not a list. JSON input may be incorrect."))
    }

    tryCatch({
        # Call R function using the data extracted from the request
        result <- pixel_ID_WebService(
            input_new_file = data$input_new_file,
            polygons_new_file = data$polygons_new_file,
            alignment_data = data$alignment_data,
            input_master_map = data$input_master_map,
            polygons_master_map = data$polygons_master_map,
            Do_Plots = data$Do_Plots %||% FALSE,  # Defaulting Do_Plots if not provided
            matrix_division = data$matrix_division %||% 32  # Default matrix division
        )
        res$status <- 200
        return(list(result = result))
    }, error = function(e) {
        res$status <- 500
        return(list(error = "Error processing data", message = e$message))
    })
}


