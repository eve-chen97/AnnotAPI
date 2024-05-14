library("sp")
library("stringr")
library("plot3D")
library("jsonlite")
# if (!require("scattermore")) install.packages("scattermore")
#     library("scattermore")
# library("alphahull") # 

pixel_ID_WebService <- function(input_new_file_path, # user upload
                                polygons_new_file_path, # user upload
                                alignment_data_path, # pass from the alignment function
                                input_master_map_path, 
                                polygons_master_map_path, # this is not mandatory 
                                Do_Plots = F, 
                                matrix_division = 32){

    if (Do_Plots){
        png(paste0("results/Annotation/master_map_original.png"),width = 600, height = 600+0)
        x_range <- c(0,1)
        y_range <- c(0,1)
        plot_dens_wrapper(input_master_map[,1:2], xlim=x_range, ylim=y_range)
        points(input_master_map[,1:2], col=input_master_map[,3], pch=19, cex=0.2)
        for (j in 1:length(polygons_master_map)){ lines(polygons_master_map[[j]], col= j)}
        dev.off()
    }
    
    # polygons_new_file = label_to_polygons_alphahull(df=input_new_file)
    # if(length(polygons_new_file) == 0){
    #     print("alpha hull did not work")
    #     polygons_new_file = label_to_polygons_convexhull(gated_df = input_new_file, concavity_val = 5)
    # }
    # print(paste("length of polygon is: ", length(polygons_new_file), sep = ''))
    # polygons_new_file <- polygons_new_file[order(names(polygons_new_file))] # reorder to match up the colours for scatter plot
    
    if (Do_Plots){
        png(paste0("results/Annotation/input_new_file_original.png"),width = 600, height = 600+0)
        x_range <- unname(unlist(c(alignment_data["x_min"], x_max=alignment_data["x_max"])))
        y_range <- unname(unlist(c(alignment_data["y_min"], x_max=alignment_data["y_max"])))
        plot_dens_wrapper(input_new_file[,1:2], xlim=x_range, ylim=y_range)
        points(input_new_file[,1:2], col=input_new_file[,3], pch=19, cex=0.2)
        for (j in 1:length(polygons_new_file)){ lines(polygons_new_file[[j]], col= j)}
        dev.off()
    }
    
    #some extra steps for tests
    input_new_file <- read.csv(input_new_file_path)
    polygons_new_file <- fromJSON(polygons_new_file_path)
    #polygons_new_file <- lapply(polygons_new_file_path, function(x) as.data.frame(fromJSON(x))) %>% dplyr::bind_rows()
    alignment_data <- read.csv(alignment_data_path)
    input_master_map <- read.csv(input_master_map_path)
    polygons_master_map <- fromJSON(polygons_master_map_path)
    
    CIDs <- calc_pixel_ID(data=input_new_file, 
                          polygons=polygons_new_file, # user upload polygon
                          matrix_division=matrix_division, 
                          Do_Plots=Do_Plots, 
                          plot_name="new_file",
                          x_min=alignment_data["x_min"], x_max=alignment_data["x_max"], 
                          y_min=alignment_data["y_min"], y_max=alignment_data["y_max"])
    
    CIDs_master <- calc_pixel_ID(data=input_master_map, 
                                 polygons=polygons_master_map, # master map polygons
                                 matrix_division=matrix_division, 
                                 Do_Plots=Do_Plots, 
                                 plot_name="master_map",
                                 x_min=0, x_max=1, 
                                 y_min=0, y_max=1)
    
    
    
    s1_len <- length(CIDs)
    s2_len <- length(CIDs_master)
    scores <- scores_old <- matrix(0, nrow=s1_len, ncol=s2_len)
    for (index1 in 1:s1_len){
        for (index2 in 1:s2_len){
            scores[index1,index2] <- string_dist(CIDs[index1], CIDs_master[index2])
        }
    }
    
    matches <- sapply(1:nrow(scores), function(x){
        which(scores[x,] == max(scores[x,]))[1]
    })
    
    pixel_ID_scores <- scores[cbind(1:s1_len, matches)]
    pixel_ID_names <- names(polygons_master_map)[matches]
    #names(pixel_ID_scores) <- pixel_ID_names
    
    #return(pixel_ID_scores)
    
    # Convert the results to a JSON object
    json_output <- toJSON(list(names = pixel_ID_names, scores = pixel_ID_scores))
    cat(json_output)
}


