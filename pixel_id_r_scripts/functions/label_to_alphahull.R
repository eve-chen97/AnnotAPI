# Function to get polygons from population input
# Currently not used

########################################
# label_to_polygons_alphahull
########################################
label_to_polygons_alphahull <- function (df){
    colnames(df) <- c("x","y","classes")
    # Create a sample data frame of points with labels
    jitter_points <- jitter(as.matrix(df[, c("x", "y")]), amount = 1e-9)
    # Create a list to store alpha shapes for each label
    alpha_shapes <- list()
    vertices_list <- list()
    
    # Loop through unique labels
    for (label in unique(df$classes)) { # changed from df$label
        if(label == 0){next}
        # Extract points for the current label
        points_label <- jitter_points[df$classes == label, ]
        
        # Generate alpha shape for the current label with alpha = 2
        alpha_shape <- tryCatch(alphahull::ashape(points_label, alpha = 2), error = function(e) NULL)
        if(is.null(alpha_shape)){
            print("something wrong in ashape, skipping")
            return(NULL)
        }
        vertices <- tryCatch(ahull2polygon(n10=alpha_shape),  error = function(e) NULL)
        if(is.null(vertices)){
            print("something wrong in ahull2polygon, skipping")
            return(NULL)
        }else if(length(vertices) == 1 && vertices == 1){
            print("something wrong in ahull2polygon, skipping")
            return(NULL)
        }
        # Store alpha shape in the list
        alpha_shapes[[label]] <- alpha_shape
        vertices_list[[label]] <- vertices
    }
    
    return(vertices_list)
  
}

########################################
# ahull2polygon
########################################
ahull2polygon <- function(n10){
    # The return structure of ashape doesn't make it easy to get the polygon out
    # We can re-order the points using the igraph package. 
    # The ashape function returns an object that has the indexes of the start and end point of each segment. We can feed these in as an edge list for a graph object:
    # Input : n10 is an alphahull
    # Output : 
    n10g = graph.edgelist(cbind(as.character(n10$edges[, 'ind1']), as.character(n10$edges[, 'ind2'])), directed = FALSE)
    if(!is.connected(n10g)){
        return(1)
        # stop("Graph not connected")
    }
    # if(any(degree(n10g) != 2)){
    #   stop("Graph not circular")
    # }
    if(igraph::clusters(n10g)$no > 1){
        stop("Graph composed of more than one circle")
    }

    cutg = n10g - E(n10g)[1]
    # find chain end points
    ends = names(which(degree(cutg) == 1))
    path = get.shortest.paths(cutg, ends[1], ends[2])[[1]][[1]]
    # this is an index into the points
    pathX = as.numeric(V(n10g)[path]$name)
    # join the ends
    pathX = c(pathX, pathX[1])
    # now show the alpha shape plot with our poly on top
    #plot(n10, lwd = 10, col = "gray")
    # get the points from the ashape object
    #lines(n10$x[pathX, ], lwd = 2)

    return(poly = n10$x[pathX, ])
}

########################################
# label_to_polygons_convexhull
########################################
# function to get the convex hull of all gates
label_to_polygons_convexhull<-function(gated_df, concavity_val){
    colnames(gated_df) <- c("x","y","classes")
    all_classes <- unique(gated_df$classes)
    list_df_hull <- list()
    for(classes in all_classes){
        if(classes == 0){next}
        inds <- which(gated_df$classes==classes)
        df_current_classes <- gated_df[inds,]
        if(nrow(df_current_classes) > 200000){ # if there are many events I use the concave function
            df_current_classes_hull_values <- as.data.frame(concaveman(as.matrix(df_current_classes[,c(1,2)]),concavity=5))
        } else { # otherwise I use the classical convex hull
            current_hull_indices <- chull(df_current_classes[,c(1,2)]) # indices points at the border of the convex hull for the current classes
            df_current_classes_hull_values <- df_current_classes[current_hull_indices,c(1,2)]
            df_current_classes_hull_values <- rbind(df_current_classes_hull_values, df_current_classes_hull_values[1,])
        }
        vec_group <- rep(sprintf("%s",classes),nrow(df_current_classes_hull_values))
        df_current_hull <- cbind(df_current_classes_hull_values,vec_group)
        colnames(df_current_hull) <- c("x","y","group_gate")
        list_df_hull[[sprintf("%s",classes)]] <- df_current_hull[,c("x", "y")]  
    }
    return(list_df_hull)
}