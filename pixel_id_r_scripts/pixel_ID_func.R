##########################################################
# polygons_2_pixels
##########################################################
polygons_2_pixels <- function(user_polygons, x.range.temp=NULL, y.range.temp=NULL, x.range = NULL, y.range = NULL, probability=NULL, source, source_type, num_polygons = NULL, final.cluster = NULL, matrix_division = NULL, just_event_label = F, label_correction_mode = "none centroid"){

    if (length(user_polygons) == 0){
        warning("Please enter a valid user polygon list")
    }
    
    if (is.null(matrix_division)){
        matrix_division <- 5
    }
    # get the bivariate data plots in
    if (source_type == 'directory'){
        data = read.csv(source, header = TRUE, stringsAsFactors = F)
    } else if (source_type == 'csv'){
        data = source
    } else {
        warning('Please enter a valid source')
        break
    }
    colnames(data) = c("x", "y")
    
    # generate okay.cluster
    okay.cluster = c()
    if (!is.null(probability)){
        first = which(probability == max(probability))[1]
        okay.cluster = c(first, first + 1, first - 1, first + 2, first + 3, first + 4, first + 5, first + 6, first + 7)
        if (probability[first] < 75){
            second = which(probability == sort(probability, decreasing = T)[2])
            if (length(second)>1){second = second[2]}
            okay.cluster = append(okay.cluster, c(second, second-1, second+1, second+2, second+3, second+4, second+5, second+6, second+7))
        }
        okay.cluster = unique(okay.cluster)
        discard.ind = which(okay.cluster < 1 | okay.cluster > 8)
        if (length(discard.ind)!=0){okay.cluster = okay.cluster[-discard.ind]}
      
    } else {
        if (!just_event_label){
            warning('Please enter a valid probability vector')
            break
        }
    }
    
    if (is.null(x.range.temp)){x.range.temp = range(data$x)}
    if (is.null(y.range.temp)){y.range.temp = range(data$y)}
    
    # if (is.null(x.range)){x.range = abs(diff(range(data$x)))}
    # if (is.null(y.range)){y.range = abs(diff(range(data$y)))}
  
    #################################################################################
    ## getting reading for the following tests and events label
    #################################################################################
    
    #### setting up everything we need later
    # obtain x and y range of each polygon, a spatialpolygons object and a spatialpoint object for further analysis
    polygon.ranges = data.frame(x = rep(-99, length(user_polygons)), y = rep(-99, length(user_polygons)))
    ps.list = list()
    
    # obtain straight lines parallel to both axis
    straight.x = c()
    straight.y = c()
    
    #### this is for event_labels
    matrix = matrix_generation(division = matrix_division)
    
    if (is.null(matrix)){
        matrix <- read.csv("~/24l/FlowGrid24l/correction_matrix.csv", header = TRUE, stringsAsFactors = F)
    }
    grid.x.range = x.range.temp
    grid.y.range = y.range.temp
    grid = data.frame(x = grid.x.range[1] + matrix$unit_x*abs(diff(grid.x.range)), y = grid.y.range[1] + matrix$unit_y*abs(diff(grid.y.range)))
    grid = SpatialPoints(grid)
    event_labels = rep(0, length(data$x))
    spatial.pts = SpatialPoints(data)
    label_stand = data.frame(score = rep(0, length(user_polygons)), grid = rep(0, length(user_polygons)), cx = rep(-99, length(user_polygons)), cy = rep(-99, length(user_polygons)))
  
    grid.table.total <- list()
    #### begin iteration over the polygons
    for (i in 1:length(user_polygons)){
        id = paste('gamer', i)
        # get the Polygons object
        vertices = user_polygons[[i]]
        p = Polygon(vertices)
        ps = Polygons(list(p), ID = id)
        # add it in the list
        ps.list = append(ps.list, ps)
        # fill in the ranges
        polygon.ranges$x[i] = abs(diff(range(vertices[,1])))
        polygon.ranges$y[i] = abs(diff(range(vertices[,2])))
        # is there a vertical line
        occur.x = table(vertices[-length(vertices[,1]),1])
        unique.x = as.numeric(names(occur.x))
        x.index = which(occur.x>1)
        if (length(x.index) != 0){
            straight.x = append(straight.x, unique.x[x.index])
        }
        # is there a horizontal line
        occur.y = table(vertices[-length(vertices[,2]),2])
        unique.y = as.numeric(names(occur.y))
        y.index = which(occur.y>1)
        if (length(y.index) != 0){
            straight.y = append(straight.y, unique.y[y.index])
        }
    
        #### this part is for event labels
        label_stand[i,c(3,4)] = get.gc(p)
        grid.in = array(over(grid, SpatialPolygons(list(ps))))
        grid.table = table(matrix$score[which(!is.na(grid.in))])
        grid.table.total[[i]] <- names(grid.table[which(grid.table>=5)])
    }
    pixels <- grid.table.total
    return (list(pixels=pixels)) 
}

##########################################################
# matrix generation 
##########################################################

matrix_generation <- function (division, output_folder = "~/24l/FlowGrid24l/"){
    start = Sys.time()
    matrix = data.frame(unit_x = rep(-99, (5*division)^2), unit_y = rep(-99, (5*division)^2), score = rep(0, (5*division)^2))
    score.vector = c()
    x.vector = c()
    y.vector = c()
    row.score = c()
    quad.count = 5 * division
    row.length = 5*division
    unit.loc = seq(from = 0, to = 1, length.out = row.length)
    # scores = c()
    # for (i in 1:division){
    #   scores = append(scores, dec_to_duosexagesimal(as.bigz(2)^seq(from = i, by = division, length.out = division)))
    # }
    # scores <- c()
    # string.template <- paste0(rep("0", division^2), collapse = "")
    # for (i in 1:(division^2)) {
    #   string <- paste0(substr(string.template, 1, i-1), "1", substr(string.template, i+1, division^2))
    #   scores <- c(scores, string)
    # }
    scores <- 1:(division^2)
  
    for (score in scores){
        row.score = append(row.score, rep(score, 5))
        if (length(row.score)==row.length){
            for(i in 1:5){
                score.vector = append(score.vector, row.score)
                x.vector = append(x.vector, unit.loc)
                y.vector = append(y.vector, rep(unit.loc[quad.count], row.length))
                quad.count = quad.count - 1
            }
            row.score = c()
        }
    }
    matrix$unit_x = x.vector
    matrix$unit_y = y.vector
    matrix$score = score.vector
    print(Sys.time()-start)
    
    #write.csv(matrix, file = paste(output_folder,toString(division), 'by', toString(division), '_matrix.csv', sep = ''), row.names = F)
    
    return(matrix)
}

##########################################################
# get.gc 
##########################################################
get.gc <- function(polygon){
    A <- polygon@area   # Area
    C <- polygon@coords # Coords of each vertix (long, lat)
    C1 <- C[-nrow(C),]  # x_i, y_i
    C2 <- C[-1,]        # x_{i+1}, y_{i+1}
    # Compute geometric centre
    GC <- c(
        sum((C1[,1] + C2[,1]) * (C1[,1]*C2[,2] - C2[,1]*C1[,2])) / (6*A),
        sum((C1[,2] + C2[,2]) * (C1[,1]*C2[,2] - C2[,1]*C1[,2])) / (6*A))
    # If the polygon is oriented clockwise, switch the signs
    if (polygon@ringDir == 1) GC <- -GC
    return(GC)
}


##########################################################
# print_pixel_text_image 
##########################################################
print_pixel_text_image <- function(strings, nchar){
    for(string in strings){
        cat(gsub(paste0("(.{",nchar,"})"), "\\1\n", string))
        if (string != strings[length(strings)]){ cat("\n")}
    }
}


##########################################################
# string_dist 
##########################################################
string_dist <- function(string1, string2){
    if (nchar(string1) != nchar(string2)){
      message("both strings must be the same length")
      return()
    }
    vector1 <- as.numeric(strsplit(string1, split="")[[1]])
    vector2 <- as.numeric(strsplit(string2, split="")[[1]])
    vector <- vector1 + vector2
    # vector_dist <- length(which(!(vector1 == vector2)))
    len_1 <- length(which(vector == 1))
    len_2 <- length(which(vector == 2))
    vector_dist <- round(len_2 / (len_1 + len_2),digits=3)
    return(vector_dist)
}

##########################################################
# in.range 
##########################################################
in.range <- function(x, range) {
    which(x >= range[1] & x <= range[2])
}


##########################################################
# CIDs
##########################################################

calc_pixel_ID <- function(data, polygons, matrix_division, Do_Plots=T, plot_name, x_min, x_max, y_min, y_max){

    x.range.temp <- unname(unlist(c(x_min,x_max)))
    y.range.temp <- unname(unlist(c(y_min,y_max)))
    ind.x <- in.range(data[,1], x.range.temp)
    ind.y <- in.range(data[,2], y.range.temp)
    data <- data[intersect(ind.x, ind.y),]
    
    x.range = abs(diff(x.range.temp))
    y.range = abs(diff(y.range.temp))
    grid.x.range = x.range.temp
    grid.y.range = y.range.temp
    grid.lines.x <- seq(grid.x.range[1],grid.x.range[2],x.range/(matrix_division))
    grid.lines.y <- seq(grid.y.range[1],grid.y.range[2],y.range/(matrix_division))
    center.points.x <- grid.lines.x[-length(grid.lines.x)] + x.range/(matrix_division)/2
    center.points.y <- grid.lines.y[-length(grid.lines.y)] + y.range/(matrix_division)/2
    mesh.points <- mesh(center.points.x, center.points.y)
    
    pixeled_polygons <- polygons_2_pixels(user_polygons = polygons, source = data[,1:2], source_type = 'csv', just_event_label = T, matrix_division = matrix_division, x.range.temp=x.range.temp, y.range.temp=y.range.temp)
    
    pixel.loc <- c()
    string.template <- paste0(rep("0", matrix_division^2), collapse = "")
    for (i in 1:(matrix_division^2)) {
      string <- paste0(substr(string.template, 1, i-1), "1", substr(string.template, i+1, matrix_division^2))
      pixel.loc <- c(pixel.loc, string)
    }
    
    CIDs <- first.pixel <- NULL
    
    for (i in 1:length(pixeled_polygons$pixels)){
      pixel.numbers <- as.numeric(pixeled_polygons$pixels[[i]])
      CIDs[i] <- paste0(ifelse(seq(matrix_division^2) %in% pixel.numbers, "1", "0"), collapse = "") # turn pixel numbers to binary
      first.pixel[i] <- min(pixel.numbers)
    }
    
    CIDs <- CIDs[order(first.pixel)]
    
    if (Do_Plots){
        cexS <- 2
        
        suppressWarnings(dir.create(paste0("results/Annotation/")))
        
        png(paste0("results/Annotation/",plot_name,".png"),width = 600, height = 600+0)
        par(mar=c(1,1,1,1), par(family="sans") )#, oma = c(0,0,4,0)
        
        plot_dens_wrapper(data[,1:2], cex.axis = cexS, cex.lab = cexS, labels=c("",""), axes=F, frame.plot=T, xlim=x.range.temp, ylim=y.range.temp)
        
        col.palette <- c("red","limegreen","blue","purple1","lightblue","orange","palevioletred1", "brown")
        
        colours <- rep("white", matrix_division^2)
        for (i in (1:length(pixeled_polygons$pixels))){
            pixel.numbers <- as.numeric(pixeled_polygons$pixels[[order(first.pixel)[i]]]) # [order(first.pixel)] ensures colours are organized from top left to top right, then bottom.
            ind.double <- which(colours[pixel.numbers] != "white")
            colours[pixel.numbers] <- col.palette[i] 
            if (length(ind.double) > 0){
                colours[pixel.numbers][ind.double] <- "grey75"
            }
        }
        
        x.loc <- as.vector(mesh.points$x)
        y.loc <- rev(as.vector(mesh.points$y))
        x.size <- x.range/(matrix_division)/2
        y.size <- y.range/(matrix_division)/2
        rect(xleft = x.loc-x.size, ybottom = y.loc-y.size, xright = x.loc+x.size, ytop = y.loc+y.size, col = colours, border = NA)
        
        abline(v=grid.lines.x, col="grey60")
        abline(h=grid.lines.y, col="grey60")
        dev.off()
    }
    
    
    return(CIDs)
}

################################################# #
#' plot_dens_wrapper
################################################# #
#'
#' A wrapper for plotting flowFrames that has few checks (e.g. no events) and uses the faster scattermore package
#'
#' @param obj The 2-column express matrix or flowFrame to be plotted.
#' @param channels If obj is a flowFrame, provide the channels to be plotted (a vector of length=2, specifying either the channel names or indices).
#' If obj is a 2-column express matrix, channels can be a named vector with the marker names (channel_order) to print both the laser and marker names on the axes.
#' @param main The title of the plot.
#' @param xlim The range of the X axis. Will be calculated automatically if missing.
#' @param ylim The range of the Y axis. Will be calculated automatically if missing.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of cex.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of cex.
#' @param cex.main The magnification to be used for main titles relative to the current setting of cex.
#' @param show_contour Adding the contour lines to the plot (see \code{\link{flow_contour}}. Default is FALSE.
#' @param contour_col Color for contour lines. Default is 'darkgrey'.
#' @export
#' @examples
#' frame <- example_flowFrame()
#' plot_dens_wrapper(frame,c("APC-H7-A", "PerCP-Cy5-5-A"))
#'
plot_dens_wrapper <- function(obj, channels, main="", xlim, ylim, labels, cex.axis=2, cex.lab=2, cex.main=2, show_contour=FALSE, contour_col="darkgrey", axes=T, frame.plot=T){

    if (is(obj, "flowFrame")) {
        if (nrow(obj) == 0) {
            if (missing(xlim)) {
              xlim <- c(0,10)
            }
            if (missing(ylim)) {
              ylim <- c(0,10)
            }
            if (missing(labels)) {
              labels <- colnames(obj)
            }
            return(graphics::plot(-1000, -1000, # -1000 to be off the plot window
                        main=main, xlab=labels[1], ylab=labels[2],
                        xlim=xlim, ylim=ylim,
                        cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main))
        }
        if (missing(xlim)) {
          xlim <- range(flowCore::exprs(obj)[, channels[1]], na.rm = T)
        }
        if (missing(ylim)) {
          ylim <- range(flowCore::exprs(obj)[, channels[2]], na.rm = T)
        }
        if (nrow(obj) < 10) {
            labels <- make_labels(obj, channels[1], channels[2])

            scattermore::scattermoreplot(x = flowCore::exprs(obj)[,channels[1]], y = flowCore::exprs(obj)[,channels[2]], xlim = xlim, ylim = ylim, col = "black",
                                         cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, main=main, xlab=labels[1], ylab=labels[2], cex = 0, axes=axes, frame.plot=frame.plot)
        } else {
            labels <- make_labels(obj, channels[1], channels[2])

            colPalette <- grDevices::colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
            col <- grDevices::densCols(flowCore::exprs(obj)[, channels], colramp = colPalette)

            scattermore::scattermoreplot(x = flowCore::exprs(obj)[,channels[1]], y = flowCore::exprs(obj)[,channels[2]], xlim = xlim, ylim = ylim, col = col,
                                         cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, main=main, xlab=labels[1], ylab=labels[2], cex = 0, axes=axes, frame.plot=frame.plot)
            if (show_contour) {
              flow_contour(x = obj, y = channels, add=TRUE, col=contour_col, lty=1, lwd=1.5)
            }
        }
    } else {
        if (nrow(obj) == 0) {
          if (missing(xlim)) {
            xlim <- c(0,10)
          }
          if (missing(ylim)) {
            ylim <- c(0,10)
          }
          if (missing(labels)) {
              labels <- colnames(obj)
          }
          return(graphics::plot(-1000, -1000, # -1000 to be off the plot window
                      main=main, xlab=labels[1], ylab=labels[2],
                      xlim=xlim, ylim=ylim,
                      cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main))
        }
        if (missing(xlim)) {
          xlim <- range(obj[,1], na.rm = T)
        }
        if (missing(ylim)) {
          ylim <- range(obj[,2], na.rm = T)
        }
        if (nrow(obj) < 10) {
          if (missing(channels)) {
              if (missing(labels)) {
                  labels <- colnames(obj)
              }
          } else {
            xlabel <- paste("", colnames(obj)[1], "", names(which(channels == colnames(obj)[1])), sep = "")
            ylabel <- paste("", colnames(obj)[2], "", names(which(channels == colnames(obj)[2])), sep = "")
            labels <- c(xlabel, ylabel)
          }
          scattermore::scattermoreplot(x = obj[,1], y = obj[,2], xlim = xlim, ylim = ylim, col = "black",
                                       cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, main=main, xlab=labels[1], ylab=labels[2], cex = 0, axes=axes, frame.plot=frame.plot)
        } else {
          if (missing(channels)) {
              if (missing(labels)) {
                  labels <- colnames(obj)
              }
          } else {
            xlabel <- paste("", colnames(obj)[1], "", names(which(channels == colnames(obj)[1])), sep = "")
            ylabel <- paste("", colnames(obj)[2], "", names(which(channels == colnames(obj)[2])), sep = "")
            labels <- c(xlabel, ylabel)
          }
          colPalette <- grDevices::colorRampPalette(c("blue", "turquoise", "green", "yellow", "orange", "red"))
          col <- grDevices::densCols(obj, colramp = colPalette)
          scattermore::scattermoreplot(x = obj[,1], y = obj[,2], xlim = xlim, ylim = ylim, col = col,
                                       cex.axis=cex.axis, cex.lab=cex.lab, cex.main=cex.main, main=main, xlab=labels[1], ylab=labels[2], cex = 0, axes=axes, frame.plot=frame.plot)
          if (show_contour) {
            flow_contour(x = obj, add=TRUE, col=contour_col, lty=1, lwd=1.5)
          }
        }
    }
}

##################################################
# Time Function                                  #
##################################################
# Prints out the time since startTime
TimePrint <- function(startTime) {
    startTime <- as.POSIXct(startTime, format="%H:%M:%OS3")
    difft <- difftime(Sys.time(), startTime, units="secs")
    format(.POSIXct(difft, tz="GMT"), "%H:%M:%OS3")
}
