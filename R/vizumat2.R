
#' Visualize matrix data in a form of a heatmap, with categorical values legend 
#' 
#' @description  
#' Matrix data visualization in a form of a heatmap, with
#' the use of \code{ggplot2} library. Numerical values are represented as
#' categorical. Minimum user input (a matrix object) is needed to produce decent visualization output. Further
#' plot adjustments are available, including tile color change,
#' adding a title, font size change, axis label clearing and others. 
#' 
#' @param matrix.object matrix 
#' @param title plot title
#' @param base_size base font size
#' @param scale_fill_manual.values vector of legend colors for categorical values
#' @param geom_tile.colour tiles color value
#' @param clear.labels logical whether or not clear both x- and y-axis labels         
#' @param clear.x.label logical whether or not clear x-axis labels      
#' @param clear.y.label logical whether or not clear y-axis labels         
#' @param uniform.labes logical whether or not define generic short column and rows labeling:
#' * 'c1','c2',...,'cp' for columns,
#' * 'r1','r2',...,'rp' for rows;
#' might be especially useful if the matrix some long colnames and rownames already assigned      
#' @param rotate.x.labels logical whether or not rotate x-axis labels by 90 degrees 
#' @param x.lab x-axis label
#' @param y.lab y-axis label
#' @param axis.text.x.size font size of x-axis text
#' @param axis.text.y.size font size of y-axis text
#' @param axis.title.x.size font size of x-axis label
#' @param axis.title.y.size font size of y-axis label
#' @param legend.text.size  font size of legend text
#' @param legend.title.size font size of legend title
#' @param legend.title legend title    
#' @param text.font.family  font family 
#' @param remove.legend logical whether or not remove legend
#' @param factor.levels vector of values defining levels of factors 
#'                        (might be used to redefine order of variables in the legend)
#' @param axis.text.x.breaks.idx indices of x-axis elements whose thicks are kept 
#'                                and whose numerical labels are kept
#' @param axis.text.y.breaks.idx indices of y-axis elements whose thicks are kept 
#'                                and whose numerical labels are kept
#' @md
#' 
#' @return \code{ggplot2} object
#' 
#' @examples 
#' mat <- diag(30)
#' vizu.mat.factor(mat)
#' vizu.mat.factor(mat, 
#'                 title = "some title",
#'                 scale_fill_manual.values = c("white","red"),
#'                 axis.text.x.breaks.idx = seq(1,30,5),
#'                 axis.text.y.breaks.idx = seq(1,30,5))
#' vizu.mat.factor(mat, 
#'                 title = "some title: large font, legend: small font",
#'                 base_size = 20, 
#'                 legend.text.size  = 10, 
#'                 legend.title.size = 10)
#' vizu.mat.factor(mat, 
#'                 scale_fill_manual.values = c("white","red"),
#'                 clear.labels = FALSE) 
#' colnames(mat) <- paste0("col", 1:30, sample(LETTERS, 30, replace = TRUE))
#' rownames(mat) <- paste0("row", 1:30, sample(LETTERS, 30, replace = TRUE))
#' vizu.mat.factor(mat, 
#'                 clear.labels = FALSE,
#'                 rotate.x.labels = TRUE) 
#' 
#' @import reshape2
#' @import ggplot2
#' @export
#' 
vizu.mat.factor <- function(matrix.object, 
                            title = "", 
                            base_size = 12, 
                            scale_fill_manual.values = NULL,
                            geom_tile.colour = "grey90",
                            clear.labels    = TRUE, 
                            clear.x.label   = FALSE, 
                            clear.y.label   = FALSE, 
                            uniform.labes   = FALSE, 
                            rotate.x.labels = FALSE, 
                            x.lab = "",
                            y.lab = "",
                            axis.text.x.size  = base_size-2, 
                            axis.text.y.size  = base_size-2,
                            axis.title.x.size = base_size-2,
                            axis.title.y.size = base_size-2,
                            legend.text.size  = base_size-2, 
                            legend.title.size = base_size-2,
                            legend.title      = "value", 
                            text.font.family  = "Helvetica",
                            remove.legend = FALSE,
                            factor.levels = NULL,
                            axis.text.x.breaks.idx = NULL,
                            axis.text.y.breaks.idx = NULL) 
{
  n.col <- ncol(matrix.object)
  n.row <- nrow(matrix.object)
  
  # Define column names if not defined or uniform.labes == TRUE
  if (is.null(colnames(matrix.object)) | uniform.labes) {
    colnames(matrix.object) <- paste0("c", 1:n.col)
  }
  # Define row names if not defined or uniform.labes == TRUE
  if (is.null(rownames(matrix.object)) | uniform.labes) {
    rownames(matrix.object) <- paste0("r", 1:n.row)
  }
  
  # Reshape matrix object 
  matrix.object.m <- melt(matrix.object)
  matrix.object.m[, "Var1"] <- factor(matrix.object.m[, "Var1"], 
                                      levels = rev(rownames(matrix.object)))
  matrix.object.m[, "Var2"] <- factor(matrix.object.m[, "Var2"], 
                                      levels = colnames(matrix.object))
  
  # Define numerical matrix values as factors; define factor levels if factor.levels supplied
  if(!is.null(factor.levels)){
    matrix.object.m$value <- factor(matrix.object.m$value, levels = factor.levels)
  } else {
    matrix.object.m$value <- factor(matrix.object.m$value)
  }

  # ggplot object skeleton 
  plot.tmp <- 
    ggplot(matrix.object.m, aes_string("Var2", "Var1")) + 
    geom_tile(aes_string(fill = "value"), colour = geom_tile.colour) + 
    labs(x     = x.lab, 
         y     = y.lab, 
         title = title, 
         fill  = legend.title) + 
    theme_grey(base_size   = base_size,                           
               base_family = text.font.family) + 
    theme(legend.text  = element_text(size = legend.text.size),  
          legend.title = element_text(size = legend.title.size)) 
  
  # Define tiles fill color if scale_fill_manual.values supplied 
  if(!is.null(scale_fill_manual.values)){
    plot.tmp <- plot.tmp + scale_fill_manual(values = scale_fill_manual.values)
  }

  # Rotate x.axis labels labels
  if (rotate.x.labels) {
    plot.tmp <- plot.tmp + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Define labels size 
  plot.tmp <- plot.tmp + 
    theme(axis.text.x = element_text(size = axis.text.x.size),
          axis.text.y = element_text(size = axis.text.y.size),
          axis.title.x = element_text(size = axis.title.x.size),
          axis.title.y = element_text(size = axis.title.y.size))
  
  # x-axis thicks and labels management
  if (!is.null(axis.text.x.breaks.idx)){
    plot.tmp <- plot.tmp + scale_x_discrete(labels = as.character(axis.text.x.breaks.idx),
                                            breaks = colnames(matrix.object)[axis.text.x.breaks.idx])
  } else {
    if (clear.labels | clear.x.label){
      plot.tmp <- plot.tmp + scale_x_discrete(breaks = NULL)
    }
  }
  
  # y-axis thicks and labels management
  if (!is.null(axis.text.y.breaks.idx)){
    plot.tmp <- plot.tmp + scale_y_discrete(labels = as.character(axis.text.y.breaks.idx),
                                            breaks = rownames(matrix.object)[axis.text.y.breaks.idx])
  } else {
    if (clear.labels | clear.y.label){
      plot.tmp <- plot.tmp + scale_y_discrete(breaks = NULL)
    }
  }
  
  # Remove legend
  if (remove.legend){
    plot.tmp <- plot.tmp + theme(legend.position = "none")
  }
  return(plot.tmp)
}




# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------

#' Visualize matrix data in a form of a heatmap, with continuous values legend 
#' 
#' @description  
#' Matrix data visualization in a form of a heatmap, with
#' the use of \code{ggplot2} library. Minimum user input (a matrix object) is needed to produce decent visualization output. 
#' Automatic plot adjustments are implemented and used as defaults, including 
#' selecting legend color palette and legend scale limits.  Further
#' plot adjustments are available, including 
#' adding a title, font size change, axis label clearing and others. 
#' 
#' @param matrix.object matrix 
#' @param title plot title
#' @param base_size base font size
#' @param adjust.limits logical whether or not adjust legend scale limits automatically:
#' * legend scale starts / ends with 0 for matrix with non-negative / non-positive values only,
#' * legend scale is symmetric for matrix with both negative and positive values
#' @param adjust.colors logical whether or not adjust legend color automatically: 
#' * legend color palette white-red for a data matrix with non-negative values only,
#' * legend color palette blue-white for a data matrix with non-positive values only,
#' * legend color palette blue-white-red for a data matrix with both positive and negative values
#' @param fill.scale.limits 2-element vector defining legend scale limits
#' @param colors.palette legend color color palette
#' @param geom_tile.colour tiles color value
#' @param clear.labels logical whether or not clear both x- and y-axis labels         
#' @param clear.x.label logical whether or not clear x-axis labels      
#' @param clear.y.label logical whether or not clear y-axis labels         
#' @param uniform.labes logical whether or not define generic short column and rows labeling:
#' * 'c1','c2',...,'cp' for columns,
#' * 'r1','r2',...,'rp' for rows;
#' might be especially useful if the matrix some long colnames and rownames already assigned   
#' @param rotate.x.labels logical whether or not rotate x-axis labels by 90 degrees 
#' @param x.lab x-axis label
#' @param y.lab y-axis label
#' @param axis.text.x.size font size of x-axis text
#' @param axis.text.y.size font size of y-axis text
#' @param axis.title.x.size font size of x-axis label
#' @param axis.title.y.size font size of y-axis label
#' @param legend.text.size  font size of legend text
#' @param legend.title.size font size of legend title
#' @param legend.title legend title    
#' @param text.font.family  font family 
#' @param remove.legend logical whether or not remove legend
#' @param axis.text.x.breaks.idx indices of x-axis elements whose thicks are kept 
#'                                and whose numerical labels are kept
#' @param axis.text.y.breaks.idx indices of y-axis elements whose thicks are kept 
#'                                and whose numerical labels are kept
#' @md
#' 
#' @return \code{ggplot2} object
#' 
#' @examples 
#' mat <- matrix(rnorm(30*30), nrow = 30, ncol = 30)
#' vizu.mat(mat)
#' vizu.mat(mat, fill.scale.limits = c(-3,3))
#' vizu.mat(mat, fill.scale.limits = c(-10,10))
#' vizu.mat(mat, fill.scale.limits = c(-10,10), 
#'          uniform.labes = TRUE, clear.labels = FALSE)
#' colnames(mat) <- paste0("col", 1:30, sample(LETTERS, 30, replace = TRUE))
#' rownames(mat) <- paste0("row", 1:30, sample(LETTERS, 30, replace = TRUE))
#' vizu.mat(mat, fill.scale.limits = c(-10,10), 
#'          clear.labels = FALSE, 
#'          rotate.x.labels = TRUE)
#' mat.positive <- abs(mat)
#' vizu.mat(mat.positive, 
#'          title = "positive values only -> legend limits and colors automatically adjusted",
#'          clear.labels = FALSE, 
#'          rotate.x.labels = TRUE)
#' 
#' @import reshape2
#' @import ggplot2
#' @export
#' 
vizu.mat <- function(matrix.object, 
                     title = "", 
                     base_size = 12, 
                     adjust.limits     = TRUE, 
                     adjust.colors     = TRUE, 
                     fill.scale.limits = NULL,
                     colors.palette    = NULL,
                     geom_tile.colour  = "grey90",
                     clear.labels    = TRUE, 
                     clear.x.label   = FALSE, 
                     clear.y.label   = FALSE, 
                     uniform.labes   = FALSE, 
                     rotate.x.labels = FALSE, 
                     x.lab = "",
                     y.lab = "",
                     axis.text.x.size  = base_size-2, 
                     axis.text.y.size  = base_size-2,
                     axis.title.x.size = base_size-2,
                     axis.title.y.size = base_size-2,
                     legend.text.size  = base_size-2, 
                     legend.title.size = base_size-2,
                     legend.title      = "value", 
                     text.font.family  = "Helvetica",
                     remove.legend = FALSE,
                     axis.text.x.breaks.idx = NULL,
                     axis.text.y.breaks.idx = NULL) 
{
  n.col <- ncol(matrix.object)
  n.row <- nrow(matrix.object)
  
  # Define column names if not defined or uniform.labes == TRUE
  if (is.null(colnames(matrix.object)) | uniform.labes) {
    colnames(matrix.object) <- paste0("c", 1:n.col)
  }
  # Define row names if not defined or uniform.labes == TRUE
  if (is.null(rownames(matrix.object)) | uniform.labes) {
    rownames(matrix.object) <- paste0("r", 1:n.row)
  }
  
  # Reshape matrix object 
  matrix.object.m <- melt(matrix.object)
  matrix.object.m[, "Var1"] <- factor(matrix.object.m[, "Var1"], 
                                      levels = rev(rownames(matrix.object)))
  matrix.object.m[, "Var2"] <- factor(matrix.object.m[, "Var2"], 
                                      levels = colnames(matrix.object))
  
  # ggplot object skeleton 
  plot.tmp <- 
    ggplot(matrix.object.m, aes_string("Var2", "Var1")) + 
    geom_tile(aes_string(fill = "value"), colour = geom_tile.colour) + 
    labs(x     = x.lab, 
         y     = y.lab, 
         title = title, 
         fill  = legend.title) + 
    theme_grey(base_size   = base_size,                           
               base_family = text.font.family) + 
    theme(legend.text  = element_text(size = legend.text.size),  
          legend.title = element_text(size = legend.title.size)) 
  
  # Tiles fill color
  mat.r <- range(matrix.object, na.rm = T)
  if (adjust.limits) {
    if (mat.r[1] < 0 & mat.r[2] > 0) {
      limits.tmp <- c(-1, 1)*max(abs(matrix.object)) 
      colours.tmp <- c("blue", "white", "red")
    }
    else if (mat.r[2] > 0) {
      limits.tmp <- c(0, mat.r[2])
      colours.tmp <- c("white", "red")
    }
    else {
      limits.tmp <- c(mat.r[1], 0)
      colours.tmp <- c("blue", "white")
    }
  }
  else {
    limits.tmp <- NULL
  }
  if (!is.null(fill.scale.limits)) {
    limits.tmp <- fill.scale.limits
  }
  if (!is.null(colors.palette)) {
    colours.tmp <- colors.palette
  }
  else if (adjust.colors) {
    if (mat.r[1] < 0 & mat.r[2] > 0) {
      colours.tmp <- c("blue", "white", "red")
    }
    else if (mat.r[2] > 0) {
      colours.tmp <- c("white", "red")
    }
    else {
      colours.tmp <- c("blue", "white")
    }
  }
  else {
    colours.tmp <- NULL
  }
  plot.tmp <- plot.tmp + scale_fill_gradientn(colours = colours.tmp, limits = limits.tmp)
  
  # Rotate x.axis labels labels
  if (rotate.x.labels) {
    plot.tmp <- plot.tmp + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
  }
  
  # Define labels size 
  plot.tmp <- plot.tmp + 
    theme(axis.text.x = element_text(size = axis.text.x.size),
          axis.text.y = element_text(size = axis.text.y.size),
          axis.title.x = element_text(size = axis.title.x.size),
          axis.title.y = element_text(size = axis.title.y.size))
  
  # x-axis thicks and labels management
  if (!is.null(axis.text.x.breaks.idx)){
    plot.tmp <- plot.tmp + scale_x_discrete(labels = as.character(axis.text.x.breaks.idx),
                                            breaks = colnames(matrix.object)[axis.text.x.breaks.idx])
  } else {
    if (clear.labels | clear.x.label){
      plot.tmp <- plot.tmp + scale_x_discrete(breaks = NULL)
    }
  }
  
  # y-axis thicks and labels management
  if (!is.null(axis.text.y.breaks.idx)){
    plot.tmp <- plot.tmp + scale_y_discrete(labels = as.character(axis.text.y.breaks.idx),
                                            breaks = rownames(matrix.object)[axis.text.y.breaks.idx])
  } else {
    if (clear.labels | clear.y.label){
      plot.tmp <- plot.tmp + scale_y_discrete(breaks = NULL)
    }
  }
  
  # Remove legend
  if (remove.legend){
    plot.tmp <- plot.tmp + theme(legend.position = "none")
  }
  return(plot.tmp)
}



