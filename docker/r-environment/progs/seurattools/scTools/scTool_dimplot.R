scTool_dimplot <- function(obj, 
                           reduction_type, 
                           group.var=NULL, 
                           filenameout=NULL, 
                           setwidth=1200, 
                           setheight=800, 
                           pt.size=1, 
                           imgtitle=NULL){
  
  try(if(is.null(group.var)) stop("Please set group.var!"))
  try(if(is.null(filenameout)) stop("Please set filenameout!"))
  try(if(is.null(imgtitle)) stop("Please provide a title for the chart"))
  
  png(filename = filenameout, width = setwidth, height = setheight)
  print(DimPlot(object = obj, reduction.use = reduction_type, group.by = group.var, 
                pt.size = pt.size, do.return = TRUE, plot.title = imgtitle))
  dev.off()
}
