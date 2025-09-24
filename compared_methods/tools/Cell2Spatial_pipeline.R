library(Seurat)
library(dplyr)

devtools::load_all("Cell2Spatial/")

runSeurat <- function(obj.sp, obj.sc, cell.name = "mainCtype", max.cells.in.spot = 10) {
    sce <- runCell2Spatial(obj.sp, obj.sc, ctype = cell.name, res = 0.8, group.size = 30, max.cells.in.spot = max.cells.in.spot)
}

args <- commandArgs(T)
scrna_path <- args[1]
spatial_path <- args[2]
celltype_final <- args[3]
output_path <- args[4]

obj.ref <- readRDS(scrna_path)
obj.sp <- readRDS(spatial_path)
cell.names <- celltype_final
out.data.dir <- output_path

obj.sp <- runSeurat(obj.sp, obj.ref, cell.name = cell.names, max.cells.in.spot = max.cells.in.spot)
saveRDS(obj.sp, file = file.path(out.data.dir, "cell2spatial.unmatched.RDS"))
