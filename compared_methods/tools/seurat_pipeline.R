library(Seurat)
library(dplyr)

runSeurat <- function(obj.sp, obj.sc, cell.names = "mainCtype") {
    ovp.genes <- intersect(rownames(obj.sc), rownames(obj.sp))
    obj.sc <- obj.sc[ovp.genes, ] %>% SCTransform(., assay = "RNA", verbose = TRUE)
    obj.sc <- FindVariableFeatures(obj.sc)
    obj.sp <- obj.sp[ovp.genes, ] %>% SCTransform(., assay = "Spatial", verbose = TRUE)
    anchors <- FindTransferAnchors(reference = obj.sc, query = obj.sp)
    predictions <- TransferData(anchorset = anchors, refdata = obj.sc@meta.data[, cell.names])
    obj.sp <- AddMetaData(object = obj.sp, metadata = predictions)
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

obj.sp <- runSeurat(obj.sp, obj.ref, cell.names = celltype_final)
saveRDS(obj.sp, file = file.path(out.data.dir, "seurat_result.RDS"))
