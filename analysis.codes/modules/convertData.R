library(Seurat)
library(SeuratDisk)
library(dplyr)

sc_data <- "../../../0.data/stxBrain/allen_cortex.rds"
obj.sc <- readRDS(sc_data)
prefix_name_sc <- gsub("\\.rds|\\.RDS", "", basename(sc_data))

# st_data <- '../../../2.results/stxBrain/data/stxBrain.sim.0.rds'
st_data <- "stxBrain.sim.20.rds"
obj.st <- readRDS(st_data)
prefix_name_st <- gsub("\\.rds|\\.RDS", "", basename(st_data))

celltyep <- "subclass"

# SaveH5Seurat(UpdateSeuratObject(obj.sc), filename = sprintf("%s.h5Seurat", prefix_name_sc), overwrite = TRUE)
# Convert(sprintf("%s.h5Seurat", prefix_name_sc), dest = "h5ad", overwrite = TRUE)

DefaultAssay(obj.st) <- "Spatial"
coords <- GetTissueCoordinates(obj.st)
obj.st@meta.data[, c("X", "Y")] <- coords[colnames(obj.st), ]

SaveH5Seurat(obj.st, filename = sprintf("%s.h5Seurat", prefix_name_st), overwrite = TRUE)
Convert(sprintf("%s.h5Seurat", prefix_name_st), dest = "h5ad", overwrite = TRUE)
saveRDS(obj.st, file = sprintf("%s.rds", prefix_name_sc))

# write.table(GetAssayData(obj.sc, slot = 'count'), file = sprintf("%s.txt", prefix_name_sc), row.names = TRUE, col.names = NA)
# write.table(GetAssayData(obj.st, slot = 'count'), file = sprintf("%s.txt", prefix_name_st), row.names = TRUE, col.names = NA)

# cell_file <- obj.sc@meta.data[, celltyep, drop = F]
# write.table(cell_file, file = sprintf("%s.txt", prefix_name_sc), row.names = TRUE, col.names = FALSE)
