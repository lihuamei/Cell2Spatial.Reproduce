library(Seurat)
library(dplyr)

out.data.dir <- "../2.results/stxBrain/data"

labels <- "stxBrain"
obj.ref <- readRDS("../0.data/stxBrain/allen_cortex.rds")
obj.sp <- readRDS("../0.data/stxBrain/brain.st.rds")

obj.sp <- UpdateSeuratObject(obj.sp)
SeuratDisk::SaveH5Seurat(obj.sp, filename = ".tmp/obj.sp.h5Seurat", overwrite = TRUE)
SeuratDisk::Convert(".tmp/obj.sp.h5Seurat", dest = "h5ad", overwrite = TRUE)

obj.ref <- UpdateSeuratObject(obj.ref)
SeuratDisk::SaveH5Seurat(obj.ref, filename = ".tmp/obj.sc.h5Seurat", overwrite = TRUE)
SeuratDisk::Convert(".tmp/obj.sc.h5Seurat", dest = "h5ad", overwrite = TRUE)
system(sprintf("python pbs/runTangram.py .tmp/obj.sc.h5ad .tmp/obj.sp.h5ad %s %s", labels, out.data.dir))
