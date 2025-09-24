library(Seurat)
library(dplyr)

source("../modules/utils.R")
# out.data.dir <- '../2.results/stxBrain/data'

runCytoSpace <- function(obj.sp, obj.sc, cell.names = "mainCtype", out.dir = "./") {
    file.res <- generateCytoSPACE(obj.sc, obj.sp, cell.names, "cytospace", out.dir)
    cmd <- sprintf(
        "cytospace --scRNA-path %s --cell-type-path %s --st-path %s --coordinates-path %s -o %s ",
        file.res$sc,
        file.res$ct,
        file.res$st,
        file.res$pos,
        out.dir
    )
    system("echo 'source /root/miniconda3/bin/activate cytospace' > .cytoscape.sh")
    system(sprintf("echo %s >> .cytoscape.sh", cmd))
    system("bash .cytoscape.sh")
}

# labels <- 'stxBrain'
# obj.ref <- readRDS("../0.data/stxBrain/allen_cortex.rds")
# obj.sp <- readRDS("../0.data/stxBrain/brain.st.rds")

args <- commandArgs(T)
scrna_path <- args[1]
spatial_path <- args[2]
celltype_final <- args[3]
output_path <- args[4]

obj.sp <- readRDS(spatial_path)
DefaultAssay(obj.sp) <- "Spatial"
obj.ref <- readRDS(scrna_path)
labels <- gsub("\\.RDS|\\.rds", "", basename(spatial_path))

out.dir <- file.path(output_path, sprintf("cytospace.%s", labels))
dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
obj.sp <- runCytoSpace(obj.sp, obj.ref, cell.names = celltype_final, out.dir = out.dir)
