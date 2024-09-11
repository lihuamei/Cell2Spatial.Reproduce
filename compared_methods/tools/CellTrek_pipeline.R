library(Seurat)
library(dplyr)
devtools::load_all("../.tmp/CellTrek-main")

# out.data.dir <- '../2.results/stxBrain/data'

runCellTrek <- function(obj.sp, obj.sc, cell.names = "CellType", ...) {
    # obj.sc.train <- readRDS('data/end.to.end/obj.sc.train.RDS')
    obj.sc.train <- readRDS("data/allen_cortex.rds")
    # obj.sc.train@meta.data[, cell.names] <- gsub('/| ', '.', obj.sc.train@meta.data[, cell.names])
    ovp.genes <- intersect(rownames(obj.sc.train), rownames(obj.sc))
    obj.sc <- obj.sc[ovp.genes, ]
    # obj.sc.train <- obj.sc.train[ovp.genes, ]
    # obj.sc.train <- merge(obj.sc.train, y = obj.sc, add.cell.ids=NULL)
    obj.traint <- traint(st_data = obj.sp, sc_data = obj.sc, sc_assay = "RNA", cell_names = cell.names)
    obj.celltrek <- celltrek(
        st_sc_int = obj.traint,
        int_assay = "traint",
        sc_data = obj.sc,
        sc_assay = "RNA",
        reduction = "pca",
        intp = TRUE,
        intp_pnt = 50,
        intp_lin = F,
        nPCs = 30,
        ntree = 10,
        dist_thresh = 0.7,
        top_spot = 1,
        spot_n = 10,
        repel_r = 10,
        repel_iter = 10,
        keep_model = TRUE,
        ...
    )$celltrek
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
obj.ref@meta.data[, celltype_final] <- gsub("/| ", ".", obj.ref@meta.data[, celltype_final])
print(obj.ref@meta.data[, celltype_final])

obj.sp <- runCellTrek(obj.sp, obj.ref, cell.names = celltype_final)
saveRDS(obj.sp, file = file.path(output_path, sprintf("%s.celltrek.RDS", "stxBrain")))
