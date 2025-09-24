library(Seurat)
library(dplyr)
library(CARD)

# out.data.dir <- '../2.results/stxBrain/data'

runCARD <- function(obj.sp, obj.sc, cell.names = "mainCtype", num.cell = 10) {
    sp.counts <- GetAssayData(obj.sp, slot = "count")
    sc.counts <- GetAssayData(obj.sc, slot = "count")
    sc.meta <- obj.sc@meta.data
    sc.meta$cellType <- as.vector(obj.sc@meta.data[, cell.names])
    sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c("x", "y"))
    sp.counts <- sweep(sp.counts, 2, max(colSums(sp.counts)), "/") * 1e6
    sp.counts <- apply(sp.counts, 1, as.integer) %>% t()
    colnames(sp.counts) <- colnames(obj.sp)

    CARD.obj <- createCARDObject(
        sc_count = sc.counts,
        sc_meta = sc.meta,
        spatial_count = sp.counts,
        spatial_location = sp.loc,
        ct.varname = cell.names,
        ct.select = unique(sc.meta$cellType),
        sample.varname = "orig.ident",
        minCountGene = 100,
        minCountSpot = 5
    )
    CARD.dec <- CARD_deconvolution(CARD_object = CARD.obj)
    # CARD.obj <- CARD_SCMapping(CARD.dec, numCell = num.cell)
    # return(list(obj = CARD.obj, decon = CARD.dec))
    return(CARD.dec)
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

obj.sp <- runCARD(obj.sp, obj.ref, cell.names = "subclass")
saveRDS(obj.sp, file = file.path(output_path, sprintf("%s.card.RDS", labels)))
