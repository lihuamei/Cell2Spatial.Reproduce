library(spacexr)
library(Seurat)

out.data.dir <- "../2.results/cerebellum"

ref <- readRDS("../0.data/Cerebellum/cerebellum_SC.RDS")
slide.seq <- readRDS("../0.data/Cerebellum/cerebellum_ST.RDS")
ref <- UpdateSeuratObject(ref)
Idents(ref) <- "liger_ident_coarse"

counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$liger_ident_coarse)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
counts <- slide.seq[["Spatial"]]$counts
coords <- GetTissueCoordinates(slide.seq)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))

RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)
saveRDS(slide.seq, file = file.path(out.data.dir, "slide.seq.RDS"))
