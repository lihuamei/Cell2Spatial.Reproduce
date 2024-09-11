library(dplyr)
library(Seurat)
library(SeuratDisk)
source("../../modules/utils.R")

data.dir <- "../../../2.results/stxBrain/data"

brain.sc <- readRDS(file = "../../../0.data/stxBrain/allen_cortex.rds") %>% SCTransform(verbose = FALSE, assay = "RNA")
brain.st <- readRDS(file = "../../../0.data/stxBrain/brain.st.rds") %>% SCTransform(verbose = FALSE, assay = "Spatial")

simSpots <- function(obj.sp, obj.sc.train, markers, lambda = 10) {
    ovp.genes <- intersect(unlist(markers), rownames(obj.sp))
    # ovp.genes <- intersect(FindVariableFeatures(obj.sp) %>% VariableFeatures, FindVariableFeatures(obj.sc.train) %>% VariableFeatures)
    sp.ref <- obj.sp[ovp.genes, ] %>%
        GetAssayData(., slot = "data", assay = "SCT") %>%
        as.matrix()
    sc.prop <- obj.sc.train[ovp.genes, ] %>%
        GetAssayData(., slot = "data", assay = "SCT") %>%
        as.matrix()
    cor.mat <- cor(sp.ref, sc.prop)
    rand.cnt <- rpois(n = ncol(cor.mat), lambda = lambda) %>%
        `names<-`(rownames(cor.mat)) %>%
        {
            . + 2
        }
    sp.sim <- lapply(rownames(cor.mat), function(sp) {
        cor.vec <- cor.mat[sp, ]
        colnames(cor.mat)[order(cor.vec) %>% rev()][1:rand.cnt[sp]]
    }) %>% `names<-`(rownames(cor.mat))
    return(sp.sim)
}

markers <- readRDS("../../../0.data/sctBrain.markers.RDS")
Idents(brain.sc) <- brain.sc$subclass
brain.sc$cellType <- Idents(brain.sc) %>% as.vector()
obj.sc.train <- downSamplSeurat(brain.sc, cnt = 2000)
obj.st.sub <- brain.st

options(future.globals.maxSize = 5000000 * 1024^2, future.seed = TRUE)
lambda <- 5
sn <- "stxBrain"
# sim.res <- simSpots(brain.st, obj.sc.train, markers, lambda = lambda)
# ref.ctpes <- Idents(obj.sc.train)
# sim.coord <- reformatSimCoord(sim.res, brain.st)
# sim.coord$CellType <- ref.ctpes[sim.coord$CellName]
# saveRDS(sim.coord, file = file.path('./', sprintf('%s.true.coord.rds', sn)))

for (lambda in c(5, 10, 15, 20)) {
    print(sprintf("Lambda lelve: %g", lambda))
    sim.res <- simSpots(brain.st, obj.sc.train, markers, lambda = lambda)
    ref.ctpes <- Idents(obj.sc.train)
    sim.coord <- reformatSimCoord(sim.res, brain.st)
    sim.coord$CellType <- ref.ctpes[sim.coord$CellName]
    saveRDS(sim.coord, file = file.path("./", sprintf("stxBrain.%g.true.coord.rds", lambda)))

    sc.expr <- GetAssayData(obj.sc.train, slot = "count") %>% as.data.frame()
    # if (nois.level > 0) {
    idxes <- sample(nrow(sc.expr), floor(nrow(sc.expr) * 5 * 0.01))
    rand.shuffle <- sc.expr[idxes, ]
    rand.shuffle <- rand.shuffle[, sample(ncol(rand.shuffle))]
    colnames(rand.shuffle) <- colnames(sc.expr)
    sc.expr[idxes, ] <- rand.shuffle
    # }
    sim.expr <- parallel::mclapply(sim.res, function(sim.vec) {
        sc.expr[, sim.vec] %>% rowSums(.)
    }, mc.cores = 10) %>%
        do.call(cbind, .) %>%
        as.data.frame()
    colnames(sim.expr) <- gsub("-", "\\.", colnames(sim.expr))
    obj.sim <- CreateSeuratObject(counts = sim.expr, project = "Sim", assay = "Spatial") %>% SCTransform(verbose = FALSE, assay = "Spatial")
    obj.sim@images <- brain.st@images
    obj.sim@meta.data <- brain.st@meta.data

    rownames(obj.sim@meta.data) <- gsub("-", "\\.", rownames(brain.st@meta.data))
    rownames(obj.sim@images[[1]]@coordinates) <- gsub("-", ".", rownames(obj.sim@images[[1]]@coordinates))
    saveRDS(obj.sim, file = file.path("./", sprintf("%s.sim.%g.rds", sn, lambda)))
}
