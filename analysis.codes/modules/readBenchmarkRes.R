readAsProportions <- function(obj.sc, obj.st, out.dir) {
    ct.raw <- unique(obj.sc$subclass)
    ct.raw <- gsub('/| ', '\\.', ct.raw)
    cell2location <- read.table(file.path(out.dir, 'Cell2location_result.txt'), sep = ',', header = TRUE, row.names = 1)
    colnames(cell2location) <- gsub('q05cell_abundance_w_sf_', '', colnames(cell2location))
    cell2location <- sweep(cell2location, 1, rowSums(cell2location), '/')
    #cell2location <- cell2location[, ctpe.order]
    ctpe.order <- colnames(cell2location)

    rctd <- read.table(file.path(out.dir, 'RCTD_result.txt'), sep = ',', header = TRUE, row.names = 1)
    rctd[, setdiff(ctpe.order, colnames(rctd))] <- 0
    rctd <- rctd[, ctpe.order]

    SPOTlight <- read.table(file.path(out.dir, 'SPOTlight_result.txt'), sep = ',', header = TRUE, row.names = 1)[, ctpe.order] %>% `rownames<-`(rownames(cell2location)) 
    SpatialDWLS <- read.table(file.path(out.dir, 'SpatialDWLS_result.txt'), sep = ',', header = TRUE, row.names = 2) %>% .[, -1] %>% .[, ctpe.order]
    DSTG <- read.table(file.path(out.dir, 'predict_output.csv'), sep = ',', header = F) %>% 
	    `colnames<-`(ct.raw) %>% `rownames<-`(rownames(SpatialDWLS)) %>% .[, ctpe.order]
    DestVI <- read.table(file.path(out.dir, 'DestVI_result.txt'), sep = ',', header = TRUE, row.names = 1)[, ctpe.order]
    Stereoscope <- read.table(file.path(out.dir, 'Stereoscope_result.txt'), sep = ',', header = TRUE, row.names = 1)[, ctpe.order]
    SpaOTsc <- read.table(file.path(out.dir, 'SpaOTsc_decon.csv'), sep = ',', header = TRUE, row.names = 1)[, ctpe.order]
    Tangram <- read.table(file.path(out.dir, 'Tangram_result.txt'), sep = ',', header = TRUE, row.names = 1)[, ctpe.order]
    novoSpaRc <- read.table(file.path(out.dir, 'novoSpaRc_decon.csv'), sep = ',', header = TRUE, row.names = 1)[, ctpe.order] 
    card <- readRDS(file.path(out.dir, 'stxBrain.card.RDS'))@Proportion_CARD
    colnames(card) <- gsub('/| ', '\\.', colnames(card))
    card <- card[, ctpe.order]
    
    cytospace <- read.csv(file.path(out.dir, 'cytospace.stxBrain/fractional_abundances_by_spot.csv'), row.names = 1)
    colnames(cytospace) <- gsub('/| ', '\\.', colnames(cytospace))
    cytospace[, setdiff(colnames(card), colnames(cytospace))] <- 0
    cytospace <- cytospace[, colnames(card)]

    obj <- readRDS(file.path(out.dir, 'cell2spatial.RDS'))
    Cell2Spatial <- table(obj$Cell2Spatial, obj$SpotName) %>% t %>% as.data.frame.matrix
    colnames(Cell2Spatial) <- gsub('/| ', '\\.', colnames(Cell2Spatial))
    Cell2Spatial <- sweep(Cell2Spatial, 1, rowSums(Cell2Spatial), '/')
    Cell2Spatial[, setdiff(colnames(card), colnames(Cell2Spatial))] <- 0
    Cell2Spatial <- Cell2Spatial[, colnames(card)]

    obj <- readRDS(file.path(out.dir, 'stxBrain.celltrek.RDS'))
    obj <- forCellTrek(list(xx = obj), obj.sc, obj.st)
    celltrek <- table(obj$Spot, obj$subclass) %>% as.data.frame.matrix
    colnames(celltrek) <- gsub('/| ', '\\.', colnames(celltrek))
    celltrek <- sweep(celltrek, 1, rowSums(celltrek), '/')
    celltrek[, setdiff(colnames(card), colnames(celltrek))] <- 0
    celltrek <- celltrek[, colnames(card)]
    
    # Seurat 
    obj <- readRDS(file.path(out.dir, 'seurat_result.RDS'))
    seurat <- obj@meta.data[, grep('prediction.score', colnames(obj@meta.data))]
    seurat <- seurat[, -ncol(seurat)]
    colnames(seurat) <- gsub('prediction.score.', '', colnames(seurat))
    colnames(seurat) <- gsub('/| ', '\\.', colnames(seurat))
    seurat <- seurat[, ctpe.order] 
    
    dat <- rbind.data.frame(
    	cbind.data.frame(cell2location, Tool = 'Cell2location', Spot = rownames(cell2location)),
	cbind.data.frame(rctd, Tool = 'RCTD', Spot = rownames(rctd)),
	cbind.data.frame(SPOTlight, Tool = 'SPOTlight', Spot = rownames(SPOTlight)),
    	cbind.data.frame(SpatialDWLS, Tool = 'SpatialDWLS', Spot = rownames(SpatialDWLS)),
	cbind.data.frame(DSTG, Tool = 'DSTG', Spot = rownames(DSTG)),
	cbind.data.frame(DestVI, Tool = 'DestVI', Spot = rownames(DestVI)),
	cbind.data.frame(Stereoscope, Tool = 'Stereoscope', Spot = rownames(Stereoscope)),
	cbind.data.frame(SpaOTsc, Tool = 'SpaOTsc', Spot = rownames(SpaOTsc)),
	cbind.data.frame(Tangram, Tool = 'Tangram', Spot = rownames(Tangram)),
	cbind.data.frame(novoSpaRc, Tool = 'novoSpaRc', Spot = rownames(novoSpaRc)),
	cbind.data.frame(card, Tool = 'CARD', Spot = rownames(card)),
	cbind.data.frame(cytospace, Tool = 'CytoSPACE', Spot = rownames(cytospace)),
	cbind.data.frame(Cell2Spatial, Tool = 'Cell2Spatial', Spot = rownames(Cell2Spatial)),
	cbind.data.frame(celltrek, Tool = 'CellTrek', Spot = rownames(celltrek)),
	cbind.data.frame(seurat, Tool = 'Seurat', Spot = rownames(seurat))
    )    
    return(dat)
}

# Gene expression decomposition

exprDecom <- function(obj.sc, obj.st, out.dir) {
    obj <- readRDS(file.path(out.dir, 'cell2spatial.RDS'))
    sp.lst <- split(obj$Cell, obj$SpotName)	
    sc.mat <- GetAssayData(obj.sc, slot = 'count')
    
    mat <- parallel::mclapply(sp.lst, function(xx) {
	xx <- gsub('\\..*', '', xx)
	if (length(xx) > 1) rowSums(sc.mat[, xx])
	else sc.mat[, xx]
    }, mc.cores = 10) %>% do.call(cbind, .) %>% as.data.frame
    
    obj.cell2spatial <- CreateSeuratObject(count = mat) %>% SCTransform(verbose = FALSE)
    
    obj <- read.table(file.path(out.dir, 'cytospace.stxBrain/assigned_locations.csv'), sep = ',', header = TRUE)
    sp.lst <- split(obj$OriginalCID, obj$SpotID)
    mat <- parallel::mclapply(sp.lst, function(xx) {
        if (length(xx) > 1) rowSums(sc.mat[, xx])
	else sc.mat[, xx]
    }, mc.cores = 10) %>% do.call(cbind, .) %>% as.data.frame
    obj.cytospace <- CreateSeuratObject(count = mat) %>% SCTransform(verbose = FALSE)
    
    obj <- readRDS(file.path(out.dir, 'stxBrain.celltrek.RDS'))
    obj <- forCellTrek(list(xx = obj), obj.sc, obj.st)
    sp.lst <- split(as.vector(obj$id_raw), as.vector(obj$Spot))
    mat <- parallel::mclapply(sp.lst, function(xx) {	      
        if (length(xx) > 1) rowSums(sc.mat[, xx])
	else sc.mat[, xx]
    }, mc.cores = 10) %>% do.call(cbind, .) %>% as.data.frame
    obj.celltrek <- CreateSeuratObject(count = mat) %>% SCTransform(verbose = FALSE)

    obj <- read.table(file.path(out.dir, 'Tangram_map2space.xls'), sep = '\t', header = TRUE) 
    sp.lst <- split(obj$X, obj$SpotID)
    mat <- parallel::mclapply(sp.lst, function(xx) {
        if (length(xx) > 1) rowSums(sc.mat[, xx])
	else sc.mat[, xx]
    }, mc.cores = 10) %>% do.call(cbind, .) %>% as.data.frame
    obj.tangram <- CreateSeuratObject(count = mat) %>% SCTransform(verbose = FALSE)
    return(list(Cell2Spatial = obj.cell2spatial, CytoSPACE = obj.cytospace, Celltrek = obj.celltrek, Tangram = obj.tangram))
} 


readMappingData <- function(obj.sc, obj.st, sim.true, out.dir) { 	
     obj <- read.table(file.path(out.dir, 'cytospace.stxBrain/assigned_locations.csv'), sep = ',', header = TRUE)	
     obj$SpotID <- gsub('\\.', '-', obj$SpotID)
     obj$CellName <- sim.true[obj$OriginalCID, 'CellNameNew']
     cytospace <- split(obj$CellName, obj$SpotID)
     
     obj <- readRDS(file.path(out.dir, 'stxBrain.celltrek.RDS'))
     #obj <- forCellTrek(list(xx = obj), obj.sc, obj.st)
     obj$SpotID <- gsub('\\.', '-', obj$SpotID)
     celltrek <- split(as.vector(obj$CellNameNew), as.vector(obj$SpotID))
    
     obj <- read.table(file.path(out.dir, 'Tangram_map2space.xls'), sep = '\t', header = TRUE)
     obj$SpotID <- gsub('\\.', '-', obj$SpotID)
     tangram <- split(obj$CellNameNew, obj$SpotID)
     
     obj <- readRDS(file.path(out.dir, 'cell2spatial.unmatched.RDS'))
     obj$SpotName <- gsub('\\.', '-', obj$SpotName)
     #obj$CellNameNew <- gsub('Pseudo_|_XYZ.*|_F\\d+', '', obj$Cell)
     obj$CellNameNew <- gsub('Pseudo_|_XYZ.*', '', obj$Cell)
     cell2spatial <- split(obj$CellNameNew, obj$SpotName)
     return(list(Cell2Spatial = cell2spatial, CytoSPACE = cytospace, CellTrek = celltrek, Tangram = tangram))
}

forTangramSim <- function (dir.paths, obj.sc, obj.sp) 
{
    sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% 
        `colnames<-`(c("row", "col"))
    proj.res <- lapply(dir.paths, function(dir) {
        meta.info <- read.csv(dir, sep = "\t", header = TRUE)
        #meta.info$SpotID <- gsub("\\.", "-", meta.info$SpotID)
        meta.info <- cbind.data.frame(meta.info, sp.loc[meta.info$SpotID, 
            ])
        meta.info <- getRandomCords(meta.info)
    })[[1]]
    count.sc <- GetAssayData(obj.sc, slot = "count")
    count.CT <- as(count.sc[, proj.res$X], "sparseMatrix")
    rownames(proj.res) <- make.unique(proj.res$X)
    colnames(count.CT) <- rownames(proj.res)
    sce <- CreateSeuratObject(count = count.CT, meta.data = proj.res, 
        project = "Cell2Spatial", assay = "Spatial")
    sce@images <- obj.sp@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = proj.res$x.new, 
        imagecol = proj.res$y.new) %>% `rownames<-`(rownames(proj.res))
    sce@images[[1]]@scale.factors <- obj.sp@images[[1]]@scale.factors
    sce@images[[1]]@coordinates <- sce@images[[1]]@coordinates/sce@images[[1]]@scale.factors$lowres
    return(sce)
}
