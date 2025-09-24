forCell2Spatial <- function(obj.lst) {
    proj.res <- lapply(obj.lst, function(obj) {
        meta.data <- SummarizedExperiment::colData(obj)
    })
    return(proj.res)
}

forCARD <- function(obj.lst) {
    proj.res <- lapply(obj.lst, function(obj) {
        meta.data <- SummarizedExperiment::colData(obj)
    })[[1]]
    return(proj.res)
}

forCellTrek <- function(obj.lst, obj.sc, obj.sp) {
    sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c("x", "y"))
    proj.res <- lapply(1:length(obj.lst), function(kk) {
        celltrek.loc <- GetTissueCoordinates(obj.lst[[kk]])
        map.sp <- parallel::mclapply(1:nrow(celltrek.loc), function(idx) {
            coord.1 <- celltrek.loc[idx, ] %>% unlist()
            dist.vec <- apply(sp.loc, 1, function(coord.2) {
                df <- cbind.data.frame(coord.1, coord.2) %>% t()
                dist(df)
            }) %>%
                {
                    rownames(sp.loc)[order(.)[1]]
                }
        }, mc.cores = 1) %>%
            unlist() %>%
            `names<-`(rownames(celltrek.loc))

        tmp.loc <- sp.loc[map.sp, ] %>% `colnames<-`(c("centerX", "centerY"))
        meta.data <- cbind.data.frame(obj.lst[[kk]]@meta.data, Spot = map.sp, tmp.loc, celltrek.loc)
    })[[1]]
    count.sc <- GetAssayData(obj.sc, slot = "count")
    count.CT <- as(count.sc[, proj.res$id_raw], "sparseMatrix")
    rownames(proj.res) <- make.unique(proj.res$id_raw)
    colnames(count.CT) <- rownames(proj.res)

    sce <- CreateSeuratObject(count = count.CT, meta.data = proj.res, project = "Cell2Spatial", assay = "Spatial")
    sce@images <- obj.sp@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = proj.res$imagerow, imagecol = proj.res$imagecol) %>% `rownames<-`(rownames(proj.res))
    sce@images[[1]]@scale.factors <- obj.sp@images[[1]]@scale.factors
    sce@images[[1]]@coordinates <- sce@images[[1]]@coordinates / sce@images[[1]]@scale.factors$lowres
    return(sce)
}

forSeurat <- function(obj.lst) {
    proj.res <- lapply(obj.lst, function(obj) {
        meta.data <- cbind.data.frame(obj@meta.data, GetTissueCoordinates(obj))
    })[[1]]
    return(proj.res)
}

forCytoSpace <- function(dir.paths, obj.sc, obj.st) {
    proj.res <- lapply(dir.paths, function(dir) {
        meta.info <- read.table(file.path(dir, "assigned_locations.csv"), sep = ",", header = TRUE)
        meta.info <- getRandomCords(meta.info, n.workers = 1)
    })[[1]]
    count.sc <- GetAssayData(obj.sc, slot = "count")
    count.CT <- as(count.sc[, proj.res$OriginalCID], "sparseMatrix")
    colnames(count.CT) <- proj.res$UniqueCID
    rownames(proj.res) <- proj.res$UniqueCID

    sce <- CreateSeuratObject(count = count.CT, meta.data = proj.res, project = "Cell2Spatial", assay = "Spatial")
    sce@images <- obj.st@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = proj.res$x.new, imagecol = proj.res$y.new) %>% `rownames<-`(rownames(proj.res))
    sce@images[[1]]@scale.factors <- obj.st@images[[1]]@scale.factors
    sce@images[[1]]@coordinates <- sce@images[[1]]@coordinates / sce@images[[1]]@scale.factors$lowres
    return(sce)
}

forTangram <- function(dir.paths, obj.sc, obj.sp) {
    sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c("row", "col"))
    proj.res <- lapply(dir.paths, function(dir) {
        meta.info <- read.csv(dir, sep = "\t", header = TRUE)
        meta.info$SpotID <- gsub("\\.", "-", meta.info$SpotID)
        meta.info <- cbind.data.frame(meta.info, sp.loc[meta.info$SpotID, ])
        meta.info <- getRandomCords(meta.info, n.workers = 1)
    })[[1]]
    count.sc <- GetAssayData(obj.sc, slot = "count")
    count.CT <- as(count.sc[, proj.res$X], "sparseMatrix")
    rownames(proj.res) <- make.unique(proj.res$X)
    colnames(count.CT) <- rownames(proj.res)

    sce <- CreateSeuratObject(count = count.CT, meta.data = proj.res, project = "Cell2Spatial", assay = "Spatial")
    sce@images <- obj.sp@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = proj.res$x.new, imagecol = proj.res$y.new) %>% `rownames<-`(rownames(proj.res))
    sce@images[[1]]@scale.factors <- obj.sp@images[[1]]@scale.factors
    sce@images[[1]]@coordinates <- sce@images[[1]]@coordinates / sce@images[[1]]@scale.factors$lowres

    return(sce)
}

forSeuratNew <- function(obj.lst, obj.sp) {
    obj <- obj.lst[[1]]
    sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c("row", "col"))
    proj.res <- lapply(obj.lst, function(obj) {
        meta.info <- obj@meta.data
        meta.info$SpotID <- meta.info$predicted.id
        meta.info <- cbind.data.frame(meta.info, sp.loc[meta.info$SpotID, ])
        meta.info <- getRandomCords(meta.info)
        # meta.data <- cbind.data.frame(obj@meta.data, GetTissueCoordinates(obj))
    })[[1]]
    count.sc <- GetAssayData(obj, slot = "count")
    sce <- CreateSeuratObject(count = count.sc, meta.data = proj.res, project = "Seurat", assay = "Spatial")
    sce@images <- obj.sp@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = proj.res$x.new, imagecol = proj.res$y.new) %>% `rownames<-`(rownames(proj.res))
    sce@images[[1]]@scale.factors <- obj.sp@images[[1]]@scale.factors
    sce@images[[1]]@coordinates <- sce@images[[1]]@coordinates / sce@images[[1]]@scale.factors$lowres
    return(sce)
}
