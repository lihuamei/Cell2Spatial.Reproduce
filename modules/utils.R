#' downSamplSeurat

#' Down-sampling for Seurat obeject.
#' @param obj Seurat object
#' @param cnt Sample size for each ident, default: 2000.
#' @param seed Randon seed, default: 123.
#' @return Subset of seurat object.
#' @export

downSamplSeurat <- function(obj, cnt = 2000, seed = 123, percent = NULL) {
    set.seed(seed)
    cells <- Idents(obj) %>% table()
    sub.cells <- sapply(names(cells), function(xx) {
        sub.cells <- Idents(obj)[Idents(obj) == xx] %>% names()
        cnt <- ifelse(is.null(percent), cnt, length(sub.cells) * percent)
        if (length(sub.cells) > cnt) sub.cells <- sample(sub.cells, cnt, replace = FALSE)
        return(sub.cells)
    }) %>% unlist(use.names = F)
    subset(obj, cells = sub.cells)
}

#' generateCytoSPACE

generateCytoSPACE <- function(obj.sc, obj.st, ctype, prefix, out.dir) {
    mat.sc <- GetAssayData(obj.sc, slot = "count") %>% as.data.frame()
    c.type <- cbind.data.frame(`cell IDs` = obj.sc@meta.data %>% rownames(), CellType = obj.sc@meta.data[, ctype])
    coord.df <- GetTissueCoordinates(obj.st) %>%
        cbind.data.frame(SpotID = colnames(obj.st), .) %>%
        `colnames<-`(c("SpotID", "row", "col"))
    mat.st <- GetAssayData(obj.st, slot = "count") %>% as.data.frame()
    mat.st <- cbind.data.frame(GENES = rownames(mat.st), mat.st)
    mat.sc <- cbind.data.frame(GENES = rownames(mat.sc), mat.sc)

    sc.file <- file.path(out.dir, paste0(prefix, "_sc.txt"))
    data.table::fwrite(mat.sc, file = sc.file, sep = "\t", row.names = FALSE, col.names = TRUE)

    st.file <- file.path(out.dir, paste0(prefix, "_st.txt"))
    data.table::fwrite(mat.st, file = st.file, sep = "\t", row.names = FALSE, col.names = TRUE)

    ctype.file <- file.path(out.dir, paste0(prefix, "_ctype.txt"))
    data.table::fwrite(c.type, file = ctype.file, sep = "\t", row.names = FALSE, col.names = TRUE)

    coord.file <- file.path(out.dir, paste0(prefix, "_coord.txt"))
    data.table::fwrite(coord.df, file = coord.file, sep = "\t", row.names = FALSE, col.names = TRUE)
    return(list(sc = sc.file, st = st.file, ct = ctype.file, pos = coord.file))
}

#' Generate randon position

getRandomCords <- function(cyto.res, seed = 123456, n.workers = 4) {
    cyto.res.sub <- cyto.res[!duplicated(cyto.res$SpotID), ]
    num.cells <- table(cyto.res$SpotID)[cyto.res.sub[, "SpotID"]] %>% as.list()
    sp.coord <- cyto.res.sub[, c("row", "col")] %>% `rownames<-`(cyto.res.sub$SpotID)
    set.seed(seed)
    future::plan("multicore", workers = n.workers)
    if (is.numeric(num.cells)) num.cells <- rep(num.cells, nrow(sp.coord)) %>% as.list()
    sp.ED <- fields::rdist(as.matrix(sp.coord))
    min.dist <- apply(sp.ED, 1, function(xx) order(xx)[2] %>% xx[.]) %>%
        {
            median(.) / 2
        }
    sp.coord.new <- future.apply::future_lapply(1:nrow(sp.coord), function(idx) {
        circle <- spatstat.random::runifdisc(
            num.cells[[idx]],
            radius = min.dist,
            centre = c(sp.coord[idx, 1], sp.coord[idx, 2]),
            nsim = 1,
            drop = TRUE
        )
        tmp.cord <- data.frame(x = circle$x, y = circle$y) %>% `colnames<-`(c("x", "y"))
        tmp.cord$centerSPOT <- paste0(sp.coord[idx, 1], "x", sp.coord[idx, 2])
        tmp.cord$centerX <- sp.coord[idx, 1]
        tmp.cord$centerY <- sp.coord[idx, 2]
        tmp.cord$SpotName <- rownames(sp.coord)[idx]
        return(tmp.cord)
    }, future.seed = TRUE) %>%
        do.call(rbind, .) %>%
        as.data.frame()
    sp.coord.new <- sp.coord.new[!duplicated(paste0(sp.coord.new$x, "x", sp.coord.new$y)), ] %>%
        `rownames<-`(paste0(sp.coord.new$x, "x", sp.coord.new$y))

    cyto.res$x.new <- 0
    cyto.res$y.new <- 0
    for (sp in unique(cyto.res$SpotID)) {
        idxes.1 <- which(cyto.res$SpotID == sp)
        idxes.2 <- which(sp.coord.new$SpotName == sp)
        cyto.res[idxes.1, "x.new"] <- sp.coord.new[idxes.2, "x"]
        cyto.res[idxes.1, "y.new"] <- sp.coord.new[idxes.2, "y"]
    }
    return(cyto.res)
}

reformatSimCoord <- function(sim.res, obj.sp) {
    coord.df <- GetTissueCoordinates(obj.sp) %>% `colnames<-`(c("row", "col"))
    sim.df <- lapply(names(sim.res), function(xx) {
        tmp.xx <- sim.res[[xx]]
        cbind.data.frame(CellName = tmp.xx, SpotID = xx)
    }) %>%
        do.call(rbind, .) %>%
        as.data.frame()

    sim.df <- cbind.data.frame(sim.df, coord.df[sim.df$SpotID, ])
    sim.df <- getRandomCords(sim.df)
    return(sim.df)
}

covertSeurat <- function(proj.res, obj.st, obj.sc) {
    count.sc <- GetAssayData(obj.sc, slot = "count")
    count.CT <- as(count.sc[, proj.res$OriginalCID], "sparseMatrix")
    colnames(count.CT) <- rownames(proj.res)

    sce <- CreateSeuratObject(count = count.CT, meta.data = proj.res, project = "Cell2Spatial", assay = "Spatial")
    sce@images <- obj.st@images
    sce@images[[1]]@assay <- DefaultAssay(sce)
    sce@images[[1]]@coordinates <- data.frame(imagerow = proj.res$x.new, imagecol = proj.res$y.new) %>% `rownames<-`(rownames(proj.res))
    sce@images[[1]]@scale.factors <- obj.st@images[[1]]@scale.factors
    return(sce)
}

createMonocleObject <- function(obj.seu) {
    sce <- obj.seu
    sample.ann <- sce@meta.data
    gene.ann <- data.frame(
        gene_short_name = rownames(sce@assays$RNA),
        row.names = rownames(sce@assays$RNA)
    )
    pd <- new("AnnotatedDataFrame", data = sample.ann)
    fd <- new("AnnotatedDataFrame", data = gene.ann)
    ct <- as.data.frame(sce@assays$RNA@counts)
    sc.cds <- newCellDataSet(
        as.matrix(ct),
        phenoData = pd,
        featureData = fd,
        expressionFamily = negbinomial.size(),
        lowerDetectionLimit = 1
    )
    return(sc.cds)
}

monocle2Analysis <- function(cds, cut.expr = 0.1) {
    cds <- estimateSizeFactors(cds)
    cds <- estimateDispersions(cds)
    cds <- detectGenes(cds, min_expr = cut.expr)

    disp_table <- dispersionTable(cds)
    unsup_clustering_genes <- subset(disp_table, mean_expression >= cut.expr)
    cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)
    cds <- reduceDimension(cds, method = "DDRTree")
    cds <- orderCells(cds)
    return(cds)
    # plot_cell_trajectory(cds, color_by = "Pseudotime")
}
