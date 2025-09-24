runSeurat <- function(obj.sp, obj.sc) {
    ovp.genes <- intersect(rownames(obj.sc), rownames(obj.sp))
    obj.sc <- obj.sc[ovp.genes, ]
    obj.sc <- NormalizeData(obj.sc)
    obj.sc <- FindVariableFeatures(obj.sc)
    obj.sc <- ScaleData(obj.sc, vars.to.regress = c("S.Score", "G2M.Score", "orig.ident", "percent.mt", "percent.ribo"))

    obj.sp <- obj.sp[ovp.genes, ]
    obj.sp <- NormalizeData(obj.sp)
    obj.sp <- FindVariableFeatures(obj.sp)
    obj.sp <- ScaleData(obj.sp)
    anchors <- FindTransferAnchors(reference = obj.sc, query = obj.sp)

    predictions <- TransferData(anchorset = anchors, refdata = obj.ref$Anno.Level.Fig.1)
    obj.sp <- AddMetaData(object = obj.sp, metadata = predictions)
}

runCellTrek <- function(obj.sp, obj.sc, cell.names = "CellType", ...) {
    brain.traint <- CellTrek::traint(st_data = obj.sp, sc_data = obj.sc, sc_assay = "RNA", cell_names = cell.names)
    brain.celltrek <- CellTrek::celltrek(
        st_sc_int = brain.traint,
        int_assay = "traint",
        sc_data = obj.sc,
        sc_assay = "RNA",
        reduction = "pca",
        intp = TRUE,
        intp_pnt = 5000,
        intp_lin = F,
        nPCs = 30,
        ntree = 1000,
        dist_thresh = 0.55,
        top_spot = 5,
        spot_n = 5,
        repel_r = 20,
        repel_iter = 20,
        keep_model = TRUE,
        ...
    )$celltrek
}

runCytoSpace <- function(obj.sp, obj.sc, cell.names = "CellType", out.dir = "./") {
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
    # cyto.res <- read.table('.tmp/assigned_locations.csv', sep = ',', header = TRUE)
    # 	cyto.res <- getRandomCords(cyto.res)
}

runCARD <- function(obj.sp, obj.sc, cell.names = "CellType", num.cell = 10) {
    sp.counts <- GetAssayData(obj.sp, slot = "count")
    sc.counts <- GetAssayData(obj.sc, slot = "count")
    sc.meta <- obj.sc@meta.data
    sc.meta$cellType <- as.vector(sc.meta@meta.data[, cell.names])
    sp.loc <- GetTissueCoordinates(object = obj.sp@images[[1]]) %>% `colnames<-`(c("x", "y"))
    CARD.obj <- createCARDObject(
        sc_count = sc.counts,
        sc_meta = sc.meta,
        spatial_count = sp.counts,
        spatial_location = sp.loc,
        ct.varname = cell.names,
        ct.select = unique(sc.meta@meta.data[, cell.names]),
        sample.varname = "orig.ident",
        minCountGene = 100,
        minCountSpot = 5
    )
    CARD.obj <- CARD_deconvolution(CARD_object = CARD.obj)
    CARD.obj <- CARD_SCMapping(CARD.obj, numCell = num.cell)
    return(CARD.obj)
}
