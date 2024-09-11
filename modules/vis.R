plotSc2Spatial <- function(sp.obj, pseud.coord, x.lab, y.lab, ct.lab = "CellType", tar.ctypes = NULL, pt.size = 1.0, colors = NULL, shape = 21, percent = NULL) {
    pseud.coord <- rename(pseud.coord, x = x.lab, y = y.lab, CellType = ct.lab)
    img <- sp.obj@images[[names(sp.obj@images)]]@image
    coord.nn <- GetTissueCoordinates(sp.obj)
    if (!is.null(tar.ctypes)) {
        tar.ctypes.sub <- intersect(tar.ctypes, pseud.coord$CellType %>% unique())
        if (length(tar.ctypes.sub) == 0) println("No cell types found in sce data, exiting...", status = "ERROR")
        pseud.coord <- subset(pseud.coord, CellType %in% tar.ctypes.sub)
    }
    tar.ctypes <- pseud.coord$CellType %>% unique()
    if (!is.null(colors) && is.null(names(colors))) {
        colors <- colors[1:length(tar.ctypes)] %>% `names<-`(tar.ctypes)
    } else {
        colors <- scales::hue_pal()(tar.ctypes %>% length()) %>% `names<-`(tar.ctypes)
    }
    colors <- colors[tar.ctypes]

    if (!is.null(percent) && percent > 0 && percent <= 1) {
        pseud.coord <- lapply(pseud.coord$CellType %>% unique(), function(ct) {
            df <- subset(pseud.coord, CellType == ct)
            int.idx <- sample(nrow(df), floor(nrow(df) * percent))
            df[int.idx, ]
        }) %>%
            do.call(rbind, .) %>%
            as.data.frame()
    }
    gp <- ggplot(pseud.coord, aes(x = y, y = x, fill = CellType)) +
        annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
        geom_point(shape = shape, size = pt.size) +
        ggplot2::scale_y_reverse() +
        theme_bw() +
        theme(
            panel.background = element_rect(fill = "black"),
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.title = element_text(hjust = 0.5)
        ) +
        scale_fill_manual(values = colors) +
        guides(fill = guide_legend(override.aes = list(size = 4))) +
        ggtitle(sp.obj@images %>% names()) +
        xlab("") +
        ylab("") #+
    # xlim(0, max(coord.nn[, 2])) + ylim(0, max(coord.nn[, 1]))
    return(gp)
}

projUMAP <- function(proj.res, image, cell.colors, f.ctypes = NULL, split = TRUE, nrow = 4, ncol = 5) {
    if (split) {
        gps <- lapply(f.ctypes, function(ct) {
            ggplot(subset(proj.res, CellType == ct), aes(x = y.new, y = x.new, color = CellType)) +
                annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
                geom_point(size = 1) +
                ggplot2::scale_y_reverse() +
                xlim(76.41996724, 492.237532229) +
                ylim(138.640278405, 540.567997997) +
                theme_bw() +
                # theme_void() +
                theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
                ggtitle(ct) +
                scale_color_manual(values = cell.colors[ct]) #+
            xlim(76.41996724, 492.237532229) + ylim(138.640278405, 540.567997997)
        })
        gp <- ggarrange(plotlist = gps, nrow = nrow, ncol = ncol)
    } else {
        gp <- ggplot(cyto.res, aes(x = y.new, y = x.new, color = CellType)) +
            annotation_raster(img, xmin = 0, xmax = ncol(img), ymin = 0, ymax = -nrow(img)) +
            geom_point(size = 1) +
            ggplot2::scale_y_reverse() +
            theme_bw() +
            theme_void() +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
            scale_color_manual(values = cell.colors)
    }
    return(gp)
}
