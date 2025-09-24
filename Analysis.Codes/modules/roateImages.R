rotimat <- function(foo, rotation) {
    if (!is.matrix(foo)) {
        cat("Input is not a matrix")
        return(foo)
    }
    if (!(rotation %in% c("180", "Hf", "Vf", "R90", "L90"))) {
        cat("Rotation should be either L90, R90, 180, Hf or Vf\n")
        return(foo)
    }
    if (rotation == "180") {
        foo <- foo %>%
            .[, dim(.)[2]:1] %>%
            .[dim(.)[1]:1, ]
    }
    if (rotation == "Hf") {
        foo <- foo %>%
            .[, dim(.)[2]:1]
    }

    if (rotation == "Vf") {
        foo <- foo %>%
            .[dim(.)[1]:1, ]
    }
    if (rotation == "L90") {
        foo <- t(foo)
        foo <- foo %>%
            .[dim(.)[1]:1, ]
    }
    if (rotation == "R90") {
        foo <- t(foo)
        foo <- foo %>%
            .[, dim(.)[2]:1]
    }
    return(foo)
}

rotateSeuratImage <- function(seuratVisumObject, slide = "slice1", rotation = "Vf") {
    if (!(rotation %in% c("180", "Hf", "Vf", "L90", "R90"))) {
        cat("Rotation should be either 180, L90, R90, Hf or Vf\n")
        return(NULL)
    } else {
        seurat.visium <- seuratVisumObject
        ori.array <- (seurat.visium@images)[[slide]]@image
        img.dim <- dim(ori.array)[1:2] / (seurat.visium@images)[[slide]]@scale.factors$lowres
        new.mx <- c()
        # transform the image array
        for (rgb_idx in 1:3) {
            each.mx <- ori.array[, , rgb_idx]
            each.mx.trans <- rotimat(each.mx, rotation)
            new.mx <- c(new.mx, list(each.mx.trans))
        }

        # construct new rgb image array
        new.X.dim <- dim(each.mx.trans)[1]
        new.Y.dim <- dim(each.mx.trans)[2]
        new.array <- array(
            c(
                new.mx[[1]],
                new.mx[[2]],
                new.mx[[3]]
            ),
            dim = c(new.X.dim, new.Y.dim, 3)
        )

        # swap old image with new image
        seurat.visium@images[[slide]]@image <- new.array

        ## step4: change the tissue pixel-spot index
        img.index <- (seurat.visium@images)[[slide]]@coordinates

        # swap index
        if (rotation == "Hf") {
            seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2] - img.index$imagecol
        }

        if (rotation == "Vf") {
            seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1] - img.index$imagerow
        }

        if (rotation == "180") {
            seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[1] - img.index$imagerow
            seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[2] - img.index$imagecol
        }

        if (rotation == "L90") {
            seurat.visium@images[[slide]]@coordinates$imagerow <- img.dim[2] - img.index$imagecol
            seurat.visium@images[[slide]]@coordinates$imagecol <- img.index$imagerow
        }

        if (rotation == "R90") {
            seurat.visium@images[[slide]]@coordinates$imagerow <- img.index$imagecol
            seurat.visium@images[[slide]]@coordinates$imagecol <- img.dim[1] - img.index$imagerow
        }

        return(seurat.visium)
    }
}
