generateRandomSpots <- function(...) {
    set.seed(123456)
    params <- list(...)
    rand.cnt <- rpois(n = 1000, lambda = params[[2]]) + 1
    cell.names <- colnames(params[[1]])
    count.mat <- GetAssayData(params[[1]], slot = "count")
    rand.spots <- lapply(1:length(rand.cnt), function(idx) {
        cnt <- rand.cnt[idx]
        rand.df <- sample(cell.names, cnt) %>%
            count.mat[, .] %>%
            as.data.frame() %>%
            rowSums()
        return(rand.df)
    }) %>%
        do.call(cbind, .) %>%
        as.data.frame() %>%
        `rownames<-`(rownames(count.mat))
    return(list(cnt = rand.cnt, spots = rand.spots))
}

getDensity <- function(x, y, ...) {
    dens <- MASS::kde2d(x, y, ...)
    ix <- findInterval(x, dens$x)
    iy <- findInterval(y, dens$y)
    ii <- cbind(ix, iy)
    return(dens$z[ii])
}
