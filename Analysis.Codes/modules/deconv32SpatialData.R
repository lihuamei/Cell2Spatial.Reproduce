library(devtools)
library(Seurat)
load_all("../Github/Cell2Spatial/")

obj.st <- readRDS('../0.data/SimualtedSpatalData/st.obj.lst.RDS')
obj.sc <- readRDS('../0.data/SimualtedSpatalData/sc.obj.lst.RDS')

out.dir <- './prop.sct'
dir.create(out.dir, showWarnings = FALSE)
data.names <- names(obj.st)
prop.lst <- parallel::mclapply(data.names, function(dta) {
    message(dta)		   
    sc.obj <- obj.sc[[dta]]
    sp.obj <- obj.st[[dta]]
    Idents(sc.obj) <- sc.obj$celltype_final
    est.pop <- runCell2Spatial(sp.obj, sc.obj, knn.spots = 0, n.workers = 1, normalize.method = "LogNormalize", adjust.deconv = "NONE")
    saveRDS(est.pop, file = file.path(out.dir, sprintf('%s.RDS', dta)) )
}, mc.cores = 4)

saveRDS(prop.lst, file = 'prop.lst.RDS')
