import squidpy as sq
import skimage.io
import numpy as np

adata = sq.read.visium('./', counts_file = 'filtered_feature_bc_matrix.h5')
lib = next(iter(adata.uns["spatial"]))                      
sp = adata.uns["spatial"][lib]
hires = sp["images"]["hires"]                                   

img = sq.im.ImageContainer(np.asarray(hires), layer="image", scale=sp["scalefactors"]["tissue_hires_scalef"])

sq.im.process(
    img=img,
    layer="image",
    method="smooth",
    layer_added="image_smooth",
)

sq.im.segment(
    img=img,
    layer="image",
    method="watershed",
    channel=2,
    layer_added="segmented_watershed",
)

features_kwargs = {"segmentation": {"label_layer": "segmented_watershed"}}

sq.im.calculate_image_features(
    adata=adata,
    img=img,
    layer="image",        
    features="segmentation",
    features_kwargs=features_kwargs,
    key_added="features_segmentation",
    library_id=lib,                      
    n_jobs=1,
)

adata.obsm['features_segmentation'].to_csv('count_res.xls', index = True, header = True, sep = '\t')
