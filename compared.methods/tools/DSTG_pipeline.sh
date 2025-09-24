set -ex

scrna_path=$1
spatial_path=$2
celltype_key=$3
output_path=$4
prefix=$5

echo $scrna_path
echo $spatial_path
echo $celltype_key

BACK_DIR=pwd/output_path
cd tools/DSTG/DSTG

#Rscript convert_data_h5ad.r $scrna_path $celltype_key $spatial_path
Rscript convert_data_RDS.R ../../../$scrna_path ../../../$spatial_path $celltype_key $prefix

Rscript convert_data.R scRNAseq_data_$prefix.RDS  spatial_data_$prefix.RDS  scRNAseq_label_$prefix.RDS 

python  train.py

cp ./DSTG_Result/predict_output.csv $BACK_DIR/DSTG_output.csv

