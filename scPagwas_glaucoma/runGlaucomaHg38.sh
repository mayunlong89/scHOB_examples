#!/bin/sh
#PBS -N runGlaucomaHg38
#PBS -q fat
#PBS -l nodes=fat03
#PBS -l ncpus=10
#PBS -l mem=300gb

source activate scRNA

shDir='/share2/pub/zhuzb/zhuzb/CodeResp/scPagwas/'
dir='/share2/pub/zhuzb/zhuzb/AMDScpagwas/'
seuDir='/share2/pub/zhuzb/zhuzb/AMDScpagwas/seurat_obj_final.rds'
gwasDes='/share2/pub/zhuzb/zhuzb/AMDScpagwas/data/glaucoma_scPagwasHg38.csv'
trait='glaucoma'

sh $shDir/lib/runscPagwas.sh --dir $dir --seuDir $seuDir --gwasDes $gwasDes --trait $trait
