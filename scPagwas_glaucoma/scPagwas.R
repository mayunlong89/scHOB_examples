library(scPagwas)
library(readr)
library(dplyr)
library(Seurat)
library(tidyverse)
# library(rhdf5)
library(ggplot2)
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(reshape2)
library(ggdendro)
library(grImport2)
library(gridExtra)
library(grid)
library(sisal)
require(RColorBrewer)
require(ggsci)
require(ggpubr)
source(system.file("extdata", "plot_scpathway_contri_dot.R", package = "scPagwas"))
library(optparse)

option_list <- list(
  make_option(c("-s", "--seuDir"), type = "character", default = FALSE,
              help = "data destination"
  ),
  make_option(c("-g", "--gwasDes"), type = "character", default = FALSE,
              help = "gwas destinatation"
  ), 
  make_option(c("-d", "--dir"), type = "character", default = FALSE,
              help = 'the work directory'
  ), 
  make_option(c("-t", "--trait"), type = "character", default = FALSE,
              help = 'trait'
  )
)

opt <- parse_args(OptionParser(option_list = option_list))
work_dir <- opt$dir
trait <- opt$trait
seuDir <- opt$seuDir
gwas <- opt$gwasDes
print(trait)
print(seuDir)
print(gwas)

result_save_path <- paste0('./results/scPagwas_', trait, '/')
figure_path <- paste0(result_save_path, '/figures/')
if(!dir.exists(result_save_path)){
  dir.create(result_save_path, recursive = TRUE)
}
if(!dir.exists(figure_path)){
  dir.create(figure_path, recursive = TRUE)
}


setwd(work_dir)
seu <- readRDS(seuDir)
Idents(seu) <- seu@meta.data$annotation2
table(Idents(seu))


Pagwas <- scPagwas_main(Pagwas = NULL,
                        gwas_data = gwas,
                        Single_data = seu,
                        output.prefix = trait,
                        output.dirs = result_save_path,
                        singlecell = T,
                        celltype = T,
                        Pathway_list = Genes_by_pathway_kegg,
                        assay = "RNA",
                        block_annotation = block_annotation,
                        seurat_return = T,
                        # n.cores = 15,
                        chrom_ld = chrom_ld
)

saveRDS(Pagwas,file = paste0(result_save_path, trait, 'scPagwas.rds'))
color26 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8",
                      "#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186",
                      "#558776","#E99497","#FFBD9B","#0A1D37")
                      
pdf(file = paste0(figure_path, 'bootstrap_p_Plot.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
Bootstrap_P_Barplot(p_results=Pagwas@misc$bootstrap_results$bp_value[-1],
                    p_names=rownames(Pagwas@misc$bootstrap_results)[-1],
                    figurenames = "Bootstrap_P_Barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = paste0(trait, '_celltype')
)
dev.off()
pdf(file = paste0(figure_path, 'plot_bar_positie_nagtive1.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="celltype",
                         var_group="positiveCells",
                         vec_group_colors=c("#E8D0B3","#7EB5A6"),
                         do_plot = T)
dev.off()
pdf(file = paste0(figure_path, 'plot_bar_positie_nagtive2.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="positiveCells",
                         var_group="celltype",
                         p_thre = 0.01,
                         vec_group_colors=NULL,
                         f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                         do_plot = T)
dev.off()

ct <-  as.vector(unique(Pagwas@meta.data$celltype))
pdf(file = paste0(figure_path, 'plot_scpathway_dot.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
plot_scpathway_dot(Pagwas=Pagwas,
                   celltypes=ct,
                   topn_path_celltype=5,
                   filter_p=0.05,
                   max_logp=15,
                   display_max_sizes=F,
                   size_var ="logrankPvalue" ,
                   col_var="proportion",
                   shape.scale = 8,
                   cols.use=c("lightgrey", "#E45826"),
                   dend_x_var = "logrankPvalue",
                   dist_method="euclidean",
                   hclust_method="ward.D",
                   do_plot = T,
                   #figurenames = "Pathway_plot.pdf",
                   width = 7,
                   height = 7)
dev.off()

pdf(file = paste0(figure_path, 'heritability_cor_scatterplot.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
heritability_cor_scatterplot(gene_heri_cor=Pagwas@misc$gene_heritability_correlation,
                             topn_genes_label=10,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width = 7,
                             height = 7)
dev.off()

top5genes<-rownames(Pagwas@misc$gene_heritability_correlation)[order(Pagwas@misc$gene_heritability_correlation,decreasing = T)[1:5]]
pdf(file = paste0(figure_path, 'plot_vln_Corgenes.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
plot_vln_Corgenes(seurat_obj=Pagwas,
                  assay="RNA", 
                  slot="data",
                  var_group="celltype",
                  vec_features=top5genes,
                  vec_group_colors= color26,
                  do_plot = T
)
dev.off()

pdf(file = paste0(figure_path, 'plot_gpas_trs_pval.pdf'), width = unit(10, 'inch'), height = unit(10, 'inch'))
scPagwas_Visualization(Single_data=Pagwas,
                       p_thre = 0.05,
                       FigureType = "umap",
                       width = 7,
                       height = 7,
                       lowColor = "white", 
                       highColor = "red",
                       output.dirs="figure",
                       size = 0.5,
                       do_plot = T)
dev.off()








