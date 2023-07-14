library(SingleCellExperiment)
library(scater)
library(scran)
library(AnnotationHub)
library(biomaRt)
library(pheatmap)
library(data.table)
library(gplots)
library(gridExtra)
library(latex2exp)
library(tidyr)
library(AnnotationDbi)
library(scry)
library(scDblFinder)
library(ggVennDiagram)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(ggpubr)
library(PCAtools)
library(intrinsicDimension)
library(bluster)
library(Seurat)
library(SC3)
library(cluster)
library(dplyr)
library(viridisLite)
library(viridis)
library(sctransform)
library(glmGamPoi)
library(parallel)

# preparetion ##################################################################
# set working directory
setwd('/home/s2321661/Test')

# set random seed
seed <- 1000



# load data ####################################################################
# read scRNA seq data (raw count) and cell annotation
raw_count <- fread('counts.txt', skip=1, header=TRUE, data.table = FALSE)
rownames(raw_count) <- raw_count$Geneid
raw_count <- raw_count[,7:ncol(raw_count)]
raw_count <- as.matrix(raw_count)
anno <- read.csv('annotation.csv')
idx <- match(colnames(raw_count),anno$Dir)
colnames(raw_count) <- anno$ArrayExpress_Accession[idx]
anno <- anno[idx,]
coldata <- anno[,c('Cell.Type','Source.Name')]
rownames(coldata) <- anno$ArrayExpress_Accession
dim(raw_count)
sce <- SingleCellExperiment(assay=list(counts = raw_count),colData=coldata)
sce
sce_origin <- sce # backup
rm(idx)
rm(anno)
rm(coldata)
rm(raw_count)

# gene name annotation
ah <- AnnotationHub()
#display(ah)  # search 'GRCm39' and 'EnsDb' to find v109 ID
ens.mm.v109 <- AnnotationHub()[["AH109655"]]
columns(ens.mm.v109)
gene_anno <- AnnotationDbi::select(ens.mm.v109, keys=rownames(sce), 
                                   keytype='GENEID', column='GENENAME')
head(gene_anno)
sum(gene_anno$GENENAME == '')  # the number of failures in gene name annotation
rowData(sce)$GENENAME <- make.names(gene_anno$GENENAME[match(rownames(sce),gene_anno$GENEID)],
                                    unique = TRUE)
sum(grepl('^NA',rowData(sce)$GENENAME))  
sum(grepl('^X\\.|^X$',rowData(sce)$GENENAME) == 1)  # the number of failures in gene name annotation
rowData(sce)$ENSEMBL <- rownames(sce)

# ERCC name annotation
is.spike <- grepl("^ERCC", rowData(sce)$ENSEMBL)
count(is.spike)
idx <-  which(is.spike)
rowData(sce)$GENENAME[idx] <- rowData(sce)$ENSEMBL[idx]
rowData(sce)$GENENAME[idx]

# change sce rownames from ENSEMBL IDs to gene names and ERCC names
rownames(sce) <- rowData(sce)$GENENAME
sce
head(rowData(sce))
length(grep('^NA',rownames(sce)))
sum(colSums(assays(sce[is.spike, ])$counts) == 0)  # check how many cells do not have spike-ins
sce_origin <- sce  #backup
rm(idx)
rm(ah)
rm(gene_anno)
rm(ens.mm.v109)
save(list = ls(), file = 'sce_preperation.RData')
#load("sce_preperation.RData")

# Get Mitochondrial genes' ENSEMBL ID from Biomart
#listEnsembl()
#ensembl <- useEnsembl(biomart = "genes")
#head(listDatasets(ensembl))
#searchDatasets(mart = ensembl, pattern = "GRCm39")
#ensembl <- useDataset(dataset='mmusculus_gene_ensembl', mart=ensembl)
#head(listFilters(ensembl))
#head(listAttributes(ensembl))
#MT_genes <- getBM(filters = 'chromosome_name', values='MT',
#                  attributes = c('ensembl_gene_id','external_gene_name'), 
#                  mart = ensembl)
#head(MT_genes)
#MT_id <- paste(MT_genes$ensembl_gene_id, collapse='|')
#head(MT_id)
#is.mito <- grepl(MT_id, rowData(sce)$ENSEMBL)
#count(is.mito)
#rm(ensembl)
#rm(MT_id)

# get mitochondrial genes from annotation files got from Biomart Website
MT_genes <- read.csv("MT_genes.csv")
head(MT_genes)
MT_id <- paste(MT_genes$Gene.stable.ID, collapse='|')
head(MT_id)
is.mito <- grepl(MT_id, rowData(sce)$ENSEMBL)
count(is.mito)
rm(MT_id)



# # Doublet removal###############################################################
# set.seed(seed) # different seed varies a bit (can vary more than 10 cells...)
# dbl.out <- capture.output(sce_dbl <- scDblFinder(sce) , type = "message")
# thres <- as.numeric(gsub("[^0-9.]", "", dbl.out[grepl('Threshold', dbl.out)]))
# thres
# table(sce_dbl$scDblFinder.class)
# 
# # histogram of doublet scores
# ggplot(data.frame(class=sce_dbl$scDblFinder.class, score=sce_dbl$scDblFinder.score)) +
#   geom_histogram(aes(x=score),bins = 100, color='grey80') +
#   geom_vline(aes(xintercept=thres), color = 'red', linetype="dashed") +
#   xlab('Doublet score')
# ggsave('figures/doublets_hist.jpg',device='jpg', width = 8, height = 5)
# 
# # usually doublets show higher levels of expression
# stats <- perCellQCMetrics(sce_dbl)
# sce_dbl$total_counts <- stats$sum/1e6
# plotColData(sce_dbl, x="Cell.Type", y="total_counts", colour_by="scDblFinder.class") +
#   ggtitle("Library sizes") + ylab('Total Reads Count') + xlab('Cell Type') +
#   theme(plot.title = element_text(size=12),
#         axis.text.x = element_text(size = 12),
#         axis.text.y = element_text(size = 12),
#         axis.title.x=element_text(size = 12),
#         axis.title.y=element_text(size = 12),
#         legend.text=element_text(size=10)) +
#   guides(color=guide_legend("scDblFinder", title.theme = element_text(size = 12))) +
#   scale_y_continuous(labels = c(0,
#                                 sapply(sort(seq(5,max(sce_dbl$total_counts),5)),
#                                        function(x){latex2exp::TeX(paste0('$',x,'\\times 10^6$'))})),
#                      breaks = sort(c(0,seq(5,max(sce_dbl$total_counts),5))))
# ggsave('figures/doublets_lib_size_scatter.pdf',device='pdf', width = 8, height = 8)
# 
# # usually doublets would appear at the periphery of the singlet cluster
# sce_dbl <- computeSumFactors(sce_dbl, cluster = quickCluster(sce_dbl))
# sce_dbl <- logNormCounts(sce_dbl)
# var.dbl <- modelGeneVar(sce_dbl)
# hvg.dbl <- getTopHVGs(var.dbl,n=2000)
# set.seed(seed)
# sce_dbl <- runPCA(sce_dbl, subset_row = hvg.dbl)
# set.seed(seed)
# sce_dbl <- runTSNE(sce_dbl, perplexity=20, dimred="PCA")
# plotTSNE(sce_dbl, colour_by = 'scDblFinder.class', shape_by = 'Cell.Type')
# ggsave('figures/doublets_tsne.jpg', device='jpg', width = 6, height = 5)
# 
# # remove doublets
# sce[,sce_dbl$scDblFinder.class != 'doublet']
# rm(sce_dbl)
# rm(dbl.out)
# rm(thres)
# rm(stats)
# rm(var.dbl)
# rm(hvg.dbl)



# cell QC ######################################################################
QC_stats.cell <- perCellQCMetrics(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
# or we can store the stats directly to colData 
# by addPerCellQC(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
QC_stats.cell
sce$total_counts <- QC_stats.cell$sum/1e6
sce$detected_genes <- QC_stats.cell$detected
sce$mito_percent <- QC_stats.cell$subsets_Mt_percent
sce$ERCC_percent <- QC_stats.cell$subsets_ERCC_percent

# check QC failed cells (filter cells)
cell_filter <- perCellQCFilters(QC_stats.cell,sub.fields=c("subsets_Mt_percent",
                                                      "subsets_ERCC_percent"))
# or use quickPerCellQC(QC_stats.cell,sub.fields=c("subsets_Mt_percent", "subsets_ERCC_percent"))
cell_filter  
colSums(as.matrix(cell_filter))  # check how many cells were dropped 
summary(cell_filter$discard)  # check total discarded cells count
sce$discard <- cell_filter$discard
sce$low_lib_size <- cell_filter$low_lib_size
sce$low_detected_genes <- cell_filter$low_n_features
sce$high_mito_percent <- cell_filter$high_subsets_Mt_percent
sce$high_ERCC_percent <- cell_filter$high_subsets_ERCC_percent
sce_origin <- sce  # backup
sce

# Scatter plots that summarize cell QC
thres_total_counts <- attributes(cell_filter$low_lib_size)$thresholds[1]/1e6
thres_detected_genes <- attributes(cell_filter$low_n_features)$thresholds[1]
thres_mito_percent <- attributes(cell_filter$high_subsets_Mt_percent)$thresholds[2]
thres_ERCC_percent <- attributes(cell_filter$high_subsets_ERCC_percent)$thresholds[2]
QC_scatter <- gridExtra::grid.arrange(
  plotColData(sce, x="Cell.Type", y="total_counts", colour_by="discard") +
    ggtitle("Library sizes") + ylab('Total Reads Count') + xlab('Cell Type') + 
    geom_hline(yintercept = thres_total_counts,
               linetype = "dashed", color = "red") + 
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(labels = c(0,
                                  sapply(sort(seq(5,max(sce$total_counts),5)), 
                                         function(x){latex2exp::TeX(paste0('$',x,'\\times 10^6$'))})),
                       breaks = sort(c(0,seq(5,max(sce$total_counts),5)))) +
    annotate('text', x = 3.6, y = thres_total_counts,
              label = latex2exp::TeX(paste0('$',round(thres_total_counts,2),'\\times 10^6$')), 
              vjust = 1.2, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce, x="Cell.Type", y="detected_genes", colour_by="discard")+
    ggtitle("Detected genes") + ylab('Number of detected genes') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_detected_genes,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce$detected_genes),2000)) +
    annotate('text',x = 3.5, y = thres_detected_genes,
              label = round(thres_detected_genes,2), 
              vjust = 1.7, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce, x="Cell.Type", y="mito_percent", colour_by="discard")+
    ggtitle("Mito percent") + ylab('Mitochondrial proportion (%)') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_mito_percent,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce$mito_percent),10)) +
    annotate('text',x = 3.5, y = thres_mito_percent,
              label = round(thres_mito_percent,2), 
              vjust = -0.8, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce, x="Cell.Type", y="ERCC_percent", colour_by="discard")+
    ggtitle("ERCC percent") + ylab('ERCC proportion (%)') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_ERCC_percent,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce$ERCC_percent),10)) +
    annotate('text',x = 3.5, y = thres_ERCC_percent,
              label = round(thres_ERCC_percent,2), 
              vjust = -0.8, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  ncol=2,
  nrow=2
)
ggsave('figures/QC_summary.pdf',QC_scatter, device='pdf', width = 15, height = 15)

# Table summary of cell QC
total.cells <- table(sce$Cell.Type)
total.cells <- c(total.cells,Total=sum(total.cells))
# lib size
low_lib_size.cells <- table(colData(sce)[sce$low_lib_size,]$Cell.Type)
low_lib_size.cells <- c(low_lib_size.cells, Total=sum(low_lib_size.cells))
low_lib_size.prop <- round(low_lib_size.cells/total.cells,4)*100
# detected genes
low_detected_genes.cells <- table(colData(sce)[sce$low_detected_genes,]$Cell.Type)
low_detected_genes.cells <- c(low_detected_genes.cells, Total=sum(low_detected_genes.cells))
low_detected_genes.prop <- round(low_detected_genes.cells/total.cells,4)*100
# mito percent
high_mito_percent.cells <- table(colData(sce)[sce$high_mito_percent,]$Cell.Type)
high_mito_percent.cells <- c(high_mito_percent.cells, Total=sum(high_mito_percent.cells))
high_mito_percent.prop <- round(high_mito_percent.cells/total.cells,4)*100
# ERCC percent
high_ERCC_percent.cells <- table(colData(sce)[sce$high_ERCC_percent,]$Cell.Type)
high_ERCC_percent.cells <- c(high_ERCC_percent.cells, Total=sum(high_ERCC_percent.cells))
high_ERCC_percent.prop <- round(high_ERCC_percent.cells/total.cells,4)*100
# discard
discard.cells <- table(colData(sce)[sce$discard,]$Cell.Type)
discard.cells <- c(discard.cells, Total=sum(discard.cells))
discard.prop <- round(discard.cells/total.cells,4)*100
# keep
keep.cells <- table(colData(sce)[!sce$discard,]$Cell.Type)
keep.cells <- c(keep.cells, Total=sum(keep.cells))
keep.prop <- round(keep.cells/total.cells,4)*100
# create table
QC_table <- cbind(total.cells, 
                  low_lib_size.cells,
                  low_lib_size.prop,
                  low_detected_genes.cells,
                  low_detected_genes.prop,
                  high_mito_percent.cells,
                  high_mito_percent.prop,
                  high_ERCC_percent.cells,
                  high_ERCC_percent.prop,
                  discard.cells, 
                  discard.prop,
                  keep.cells,
                  keep.prop)
QC_table <- as.data.frame(QC_table)
QC_table
# lib size
QC_table['low_lib_size.prop'] <- sapply(QC_table['low_lib_size.prop'],
                                              function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, low_lib_size, 
                         c(low_lib_size.cells, low_lib_size.prop), sep=' ')
# detected genes
QC_table['low_detected_genes.prop'] <- sapply(QC_table['low_detected_genes.prop'],
                                              function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, low_detected_genes, 
             c(low_detected_genes.cells, low_detected_genes.prop), sep=' ')
# mito percent
QC_table['high_mito_percent.prop'] <- sapply(QC_table['high_mito_percent.prop'],
                                              function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, high_mito_percent, 
                         c(high_mito_percent.cells, high_mito_percent.prop), sep=' ')
# ERCC percent
QC_table['high_ERCC_percent.prop'] <- sapply(QC_table['high_ERCC_percent.prop'],
                                              function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, high_ERCC_percent, 
                         c(high_ERCC_percent.cells, high_ERCC_percent.prop), sep=' ')
# discard
QC_table['discard.prop'] <- sapply(QC_table['discard.prop'],
                                             function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, discard, 
                         c(discard.cells, discard.prop), sep=' ')
# keep 
QC_table['keep.prop'] <- sapply(QC_table['keep.prop'],
                                   function(x){paste0('(',x,'%)')})
QC_table <- tidyr::unite(QC_table, keep, 
                         c(keep.cells, keep.prop), sep=' ')
# add index column
QC_table <- cbind(rownames(QC_table),QC_table)
# change column names
colnames(QC_table) <- c('Cell Types', 
                        'Total', 
                        'Low Library Size',
                        'Low Detected Genes',
                        'High Mito Percent',
                        'High ERCC Percent',
                        'Total Discard', 
                        'Keep')
QC_table
write.csv(QC_table ,'figures/QC_table.csv', row.names = FALSE)
rm(total.cells)
rm(low_lib_size.cells)
rm(low_lib_size.prop)
rm(low_detected_genes.cells)
rm(low_detected_genes.prop)
rm(high_mito_percent.cells)
rm(high_mito_percent.prop)
rm(high_ERCC_percent.cells)
rm(high_ERCC_percent.prop)
rm(discard.cells)
rm(discard.prop)
rm(keep.cells)
rm(keep.prop)

# Venn Diagram summary of failures in cell QC
venn <- as.data.frame(cell_filter)
QC_venn <- ggVennDiagram(apply(venn[1:4], 2, function(x) which(x == TRUE)),
                         label_alpha=0, 
                         set_color = c("deepskyblue2","darkolivegreen","darkorange","darkorchid"), 
                         category.names = c(paste("Low Library Size (",sum(venn[1]),")",sep = ''),
                                            paste("Low Detected Genes (",sum(venn[2]),")",sep = ''),
                                            paste("High Mito Percent (",sum(venn[3]),")",sep = ''),
                                            paste("High ERCC percent (",sum(venn[4]),")",sep = ''))) + 
  scale_fill_gradient(low="white",high = "coral1") +
  ggtitle('Summary of failures in cell QC') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_color_manual(values = c("deepskyblue2","darkolivegreen","darkorange","darkorchid")) +
  annotate('table',label = QC_table[,c(1,2,7,8)],x = 0, y = 0,vjust = 0.1, hjust = -0.4) +
  theme(legend.position = c(0.95, 0.5))
QC_venn
ggsave('figures/QC_venn.jpg',QC_venn, device='jpg', width = 8, height = 8)
rm(venn)

# Histogram summary of cell QC
QC_hist <- gridExtra::grid.arrange(
  ggplot(as.data.frame(QC_stats.cell), aes(x=sum/1e6))+
    geom_histogram(color='grey80',bins=50) +
    geom_density(alpha=.2, fill="blue") +
    xlab('Library sizes (millions)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_total_counts), color = 'red', linetype="dashed") +
    annotate('text',x=thres_total_counts-2,y=100, label=round(thres_total_counts,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(QC_stats.cell), aes(x=detected))+
    geom_histogram(color='grey80',bins=50) +
    xlab('Number of detected genes') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_detected_genes), color = 'red', linetype="dashed") +
    annotate('text',x=thres_detected_genes-1400,y=60, label=round(thres_detected_genes,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(QC_stats.cell), aes(x=subsets_Mt_percent))+
    geom_histogram(color='grey80',bins=50) +
    xlab('Mitochondrial proportion (%)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_mito_percent), color = 'red', linetype="dashed") +
    annotate('text',x=thres_mito_percent+6,y=200, label=round(thres_mito_percent,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(QC_stats.cell), aes(x=subsets_ERCC_percent))+
    geom_histogram(color='grey80',bins=50) +
    xlab('ERCC proportion (%)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_ERCC_percent), color = 'red', linetype="dashed") +
    annotate('text',x=thres_ERCC_percent+5,y=300, label=round(thres_ERCC_percent,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ncol=2
)
ggsave('figures/QC_hist.jpg',QC_hist, device='jpg', width = 10, height = 8)

# drop cells 
# cell QC fail
sce <- sce[, !cell_filter$discard]

# remove cells with no expression of Oct4 and Sox2 (cell QC 2)
# check how many cells are not expressing Oct4 in each cell type
Oct4_fail <- c(
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == '2i'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'a2i'] == 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'serum'] == 0)
)
Oct4_fail
# check how many cells are not expressing Sox2 in each cell type
Sox2_fail <- c(
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == '2i'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'a2i'] == 0),
  sum(assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'serum'] == 0)
)
Sox2_fail
# summary
Cell_Types <- c('2i','a2i','serum')
Keep <- c(
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == '2i'] != 0 &
        assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == '2i'] != 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'a2i'] != 0 &
        assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'a2i'] != 0),
  sum(assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),colData(sce)$Cell.Type == 'serum'] != 0 &
        assays(sce)$counts[grepl('^Sox2$',rownames(sce)),colData(sce)$Cell.Type == 'serum'] != 0)
)
Keep
Total <- c(
  sum(colData(sce)$Cell.Type == '2i'),
  sum(colData(sce)$Cell.Type == 'a2i'),
  sum(colData(sce)$Cell.Type == 'serum')
)
QC_table_Oct4_Sox2 <- data.frame(Cell_Types, Total, Oct4_fail, Sox2_fail, Keep)
QC_table_Oct4_Sox2
write.csv(QC_table_Oct4_Sox2 ,'figures/QC_table_Oct4_Sox2.csv', row.names = FALSE)
rm(Oct4_fail)
rm(Sox2_fail)
rm(Cell_Types)
rm(Keep)
rm(Total)

# drop cells
sce <- sce[,assays(sce)$counts[grepl('^Pou5f1$',rownames(sce)),] != 0]
sce <- sce[,assays(sce)$counts[grepl('^Sox2$',rownames(sce)),] != 0]
sce



# gene filtering ###############################################################
# remove spikes
sce <- splitAltExps(sce, is.spike)
altExpNames(sce) <- 'spikes'
sce

# remove low expression genes (express in less than 3 cells)
QC_stats.gene <- perFeatureQCMetrics(sce)
detected_cell_prop.hist <- ggplot(as.data.frame(QC_stats.gene), aes(x=detected)) +
  geom_histogram(color='grey80',bins=100) +
  xlab('Detected %') +
  ylab('Counts') + 
  ggtitle('Genes with expression >0')
detected_cell_prop.hist
ggsave('figures/detected_cell_prop.jpg',detected_cell_prop.hist, 
       device='jpg', width = 8, height = 6)
rowData(sce)$mean <- QC_stats.gene$mean
rowData(sce)$detected_prop <- QC_stats.gene$detected
sum(rowSums(assays(sce)$counts > 0) >= 3)  # check how many genes are left 
# remove genes
sce <- sce[rowSums(assays(sce)$counts > 0) >= 3,]
sce



# Normalisation ################################################################
sce <- computeSumFactors(sce, cluster = quickCluster(sce))
sce <- logNormCounts(sce)
sce  # assays(sce) has 'logcounts'



# check oct4 variation #########################################################
# raw counts
oct4_variation_rawcounts <- plotExpression(sce,'Pou5f1',x='Cell.Type',exprs_values = "counts") +
  xlab('Cell Type') + ylab('Raw Counts') + 
  ggtitle('Oct4 expression level (raw counts)')
oct4_variation_rawcounts
ggsave('figures/oct4_variation_rawcounts.jpg',oct4_variation_rawcounts, device='jpg', width = 8, height = 10)

# log counts
oct4_variation_logcounts <- plotExpression(sce,'Pou5f1',x='Cell.Type',exprs_values = "logcounts") +
  xlab('Cell Type') + ylab('Log Counts') + 
  ggtitle('Oct4 expression level after Normalisation (log counts)')
oct4_variation_logcounts
ggsave('figures/oct4_variation_logcounts.jpg',oct4_variation_logcounts, device='jpg', width = 8, height = 10)



# famous correlated & non-correlated genes(scatter plot) #######################
# positive control
targets <- c('Nanog', 'Sox2','Klf4','Zfp42','Utf1','Esrrb')
#rownames(sce)[grepl('Esrrb',rownames(sce))]

# raw counts
oct4_corr_scatter_rawcounts <- plotExpression(sce,targets,x='Pou5f1',exprs_values = "counts",color_by = 'Cell.Type') + 
  xlab('Oct4 expression level') + 
  ylab('Expression level (raw counts)') +
  ggtitle('Expression level of Oct4 and some of its target genes') +
  guides(color=guide_legend("Cell Type"))
oct4_corr_scatter_rawcounts
ggsave('figures/oct4_corr_scatter_rawcounts.jpg',oct4_corr_scatter_rawcounts, device='jpg', width = 8, height = 10)

# log counts
oct4_corr_scatter_logcounts <- plotExpression(sce,targets,x='Pou5f1',exprs_values = "logcounts",color_by = 'Cell.Type') + 
  xlab('Oct4 expression level') + 
  ylab('Expression level (log counts after normalisation)') +
  ggtitle('Expression level of Oct4 and some of its target genes (after normalisation)') +
  guides(color=guide_legend("Cell Type"))
oct4_corr_scatter_logcounts
ggsave('figures/oct4_corr_scatter_logcounts.jpg',oct4_corr_scatter_logcounts, device='jpg', width = 8, height = 10)
rm(targets)

# negative control
hk_genes <- c('Actb','Tbp','Pgk1','Ppia','Rpl38','Hmbs')

# raw counts
oct4_corr_scatter_rawcounts.hk <- plotExpression(sce,hk_genes,x='Pou5f1',exprs_values = "counts",color_by = 'Cell.Type') + 
  xlab('Oct4 expression level') + 
  ylab('Expression level (raw counts)') +
  ggtitle('Expression level of Oct4 and some housekeeping genes') +
  guides(color=guide_legend("Cell Type"))
oct4_corr_scatter_rawcounts.hk
ggsave('figures/oct4_corr_scatter_rawcounts_hk.jpg',oct4_corr_scatter_rawcounts.hk, device='jpg', width = 8, height = 10)

# log counts
oct4_corr_scatter_logcounts.hk <- plotExpression(sce,hk_genes,x='Pou5f1',exprs_values = "logcounts",color_by = 'Cell.Type') + 
  xlab('Oct4 expression level') + 
  ylab('Expression level (log counts after normalisation)') +
  ggtitle('Expression level of Oct4 and some housekeeping genes (after normalisation)') +
  guides(color=guide_legend("Cell Type"))
oct4_corr_scatter_logcounts.hk
ggsave('figures/oct4_corr_scatter_logcounts_hk.jpg',oct4_corr_scatter_logcounts.hk, device='jpg', width = 8, height = 10)



# fearture selection ###########################################################
assays(sce)$norm_counts <- 2^(logcounts(sce))-1
# scran HVG (model gene variance by trend)
#some how modelGeneVarWithSpikes(sce,'spikes') failed here...
var.out <- modelGeneVar(sce)
features_hvg <- getTopHVGs(var.out,n=50)
deconvolution_HVG <- ggplot(as.data.frame(var.out), aes(x=mean, y=total)) + geom_point() +
  geom_point(data=as.data.frame(var.out[features_hvg,c('mean','total')]), color='red') +
  geom_line(aes(x=mean,y=tech), color='skyblue') +
  geom_text_repel(label=features_hvg, 
            data=as.data.frame(var.out[features_hvg,c('mean','total')]),
            color='coral1', size=3,box.padding = unit(0.2, "lines"),
            point.padding = unit(0.25, "lines"), max.overlaps = 20) +
  ggtitle('Technicle variation (blue line) and top 50 HVG') +
  xlab('Mean log-expression') +
  ylab('Total Variance of log-expression')
deconvolution_HVG
ggsave('figures/deconvolution_HVG.jpg', deconvolution_HVG, device='jpg', width = 10, height = 6)

# log counts of HVG
features_hvg <- getTopHVGs(var.out,n=100)
plotExpression(sce,features_hvg,exprs_values = "logcounts") + ylim(0,20)
ggsave('figures/top100_hvg_log.jpg', device='jpg', width = 20, height = 6)

# normalised raw counts of HVG
plotExpression(sce,features_hvg,exprs_values = "norm_counts") + ylim(0,2e+5)
ggsave('figures/top100_hvg_raw.jpg',device='jpg', width = 20, height = 6)


# Deviance (select genes that are both highly expressed and highly variable)
dev <- rowData(devianceFeatureSelection(sce,assay = "counts", fam='binomial'))$binomial_deviance

# log counts of deviance
features_hvg <- names(dev[order(dev,decreasing=TRUE)])[1:100]
plotExpression(sce,features_hvg,exprs_values = "logcounts")+ ylim(0,20)
ggsave('figures/top100_dev_log.jpg',device='jpg', width = 20, height = 6)

# normalised raw counts of deviance
plotExpression(sce,features_hvg,exprs_values = "norm_counts") + ylim(0,2e+5)
ggsave('figures/top100_dev_raw.jpg', device='jpg', width = 20, height = 6)


# check intersection between features selected by HVG and DEV
features_hvg <- getTopHVGs(var.out,n=2000)
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:2000]
features_venn <- ggVennDiagram(data.frame(features_dev=features_dev,features_hvg=features_hvg),
                               label = "count", label_size = 5) +
  theme(legend.position = "none") +
  ggtitle('intersection of top 2000 features selected by Dev VS. HVG') +
  theme(plot.title = element_text(hjust = 0.5))
features_venn
ggsave("figures/features_venn.jpg",features_venn,device='jpg', width = 6, height = 6)

# mean-variance plot
meanvar_plotdata <- function(sce,G=2000){
  #this function is copied from: https://github.com/willtownes/scrna2019
  #sce=SingleCellExperiment with assays "counts" and "logcounts"
  #logcounts = log2(1+counts/scran size factor)
  #G=number of genes to be "highly informative"
  #assumes rowData(sce) contains columns "dev","hvg" with ranks
  #rank=1 means gene is most informative according to the criterion
  #returns a data frame to use for making a mean/variance plot
  gm<-as.data.frame(rowData(sce))
  rk<-gm[,c("dev","hvg")]
  rk<-rk[rownames(sce),]
  Y<-2^(logcounts(sce))-1
  pd<-data.frame(m=rowMeans(Y),v=apply(Y,1,var))
  pd$vmr=pd$v/pd$m
  pd<-cbind(pd,rk)
  pd<-subset(pd,vmr>0)
  xt<-table(pd$dev<=G,pd$hvg<=G)
  #formatting plot labels function
  f<-function(s,n){paste0(s," (",prettyNum(n,big.mark=","),")")}
  pd$criteria<-f("neither",xt[1,1])
  pd$criteria[pd$dev<=G & pd$hvg<=G]<-f("both",xt[2,2])
  pd$criteria[pd$dev<=G & pd$hvg>G]<-f("high deviance",xt[2,1])
  pd$criteria[pd$dev>G & pd$hvg<=G]<-f("highly variable",xt[1,2])
  #sort so that the "neither" category doesn't obscure the "highly deviant"
  pd[order(pd$criteria,decreasing=TRUE),]
}
test.sce <- sce
rowData(test.sce)$hvg <- rank(-var.out$bio)
rowData(test.sce)$dev <- rank(-dev)
#View(as.data.frame(rowData(test.sce)))
df <- meanvar_plotdata(test.sce, G=2000)
features_mean_var <- ggplot(df,aes(x=m,y=vmr,colour=criteria))+
  geom_point(alpha=.9)+xlab("mean normalized expression (no log transformation)") +
  ylab(latex2exp::TeX('$CV^2$ (variance/mean)')) +
  theme(legend.position=c(0.2,.8),
        plot.margin = margin(0.5, 0.5, 0.1, 0.1, "cm"))+
  scale_color_manual(values=c("orange","red","blue","gray"))+
  scale_x_log10()+scale_y_log10()+
  ggtitle('mean-variance plot of top 2000 features selected by Dev VS. HVG')
features_mean_var
ggsave("figures/features_mean_var.jpg",features_mean_var,device='jpg', width = 8, height = 7)
rm(test.sce)
rm(meanvar_plotdata)
rm(df)



# tSNE perplexity test #########################################################
# HVG:1000-4000 + perplexity:5,10,15,20,25,30,40,50 (8 rows, 4 cols) 
perplex_test.func <- function(sce,feat_num, perplex_num, ref, seed){
  # feat_num and perplex_num are both arrays of numbers
  r <- list()
  for(i in feat_num){
    print(i)
    features_hvg <- getTopHVGs(ref,n=i)
    p <- list()
    for(j in perplex_num){
      print(j)
      set.seed(seed)
      sce_hvg <- runPCA(sce,subset_row = features_hvg)
      set.seed(seed)
      sce_hvg <- runTSNE(sce_hvg,perplexity=j, dimred="PCA")
      p[[as.character(j)]] <- plotReducedDim(sce_hvg, "TSNE", colour_by = 'Cell.Type') +
                              ggtitle(paste0('Perplexity: ', j))
    }
    r[[as.character(i)]] <- annotate_figure(ggarrange(plotlist=p, nrow=length(perplex_num), legend='none'),
                    top=text_grob(paste0('top ',i,' HVG'), color = 'red',face = "bold"))
  }
  plot <- ggarrange(plotlist = r, ncol=length(feat_num),common.legend = TRUE, legend="right",
            legend.grob = get_legend(plotReducedDim(sce_hvg, "TSNE", colour_by = 'Cell.Type'))) +
          theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
  return(plot)
}
perplex_test.func(sce,c(1000,2000,3000,4000), c(5,10,15,20,25,30,40,50), var.out, seed=seed)
ggsave('figures/perplex_test.jpg',device='jpg', width = 20, height = 30)



# number of features test ######################################################
# by hvg
reduce_dim_hvg.func <- function(sce,num,ref, seed){
  # num is an array of numbers
  p <- list()
  for(n in num){
    print(n)
    features_hvg <- getTopHVGs(ref,n=n)
    set.seed(seed)
    sce_hvg <- runPCA(sce,subset_row = features_hvg)
    set.seed(seed)
    sce_hvg <- runTSNE(sce_hvg,perplexity=20, dimred="PCA")
    set.seed(seed)
    sce_hvg <- runUMAP(sce_hvg, subset_row = features_hvg)
    p[[as.character(n)]] <- annotate_figure(ggarrange(
      plotReducedDim(sce_hvg, "PCA", colour_by = 'Cell.Type')
      +ggtitle('PCA'),
      plotReducedDim(sce_hvg, "TSNE", colour_by = 'Cell.Type')
      +ggtitle('TSNE'),
      plotReducedDim(sce_hvg, "UMAP", colour_by = 'Cell.Type')
      +ggtitle('UMAP'),
      nrow=3,
      #common.legend = TRUE, 
      legend="none"
    ),
    top=text_grob(paste0('top ',n,' features'), color = 'red',face = "bold"))
  }
  r <- list()
  r[['plots']] <- p
  r[['legend']] <- get_legend(plotReducedDim(sce_hvg, "PCA", colour_by = 'Cell.Type'))
  return(r)
}
results <- reduce_dim_hvg.func(sce,c(100,500,1000,1500,2000,3000,4000,5000),ref=var.out, seed=seed)
dim_reduce_hvg <- ggarrange(plotlist=results[['plots']], 
                            ncol=4, nrow=2, common.legend = TRUE, legend="right", 
                            legend.grob = results[['legend']]) +
                  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
dim_reduce_hvg
ggsave('figures/dim_reduce_hvg.jpg',dim_reduce_hvg, device='jpg', width = 20, height = 25)
rm(results)

# by dev
reduce_dim_dev.func <- function(sce,num,ref, seed){
  # num is an array of numbers
  p <- list()
  for(n in num){
    print(n)
    features_dev <- names(ref[order(ref,decreasing=TRUE)])[1:n]
    set.seed(seed)
    sce_dev <- runPCA(sce,subset_row = features_dev)
    set.seed(seed)
    sce_dev <- runTSNE(sce_dev,perplexity=20, dimred="PCA")
    set.seed(seed)
    sce_dev <- runUMAP(sce_dev, subset_row = features_dev)
    p[[as.character(n)]] <- annotate_figure(ggarrange(
      plotReducedDim(sce_dev, "PCA", colour_by = 'Cell.Type')
      +ggtitle('PCA'),
      plotReducedDim(sce_dev, "TSNE", colour_by = 'Cell.Type')
      +ggtitle('TSNE'),
      plotReducedDim(sce_dev, "UMAP", colour_by = 'Cell.Type')
      +ggtitle('UMAP'),
      nrow=3,
      #common.legend = TRUE, 
      legend="none"
    ),
    top=text_grob(paste0('top ',n,' features'), color = 'red',face = "bold"))
  }
  r <- list()
  r[['plots']] <- p
  r[['legend']] <- get_legend(plotReducedDim(sce_dev, "PCA", colour_by = 'Cell.Type'))
  return(r)
}
results <- reduce_dim_dev.func(sce,c(100,500,1000,1500,2000,3000,4000,5000),ref=dev, seed=seed)
dim_reduce_dev <- ggarrange(plotlist=results[['plots']], 
                            ncol=4, nrow=2, common.legend = TRUE, legend="right", 
                            legend.grob = results[['legend']]) +
                  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
dim_reduce_dev
ggsave('figures/dim_reduce_dev.jpg',dim_reduce_dev, device='jpg', width = 20, height = 25)
rm(results)
# check seurat clustering part (sweep_para function) to choose the right number!



# # PC number selection methods ##################################################
# # here just select a feature number randomly (e.g. 2000)  
# # elbow
# features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:4000]
# #features_hvg <- getTopHVGs(var.out,n=4000)
# set.seed(seed)
# sce_dev <- runPCA(sce,subset_row = features_dev)
# percent.var <- attr(reducedDim(sce_dev), "percentVar")
# PC_num.elbow <- findElbowPoint(percent.var)
# PC_num.elbow
# # global maximum likelihood (based on translated Poisson mixture model)
# PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(reducedDim(sce_dev),k=10)
# PC_num.gml
# PC_num.gml <- round(PC_num.gml$dim.est,0)
# ggplot(data.frame(variance_percent = percent.var, PCs=1:length(percent.var))) + 
#   geom_point(aes(x=PCs, y=variance_percent)) +
#   geom_vline(xintercept = PC_num.elbow, color = 'blue', linetype="dashed") +
#   geom_vline(xintercept = PC_num.gml, color = 'coral', linetype="dashed") +
#   annotate('text',x=PC_num.elbow, y=5,label='elbow point',color = 'red') +
#   annotate('text',x=PC_num.gml, y=4,label='global maximum likelihood',color = 'red') +
#   ylab("Variance explained (%)")
# ggsave('figures/PC_num_sele.jpg',device='jpg', width = 8, height = 6)
# rm(percent.var)



# cluster cells (deviance feature) #############################################
# cluster by seurat package (using Leiden/Louvain algorithm)
#features_hvg <- getTopHVGs(var.out,n=5000)
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:5000]
sce_seurat <- sce_origin
sce_seurat <- sce_seurat[, !cell_filter$discard]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Pou5f1$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Sox2$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[!is.spike,]
sce_seurat <- sce_seurat[rowSums(assays(sce_seurat)$counts > 0) >= 3,]
sce_seurat
sce_seurat <- computeSumFactors(sce_seurat, cluster = quickCluster(sce_seurat))
sce_seurat <- logNormCounts(sce_seurat)
seurat <- as.Seurat(sce_seurat, counts = "counts", data = "logcounts")
#seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
#VariableFeatures(seurat)[1:100]
#VlnPlot(object = seurat, features = c('Pou5f1','Nanog'))
seurat <- ScaleData(seurat, features = rownames(seurat))

# test for resolution of leiden and features number (can also test other parameters!)
sweep_para <- function(seurat, seed, features, res, feature_num){
  para_res <- c()
  sil_score <- c()
  cluster_num <- c()
  ft_num <- c()
  PC_num <- c()
  for(j in feature_num){
    print(paste0('feature number = ', j))
    seurat <- RunPCA(seurat, features = features[1:j], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml
    PC_sele <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_sele)
    for(i in res){
      print(paste0('res = ', i))
      para_res <- c(para_res,i)
      seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
      cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
      sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
      sil.data <- as.data.frame(sil)
      score <- mean(sil.data$sil_width)
      sil_score <- c(sil_score,score)
      cluster_num <- c(cluster_num, length(table(cluster)))
      ft_num <- c(ft_num, j)
      PC_num <- c(PC_num, PC_sele)
    }
    cat('\n')
  }
  r <- data.frame(res=para_res, sil_score=sil_score, clust_num=cluster_num,
                  feature_num=ft_num, PC_num=PC_num)
  return(r)
}
out <- sweep_para(seurat, seed, features_dev, 
                  res=seq(0.1,2.0,0.1),
                  feature_num = c(100, 250, 500, seq(1000,5000,500)))
out
ggplot(out, aes(x=res,y=sil_score,label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_line(aes(color = factor(feature_num))) +
  ylab('Silhouette Score') +
  xlab('Leiden Resolution') +
  guides(color=guide_legend(title="Feature Number")) +
  scale_color_brewer(palette="Paired")
ggsave('figures/leiden_para_1_dev.jpg',device='jpg', width = 8, height = 6)
ggplot(out, aes(x = factor(feature_num), y = sil_score, 
                label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_bar(stat = "summary", fun = "var")+
  ylab('Variance of Silhouette') +
  xlab('Feature Number')
ggsave('figures/leiden_para_2_dev.jpg',device='jpg', width = 8, height = 5)
rm(out)

# comparison of different cluster results (different resolutions)
compare_para <- function(seurat, seed, resolution, feature_num){
  clusters <- list()
  for(i in c(1,2)){
    res <- resolution[i]
    ft_num <- feature_num[i]
    seurat <- RunPCA(seurat, features = features_dev[1:ft_num], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
    seurat <- FindClusters(seurat, resolution = res, algorithm = 4)
    name <- paste0('originalexp_snn_res.',res)
    cluster <- seurat@meta.data[[name]]
    clusters[[i]] <- cluster
  }
  return(clusters)
}
com_para_r <- compare_para(seurat, seed, c(0.7,0.7),c(2500, 4000))
table(com_para_r[[1]],com_para_r[[2]])
rm(com_para_r)

# select the best feature number and best PC number
feat_num_sele <- 4000
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:feat_num_sele]
seurat <- RunPCA(seurat, features = features_dev, seed.use = seed)  # here the features are chosen by deviance
#DimHeatmap(seurat, reduction = "pca",dims = 1:3)
PC_num.elbow <- findElbowPoint(seurat[['pca']]@stdev)
PC_num.elbow
PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
PC_num.gml
PC_num.gml <- round(PC_num.gml$dim.est,0)
PC_num.gml
ElbowPlot(seurat) +
  geom_vline(xintercept = PC_num.elbow, color = 'blue', linetype="dashed") +
  geom_vline(xintercept = PC_num.gml, color = 'coral', linetype="dashed") +
  annotate('text',x=PC_num.elbow, y=10,label='elbow point',color = 'red') +
  annotate('text',x=PC_num.gml, y=10,label='global maximum likelihood',color = 'red')
ggsave('figures/PC_num_sele.jpg',device='jpg', width = 10, height = 6)
seurat <- RunTSNE(seurat, seed.use = seed, perplexity = 20, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "tsne" , group.by = 'Cell.Type')
#ggsave('figures/leiden_tsne.jpg', device='jpg', width = 8, height = 7)
seurat <- RunUMAP(seurat,seed.use = seed, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "umap" , group.by = 'Cell.Type')

#  convey those dim-reduction plots to a SCE object (to use some functions in bioconductor)
sce_dev <- as.SingleCellExperiment(seurat)
plotReducedDim(sce_dev, "TSNE", colour_by="Cell.Type")  # make sure it's the same as seurat's tSNE
counts(sce_dev) <- as.matrix(counts(sce_dev))
assays(sce_dev)$norm_counts <- 2^(logcounts(sce_dev))-1
sce_dev

# select the best resolution and cluster cells using Leiden algorithm
res_sele <- 0.7
seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
seurat <- FindClusters(seurat, resolution = res_sele, algorithm = 4)  # algorithm 4 is "Leiden"; 1 is "Louvain"
DimPlot(seurat, reduction = "tsne", label = TRUE, shape.by = 'Cell.Type')

# pass the clustering result back to SCE object
cluster.leiden <- seurat@meta.data[[paste0('originalexp_snn_res.',res_sele)]]
colLabels(sce_dev) <- cluster.leiden

# plot marker genes of each cluster (heatmap)
group_marker <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
group_marker %>% dplyr::filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 10, order_by = p_val_adj) -> top10
DoHeatmap(seurat,features=top10$gene, slot = "scale.data") #+ scale_fill_virdis()
#FeaturePlot(seurat, reduction = 'TSNE', features=top10[top10['cluster']==1,]$gene[1:4])
leiden_markers <- plotHeatmap(sce_dev, exprs_values = "logcounts", 
                              order_columns_by=c("label", "Cell.Type"), 
                              features=top10$gene, cluster_rows = FALSE, center = TRUE,
                              gaps_col = cumsum(as.numeric(table(colLabels(sce_dev)))),
                              gaps_row = seq(0,50,10), 
                              color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40))
ggsave('figures/leiden_markers.jpg', leiden_markers, device='jpg', width = 10, height = 8)
#colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(100)
#colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)
#colorRampPalette(c(rev(RColorBrewer::brewer.pal(9,'Blues')), RColorBrewer::brewer.pal(9,'Reds')))(100)
#colorRampPalette(c('#fc00fc','black','#fcfc00'))(100)

# # Find markers between selected two groups (here we use 2 clsuters in 'serum' group)
# #FindMarkers(seurat, ident.1 = 1, min.pct = 0.25, only.pos = TRUE)  # compare to the rest of cells
# group_diff <- FindMarkers(seurat, ident.1 = 3, ident.2 = 5, min.pct = 0.25)
# head(group_diff)
# group_diff <- group_diff[group_diff$p_val_adj < 0.01, ]
# nrow(group_diff)
# write.csv(group_diff,'figures/group_diff_serum.csv',row.names = TRUE)
# 
# # show difference of expression on tSNE plot
# FeaturePlot(seurat, features=rownames(group_diff[1:20,]))
# ggsave('figures/group_diff_serum_1.jpg', device='jpg', width = 17, height = 15)
# 
# # seperate into two groups based on positive/negative fold change
# up_reg <- group_diff[group_diff$avg_log2FC > 0,]
# down_reg <- group_diff[group_diff$avg_log2FC < 0,]
# 
# # show difference of expression on heatmap
# #DoHeatmap(seurat,features=rownames(group_diff[1:50,]), slot = "scale.data")
# group_diff_serum_2 <- plotHeatmap(sce_dev[,sce_dev$Cell.Type=='serum'], 
#                                   exprs_values = "logcounts", order_columns_by=c("label", "Cell.Type"), 
#                                   features=rownames(rbind(up_reg[1:25,], down_reg[1:25,])), 
#                                   cluster_rows = FALSE, center = TRUE, zlim = c(-9,9), 
#                                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40), 
#                                   gaps_col = cumsum(as.numeric(table(colLabels(sce_dev)))[c(3,5)]),
#                                   gaps_row = c(25,50))
# ggsave('figures/group_diff_serum_2.jpg', group_diff_serum_2, device='jpg', width = 12, height = 10)

# details of Silhouette width (and tSNE plots of clustering results)
#sil <- cluster::silhouette(as.integer(cluster.ld), dist(reducedDim(sce_dev, "PCA")))
detail_sil <- function(res, seurat, PC_sele){
  for(i in res){
    print(paste0('res = ', i))
    slot_name <- paste0('originalexp_snn_res.',i)
    dir1 <- paste0('figures/leiden_res_',i,'_1.jpg')
    dir2 <- paste0('figures/leiden_res_',i,'_2.jpg')
    dir3 <- paste0('figures/leiden_res_',i,'_tsne.jpg')
    seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
    cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
    sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
    sil.data <- as.data.frame(sil)
    sil.data$closest <- factor(ifelse(sil.data$sil_width > 0, cluster, sil.data$neighbor))
    #sil.data$cluster <- cluster.ld
    ggplot(sil.data, aes(x=cluster, y=sil_width, colour=closest)) +
      ggbeeswarm::geom_quasirandom(method="smiley") + 
      ylab('Silhouette width') +
      geom_hline(aes(yintercept=0), color = 'red', linetype="dashed")
    ggsave(dir1, device='jpg', width = 6, height = 5)
    ggplot(as.data.frame(sil.data), aes(x=sil_width))+
      geom_histogram(color='grey80',bins=50) +
      geom_density(alpha=.2, fill="blue") + xlab('Silhouette Width')+
      ggtitle(paste0('avarage Silhouette Width: ',round(mean(sil.data$sil_width),4))) +
      geom_vline(aes(xintercept=0), color = 'red', linetype="dashed")
    ggsave(dir2, device='jpg', width = 5, height = 4)
    DimPlot(seurat, reduction = "tsne" ,label = TRUE, shape.by = 'Cell.Type')
    ggsave(dir3, device='jpg', width = 8, height = 6)
  }
}
detail_sil(res=c(0.7,1.2,1.7), seurat = seurat, PC_sele=PC_num.gml)

# clean up variable
rm(group_diff)
rm(group_marker)
rm(top10)
rm(group_diff_serum_2)
rm(leiden_markers)
rm(up_reg)
rm(down_reg)
rm(sce_seurat)



# # cluster by scran (using Leiden/Louvain/walktrap algorithm)
# cluster.ld <- clusterCells(sce_dev, use.dimred="PCA", BLUSPARAM=NNGraphParam(k=10, cluster.fun='leiden', cluster.args = list(resolution = 0.8)))
# table(cluster.ld)
# cluster.wt <- clusterCells(sce_dev, use.dimred="PCA", BLUSPARAM=NNGraphParam(k=10, type="rank", cluster.fun="walktrap"))
# table(cluster.wt)
# tab <- table(Walktrap=cluster.wt, Leiden=cluster.ld)  # comparison
# plot
# colLabels(sce_dev) <- cluster.ld
# plotReducedDim(sce_dev, "TSNE", colour_by="label",text_by="label")
# 
# # Find Markers by scran
# marker.info <- scoreMarkers(sce_dev, cluster.ld)
# # chose the upregulated marker genes
# #head(rownames(marker.info[[2]])[order(marker.info[[2]]$mean.AUC, decreasing=TRUE)],10)
# # chose the downregulated marker genes
# head(rownames(marker.info[[2]])[order(marker.info[[2]]$mean.AUC, decreasing=FALSE)],10)
# rm(marker.info)


# # cluster by SC3 package (k-mean based)
# sce_dev_SC3 <- sce_dev
# rowData(sce_dev_SC3)$feature_symbol <- rownames(sce_dev_SC3)
# sce_dev_SC3
# sce_dev_SC3 <- sc3(sce_dev_SC3, gene_filter=TRUE, ks = 3:8, biology = TRUE, rand_seed=seed,
#                    d_region_min=0.01)
# #sc3_plot_consensus(sce_dev_SC3, k = 6, show_pdata = c("Cell.Type", "sc3_6_clusters"))
# 
# # evaluate the best clustering number for SC3 (check the Silhouette score in picture!)
# eva_SC3_cluster <- function(sce,cluster_num){
#   for(i in cluster_num){
#     filename <- paste0('figures/SC3_sil_',i,'.png')
#     png(filename,1000,1000)
#     sc3_plot_silhouette(sce, k = i)
#     dev.off()
#   }
# }
# eva_SC3_cluster(sce_dev_SC3,3:8)
# #sc3_plot_cluster_stability(sce_dev_SC3, k = 6)
# 
# # plot cluster result
# plotTSNE(sce_dev_SC3,shape_by = 'Cell.Type', colour_by = 'sc3_6_clusters')
# ggsave('figures/SC3_tsne_6_clust.jpg', device='jpg', width = 6, height = 4.5)
# #sc3_plot_expression(sce_dev_SC3, k = 6,show_pdata = c("sc3_6_clusters","Cell.Type"))
# 
# # markers of each cluster
# #sc3_plot_de_genes(sce_dev_SC3, k = 6,show_pdata = c("sc3_6_clusters", "Cell.Type"))
# SC3_markerplot <- sc3_plot_markers(sce_dev_SC3, k = 6, auroc = 0.8, p.val = 0.01,
#                                    show_pdata = c("sc3_6_clusters", "Cell.Type"))
# ggsave('figures/SC3_markerplot_6_clust.jpg', SC3_markerplot, device='jpg', width = 10, height = 10)
# cluster.sc3 <- colData(sce_dev_SC3)$sc3_6_clusters
# colLabels(sce_dev_SC3) <- cluster.sc3
# rm(SC3_markerplot)
# 
# # Find Markers by SC3 (not as good as seurat)
# find_marker_sc3 <- function(sce, k, padj, auroc, clust){
#   clust_name <- paste0('sc3_',k,'_markers_clusts')
#   padj_name <- paste0('sc3_',k,'_markers_padj')
#   auroc_name <- paste0('sc3_',k,'_markers_auroc')
#   df <- rowData(sce)[,c(clust_name, padj_name, auroc_name)]
#   df[is.na(df[[clust_name]]) |
#      is.na(df[[padj_name]]) |
#      is.na(df[[auroc_name]]),] <- 0
#   df <- df[df[[clust_name]] == clust &
#            df[[padj_name]] < padj &
#            df[[auroc_name]] > auroc,]
#   df <- df[order(df[[padj_name]], decreasing = FALSE),]
#   markers <- rownames(df)
#   return(markers)
# }
# sc3_markers <- find_marker_sc3(sce_dev_SC3, 6, padj=0.01, auroc=0.7, clust=1)
# sc3_markers
# plotHeatmap(sce_dev_SC3, exprs_values = "logcounts",
#             order_columns_by=c("label", "Cell.Type"), features=c(sc3_markers),
#             cluster_rows = FALSE, center=TRUE)



# features selected by HVG #####################################################
features_hvg <- getTopHVGs(var.out,n=5000)
#features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:5000]
sce_seurat <- sce_origin
sce_seurat <- sce_seurat[, !cell_filter$discard]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Pou5f1$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Sox2$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[!is.spike,]
sce_seurat <- sce_seurat[rowSums(assays(sce_seurat)$counts > 0) >= 3,]
sce_seurat
sce_seurat <- computeSumFactors(sce_seurat, cluster = quickCluster(sce_seurat))
sce_seurat <- logNormCounts(sce_seurat)
seurat <- as.Seurat(sce_seurat, counts = "counts", data = "logcounts")
#seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
#VariableFeatures(seurat)[1:100]
#VlnPlot(object = seurat, features = c('Pou5f1','Nanog'))
seurat <- ScaleData(seurat, features = rownames(seurat))

# test for resolution of leiden and features number (can also test other parameters!)
sweep_para <- function(seurat, seed, features, res, feature_num){
  para_res <- c()
  sil_score <- c()
  cluster_num <- c()
  ft_num <- c()
  PC_num <- c()
  for(j in feature_num){
    print(paste0('feature number = ', j))
    seurat <- RunPCA(seurat, features = features[1:j], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml
    PC_sele <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_sele)
    for(i in res){
      print(paste0('res = ', i))
      para_res <- c(para_res,i)
      seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
      cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
      sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
      sil.data <- as.data.frame(sil)
      score <- mean(sil.data$sil_width)
      sil_score <- c(sil_score,score)
      cluster_num <- c(cluster_num, length(table(cluster)))
      ft_num <- c(ft_num, j)
      PC_num <- c(PC_num, PC_sele)
    }
    cat('\n')
  }
  r <- data.frame(res=para_res, sil_score=sil_score, clust_num=cluster_num,
                  feature_num=ft_num, PC_num=PC_num)
  return(r)
}
out <- sweep_para(seurat, seed, features_hvg, 
                  res=seq(0.1,2.0,0.1),
                  feature_num = c(100, 250, 500, seq(1000,5000,500)))
out
ggplot(out, aes(x=res,y=sil_score,label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_line(aes(color = factor(feature_num))) +
  ylab('Silhouette Score') +
  xlab('Leiden Resolution') +
  guides(color=guide_legend(title="Feature Number")) +
  scale_color_brewer(palette="Paired")
ggsave('figures/leiden_para_1_hvg.jpg',device='jpg', width = 8, height = 6)
ggplot(out, aes(x = factor(feature_num), y = sil_score, 
                label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_bar(stat = "summary", fun = "var")+
  ylab('Variance of Silhouette') +
  xlab('Feature Number')
ggsave('figures/leiden_para_2_hvg.jpg',device='jpg', width = 8, height = 5)
rm(out)

# comparison of different cluster results (different resolutions)
compare_para <- function(seurat, seed, resolution, feature_num){
  clusters <- list()
  for(i in c(1,2)){
    res <- resolution[i]
    ft_num <- feature_num[i]
    seurat <- RunPCA(seurat, features = features_hvg[1:ft_num], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
    seurat <- FindClusters(seurat, resolution = res, algorithm = 4)
    name <- paste0('originalexp_snn_res.',res)
    cluster <- seurat@meta.data[[name]]
    clusters[[i]] <- cluster
  }
  return(clusters)
}
com_para_r <- compare_para(seurat, seed, c(0.7,0.7),c(2000, 4000))
table(com_para_r[[1]],com_para_r[[2]])
rm(com_para_r)

# select the best feature number and best PC number
feat_num_sele <- 4000
features_hvg <- getTopHVGs(var.out,n=feat_num_sele)
seurat <- RunPCA(seurat, features = features_hvg, seed.use = seed)  # here the features are chosen by deviance
#DimHeatmap(seurat, reduction = "pca",dims = 1:3)
PC_num.elbow <- findElbowPoint(seurat[['pca']]@stdev)
PC_num.elbow
PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
PC_num.gml
PC_num.gml <- round(PC_num.gml$dim.est,0)
PC_num.gml
ElbowPlot(seurat) +
  geom_vline(xintercept = PC_num.elbow, color = 'blue', linetype="dashed") +
  geom_vline(xintercept = PC_num.gml, color = 'coral', linetype="dashed") +
  annotate('text',x=PC_num.elbow, y=10,label='elbow point',color = 'red') +
  annotate('text',x=PC_num.gml, y=10,label='global maximum likelihood',color = 'red')
ggsave('figures/PC_num_sele.jpg',device='jpg', width = 10, height = 6)
seurat <- RunTSNE(seurat, seed.use = seed, perplexity = 20, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "tsne" , group.by = 'Cell.Type')
#ggsave('figures/leiden_tsne.jpg', device='jpg', width = 8, height = 7)
seurat <- RunUMAP(seurat,seed.use = seed, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "umap" , group.by = 'Cell.Type')

#  convey those dim-reduction plots to a SCE object (to use some functions in bioconductor)
sce_hvg <- as.SingleCellExperiment(seurat)
plotReducedDim(sce_hvg, "TSNE", colour_by="Cell.Type")  # make sure it's the same as seurat's tSNE
counts(sce_hvg) <- as.matrix(counts(sce_hvg))
assays(sce_hvg)$norm_counts <- 2^(logcounts(sce_hvg))-1
sce_hvg

# select the best resolution and cluster cells using Leiden algorithm
res_sele <- 0.7
seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
seurat <- FindClusters(seurat, resolution = res_sele, algorithm = 4)  # algorithm 4 is "Leiden"; 1 is "Louvain"
DimPlot(seurat, reduction = "tsne", label = TRUE, shape.by = 'Cell.Type')

# pass the clustering result back to SCE object
cluster.leiden <- seurat@meta.data[[paste0('originalexp_snn_res.',res_sele)]]
colLabels(sce_hvg) <- cluster.leiden

# compare to Deviance results
plotReducedDim(sce_hvg, "TSNE", shape_by="Cell.Type", colour_by = 'label') 
ggsave('figures/leiden_tsne_hvg.jpg', device='jpg', width = 8, height = 7)
sce_hvg$dev_clust <- colLabels(sce_dev)
plotReducedDim(sce_hvg, "TSNE", shape_by="Cell.Type", colour_by = 'dev_clust')
ggsave('figures/leiden_tsne_hvg_dev.jpg', device='jpg', width = 8, height = 7)
table(dev=sce_hvg$dev_clust, hvg=colLabels(sce_hvg))

# plot marker genes of each cluster (heatmap)
group_marker <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
group_marker %>% dplyr::filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 10, order_by = p_val_adj) -> top10
DoHeatmap(seurat,features=top10$gene, slot = "scale.data") #+ scale_fill_virdis()
#FeaturePlot(seurat, reduction = 'TSNE', features=top10[top10['cluster']==1,]$gene[1:4])
leiden_markers <- plotHeatmap(sce_hvg, exprs_values = "logcounts", 
                              order_columns_by=c("label", "Cell.Type"), 
                              features=top10$gene, cluster_rows = FALSE, center = TRUE,
                              gaps_col = cumsum(as.numeric(table(colLabels(sce_hvg)))),
                              gaps_row = seq(0,50,10), 
                              color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40))
ggsave('figures/leiden_markers_hvg.jpg', leiden_markers, device='jpg', width = 10, height = 8)
#colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(100)
#colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)
#colorRampPalette(c(rev(RColorBrewer::brewer.pal(9,'Blues')), RColorBrewer::brewer.pal(9,'Reds')))(100)
#colorRampPalette(c('#fc00fc','black','#fcfc00'))(100)

# Find markers between selected two groups (here we use 2 clsuters in 'serum' group)
#FindMarkers(seurat, ident.1 = 1, min.pct = 0.25, only.pos = TRUE)  # compare to the rest of cells
group_diff <- FindMarkers(seurat, ident.1 = 2, ident.2 = 5, min.pct = 0.25)
head(group_diff)
group_diff <- group_diff[group_diff$p_val_adj < 0.01, ]
nrow(group_diff)
write.csv(group_diff,'figures/group_diff_serum_hvg.csv',row.names = TRUE)

# show difference of expression on tSNE plot
FeaturePlot(seurat, features=rownames(group_diff[1:20,]))
ggsave('figures/group_diff_serum_1.jpg', device='jpg', width = 17, height = 15)

# seperate into two groups based on positive/negative fold change
up_reg <- group_diff[group_diff$avg_log2FC > 0,]
down_reg <- group_diff[group_diff$avg_log2FC < 0,]

# show difference of expression on heatmap
#DoHeatmap(seurat,features=rownames(group_diff[1:50,]), slot = "scale.data")
group_diff_serum_2 <- plotHeatmap(sce_hvg[,sce_hvg$Cell.Type=='serum'],
                                  exprs_values = "logcounts", order_columns_by=c("label", "Cell.Type"),
                                  features=rownames(rbind(up_reg[1:25,], down_reg[1:25,])),
                                  cluster_rows = FALSE, center = TRUE, zlim = c(-9,9),
                                  color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40),
                                  gaps_col = cumsum(as.numeric(table(colLabels(sce_hvg)))[c(2,5)]),
                                  gaps_row = c(25,50))
ggsave('figures/group_diff_serum_hvg.jpg', group_diff_serum_2, device='jpg', width = 12, height = 10)

# details of Silhouette width (and tSNE plots of clustering results)
#sil <- cluster::silhouette(as.integer(cluster.ld), dist(reducedDim(sce_hvg, "PCA")))
detail_sil <- function(res, seurat, PC_sele){
  for(i in res){
    print(paste0('res = ', i))
    slot_name <- paste0('originalexp_snn_res.',i)
    dir1 <- paste0('figures/leiden_hvg_res_',i,'_1.jpg')
    dir2 <- paste0('figures/leiden_hvg_res_',i,'_2.jpg')
    dir3 <- paste0('figures/leiden_hvg_res_',i,'_tsne.jpg')
    seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
    cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
    sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
    sil.data <- as.data.frame(sil)
    sil.data$closest <- factor(ifelse(sil.data$sil_width > 0, cluster, sil.data$neighbor))
    #sil.data$cluster <- cluster.ld
    ggplot(sil.data, aes(x=cluster, y=sil_width, colour=closest)) +
      ggbeeswarm::geom_quasirandom(method="smiley") + 
      ylab('Silhouette width') +
      geom_hline(aes(yintercept=0), color = 'red', linetype="dashed")
    ggsave(dir1, device='jpg', width = 6, height = 5)
    ggplot(as.data.frame(sil.data), aes(x=sil_width))+
      geom_histogram(color='grey80',bins=50) +
      geom_density(alpha=.2, fill="blue") + xlab('Silhouette Width')+
      ggtitle(paste0('avarage Silhouette Width: ',round(mean(sil.data$sil_width),4))) +
      geom_vline(aes(xintercept=0), color = 'red', linetype="dashed")
    ggsave(dir2, device='jpg', width = 5, height = 4)
    DimPlot(seurat, reduction = "tsne" ,label = TRUE, shape.by = 'Cell.Type')
    ggsave(dir3, device='jpg', width = 8, height = 6)
  }
}
detail_sil(res=c(0.3,0.7), seurat = seurat, PC_sele=PC_num.gml)

# clean up variable
rm(group_diff)
rm(group_marker)
rm(top10)
rm(group_diff_serum_2)
rm(leiden_markers)
rm(up_reg)
rm(down_reg)
rm(sce_seurat)


# normalisation by SCtransform #################################################
# cluster by seurat package (using Leiden/Louvain algorithm)
#features_hvg <- getTopHVGs(var.out,n=5000)
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:5000]
sce_seurat <- sce_origin
sce_seurat <- sce_seurat[, !cell_filter$discard]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Pou5f1$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Sox2$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[!is.spike,]
sce_seurat <- sce_seurat[rowSums(assays(sce_seurat)$counts > 0) >= 3,]
sce_seurat
seurat <- as.Seurat(sce_seurat, data = "counts")
seurat <- SCTransform(seurat, vst.flavor = "v2", assay = 'originalexp')
#seurat2 <- PercentageFeatureSet(seurat2, pattern = "^mt.", col.name = "percent.mt")
#SCTransform(seurat2, method = "glmGamPoi", vars.to.regress = "percent.mt", assay = 'originalexp')
#seurat <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
#VariableFeatures(seurat)[1:100]
#VlnPlot(object = seurat, features = c('Pou5f1','Nanog'))

# test for resolution of leiden and features number (can also test other parameters!)
sweep_para <- function(seurat, seed, features, res, feature_num){
  para_res <- c()
  sil_score <- c()
  cluster_num <- c()
  ft_num <- c()
  PC_num <- c()
  for(j in feature_num){
    print(paste0('feature number = ', j))
    seurat <- RunPCA(seurat, features = features[1:j], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml
    PC_sele <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_sele)
    for(i in res){
      print(paste0('res = ', i))
      para_res <- c(para_res,i)
      seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
      cluster <- seurat@meta.data[[paste0('SCT_snn_res.',i)]]
      sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
      sil.data <- as.data.frame(sil)
      score <- mean(sil.data$sil_width)
      sil_score <- c(sil_score,score)
      cluster_num <- c(cluster_num, length(table(cluster)))
      ft_num <- c(ft_num, j)
      PC_num <- c(PC_num, PC_sele)
    }
    cat('\n')
  }
  r <- data.frame(res=para_res, sil_score=sil_score, clust_num=cluster_num,
                  feature_num=ft_num, PC_num=PC_num)
  return(r)
}
out <- sweep_para(seurat, seed, features_dev, 
                  res=seq(0.1,2.0,0.1),
                  feature_num = c(100, 250, 500, seq(1000,5000,500)))
out
ggplot(out, aes(x=res,y=sil_score,label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_line(aes(color = factor(feature_num))) +
  ylab('Silhouette Score') +
  xlab('Leiden Resolution') +
  guides(color=guide_legend(title="Feature Number")) +
  scale_color_brewer(palette="Paired")
ggsave('figures/leiden_para_1_sct_dev.jpg',device='jpg', width = 8, height = 6)
ggplot(out, aes(x = factor(feature_num), y = sil_score, 
                label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_bar(stat = "summary", fun = "var")+
  ylab('Variance of Silhouette') +
  xlab('Feature Number')
ggsave('figures/leiden_para_2_sct_dev.jpg',device='jpg', width = 8, height = 5)
rm(out)

# comparison of different cluster results (different resolutions)
compare_para <- function(seurat, seed, resolution, feature_num){
  clusters <- list()
  for(i in c(1,2)){
    res <- resolution[i]
    ft_num <- feature_num[i]
    seurat <- RunPCA(seurat, features = features_dev[1:ft_num], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
    seurat <- FindClusters(seurat, resolution = res, algorithm = 4)
    name <- paste0('SCT_snn_res.',res)
    cluster <- seurat@meta.data[[name]]
    clusters[[i]] <- cluster
  }
  return(clusters)
}
com_para_r <- compare_para(seurat, seed, c(0.7,0.7),c(2000, 4000))
table(com_para_r[[1]],com_para_r[[2]])
rm(com_para_r)

# select the best feature number and best PC number
feat_num_sele <- 4000
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:feat_num_sele]
seurat <- RunPCA(seurat, features = features_dev, seed.use = seed)  # here the features are chosen by deviance
#DimHeatmap(seurat, reduction = "pca",dims = 1:3)
PC_num.elbow <- findElbowPoint(seurat[['pca']]@stdev)
PC_num.elbow
PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
PC_num.gml
PC_num.gml <- round(PC_num.gml$dim.est,0)
PC_num.gml
ElbowPlot(seurat) + 
  geom_vline(xintercept = PC_num.elbow, color = 'blue', linetype="dashed") +
  geom_vline(xintercept = PC_num.gml, color = 'coral', linetype="dashed") +
  annotate('text',x=PC_num.elbow, y=10,label='elbow point',color = 'red') +
  annotate('text',x=PC_num.gml, y=10,label='global maximum likelihood',color = 'red')
#ggsave('figures/PC_num_sele.jpg',device='jpg', width = 10, height = 6)
seurat <- RunTSNE(seurat, seed.use = seed, perplexity = 20, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "tsne" , group.by = 'Cell.Type')
#ggsave('figures/leiden_tsne.jpg', device='jpg', width = 8, height = 7)
seurat <- RunUMAP(seurat,seed.use = seed, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "umap" , group.by = 'Cell.Type')

#  convey those dim-reduction plots to a SCE object (to use some functions in bioconductor)
sce_dev_sct <- as.SingleCellExperiment(seurat)
plotReducedDim(sce_dev_sct, "TSNE", colour_by="Cell.Type")  # make sure it's the same as seurat's tSNE
counts(sce_dev_sct) <- as.matrix(counts(sce_dev_sct))
assays(sce_dev_sct)$norm_counts <- 2^(logcounts(sce_dev_sct))-1
sce_dev_sct

# select the best resolution and cluster cells using Leiden algorithm
res_sele <- 0.7
seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
seurat <- FindClusters(seurat, resolution = res_sele, algorithm = 4)  # algorithm 4 is "Leiden"; 1 is "Louvain"
DimPlot(seurat, reduction = "tsne", label = TRUE, shape.by = 'Cell.Type')

# pass the clustering result back to SCE object
cluster.leiden <- seurat@meta.data[[paste0('SCT_snn_res.',res_sele)]]
colLabels(sce_dev_sct) <- cluster.leiden

# compare to deconvolution results
plotReducedDim(sce_dev_sct, "TSNE", shape_by="Cell.Type", colour_by = 'label') 
ggsave('figures/leiden_tsne_sct.jpg', device='jpg', width = 8, height = 7)
sce_dev_sct$dev_clust <- colLabels(sce_dev)
plotReducedDim(sce_dev_sct, "TSNE", shape_by="Cell.Type", colour_by = 'dev_clust')
ggsave('figures/leiden_tsne_sct_decon.jpg', device='jpg', width = 8, height = 7)
table(decon=sce_dev_sct$dev_clust,sct=colLabels(sce_dev_sct))

# plot marker genes of each cluster (heatmap)
group_marker <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
group_marker %>% dplyr::filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 10, order_by = p_val_adj) -> top10
DoHeatmap(seurat,features=top10$gene, slot = "scale.data") #+ scale_fill_virdis()
#FeaturePlot(seurat, reduction = 'TSNE', features=top10[top10['cluster']==1,]$gene[1:4])
leiden_markers <- plotHeatmap(sce_dev_sct, exprs_values = "logcounts", 
                              order_columns_by=c("label", "Cell.Type"), 
                              features=top10$gene, cluster_rows = FALSE, center = TRUE,
                              gaps_col = cumsum(as.numeric(table(colLabels(sce_dev_sct)))),
                              gaps_row = seq(0,50,10), zlim= c(-6,6),
                              color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(30))
ggsave('figures/leiden_markers_sct.jpg', leiden_markers, device='jpg', width = 10, height = 8)
#colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(100)
#colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(100)
#colorRampPalette(c(rev(RColorBrewer::brewer.pal(9,'Blues')), RColorBrewer::brewer.pal(9,'Reds')))(100)
#colorRampPalette(c('#fc00fc','black','#fcfc00'))(100)

# Find markers between selected two groups (here we use 2 clsuters in 'serum' group)
#FindMarkers(seurat, ident.1 = 1, min.pct = 0.25, only.pos = TRUE)  # compare to the rest of cells
group_diff <- FindMarkers(seurat, ident.1 = 2, ident.2 = 5, min.pct = 0.25)
head(group_diff)
group_diff <- group_diff[group_diff$p_val_adj < 0.01, ]
nrow(group_diff)
write.csv(group_diff,'figures/group_diff_serum_sct.csv',row.names = TRUE)

# show difference of expression on tSNE plot
FeaturePlot(seurat, features=rownames(group_diff[1:20,]))
ggsave('figures/group_diff_serum_1.jpg', device='jpg', width = 17, height = 15)

# seperate into two groups based on positive/negative fold change
up_reg <- group_diff[group_diff$avg_log2FC > 0,]
down_reg <- group_diff[group_diff$avg_log2FC < 0,]

# show difference of expression on heatmap
#DoHeatmap(seurat,features=rownames(group_diff[1:50,]), slot = "scale.data")
group_diff_serum_2 <- plotHeatmap(sce_dev_sct[,sce_dev_sct$Cell.Type=='serum'],
                                  exprs_values = "logcounts", order_columns_by=c("label", "Cell.Type"),
                                  features=rownames(rbind(up_reg[1:25,], down_reg[1:25,])),
                                  cluster_rows = FALSE, center = TRUE, zlim = c(-7,7),
                                  color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(30),
                                  gaps_col = cumsum(as.numeric(table(colLabels(sce_dev_sct)))[c(2,5)]),
                                  gaps_row = c(25,50))
ggsave('figures/group_diff_serum_sct.jpg', group_diff_serum_2, device='jpg', width = 12, height = 10)

# details of Silhouette width (and tSNE plots of clustering results)
#sil <- cluster::silhouette(as.integer(cluster.ld), dist(reducedDim(sce_dev_sct, "PCA")))
detail_sil <- function(res, seurat, PC_sele){
  for(i in res){
    print(paste0('res = ', i))
    slot_name <- paste0('originalexp_snn_res.',i)
    dir1 <- paste0('figures/leiden_sct_res_',i,'_1.jpg')
    dir2 <- paste0('figures/leiden_sct_res_',i,'_2.jpg')
    dir3 <- paste0('figures/leiden_sct_res_',i,'_tsne.jpg')
    seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
    cluster <- seurat@meta.data[[paste0('SCT_snn_res.',i)]]
    sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
    sil.data <- as.data.frame(sil)
    sil.data$closest <- factor(ifelse(sil.data$sil_width > 0, cluster, sil.data$neighbor))
    #sil.data$cluster <- cluster.ld
    ggplot(sil.data, aes(x=cluster, y=sil_width, colour=closest)) +
      ggbeeswarm::geom_quasirandom(method="smiley") + 
      ylab('Silhouette width') +
      geom_hline(aes(yintercept=0), color = 'red', linetype="dashed")
    ggsave(dir1, device='jpg', width = 6, height = 5)
    ggplot(as.data.frame(sil.data), aes(x=sil_width))+
      geom_histogram(color='grey80',bins=50) +
      geom_density(alpha=.2, fill="blue") + xlab('Silhouette Width')+
      ggtitle(paste0('avarage Silhouette Width: ',round(mean(sil.data$sil_width),4))) +
      geom_vline(aes(xintercept=0), color = 'red', linetype="dashed")
    ggsave(dir2, device='jpg', width = 5, height = 4)
    DimPlot(seurat, reduction = "tsne" ,label = TRUE, shape.by = 'Cell.Type')
    ggsave(dir3, device='jpg', width = 8, height = 6)
  }
}
detail_sil(res=c(0.3,0.7), seurat = seurat, PC_sele=PC_num.gml)

# clean up variable
rm(group_diff)
rm(group_marker)
rm(top10)
rm(group_diff_serum_2)
rm(leiden_markers)
rm(up_reg)
rm(down_reg)
rm(sce_seurat)



# Spearman Correlation of all Oct4 potential targets############################
# load potential Oct4 targets
oct4_tar.potent <- read.csv('oct4_targets.csv', na.strings = c("", "NA"))
oct4_tar.potent <- na.omit(oct4_tar.potent[,1:2])
oct4_tar.potent <- oct4_tar.potent$Gene.name
table(oct4_tar.potent %in% rownames(sce_dev))  
oct4_tar.potent[!oct4_tar.potent %in% rownames(sce_dev)]  # check which genes are not available
oct4_tar.potent <- c('Pou5f1',oct4_tar.potent[oct4_tar.potent %in% rownames(sce_dev)])
tar_length <- length(oct4_tar.potent)
tar_length

# calculate Spearman correlation
oct4_tar.corr <- correlatePairs(sce_dev[,sce_dev$Cell.Type == 'serum'], subset.row=oct4_tar.potent)
oct4_tar.corr <- oct4_tar.corr[apply(is.na(oct4_tar.corr),1,sum)==0,]
oct4_tar.corr

# build correlation matrix
build_corr <- function(corr){
  df <- as.data.frame(corr[c('gene1','gene2','rho')])
  df2 <- data.frame(gene1=df$gene2, gene2=df$gene1, rho=df$rho)
  df <- rbind(df,df2)
  corr.matrix <- spread(df, key = gene2, value = rho)
  rownames(corr.matrix) <- corr.matrix$gene1
  corr.matrix$gene1 <- NULL  # remove 1st column
  diag(corr.matrix) <- 1
  return(corr.matrix)
}
oct4_tar_corr.matrix <- build_corr(oct4_tar.corr)

# df <- as.data.frame(oct4_tar.corr[c('gene1','gene2','rho')])
# df2 <- data.frame(gene1=df$gene2, gene2=df$gene1, rho=df$rho)
# df <- rbind(df,df2)
# oct4_tar_corr.matrix <- spread(df, key = gene2, value = rho)
# rownames(oct4_tar_corr.matrix) <- oct4_tar_corr.matrix$gene1
# oct4_tar_corr.matrix$gene1 <- NULL  # remove 1st column
# diag(oct4_tar_corr.matrix) <- 1

# plot correlation
oct4_potent_tar <- pheatmap(oct4_tar_corr.matrix,fontsize_row=3,fontsize_col=3, 
                         main = 'Spearman correlation of potential Oct4 targets in serum',
                         breaks = seq(-1,1,0.05), 
                         color = colorRampPalette(rev(RColorBrewer::brewer.pal(11,"RdBu")))(40))
# change Oct4 text color to red
names = oct4_potent_tar$gtable$grobs[[6]]$label
color = rep('black', length(oct4_tar.potent))
color[grep('Pou5f1',names)] <- 'red'
color[grep('Sox2',names)] <- 'red'
color[grep('Nanog',names)] <- 'red'
oct4_potent_tar$gtable$grobs[[5]]$gp=grid::gpar(col=color, fontsize=rep(3,length(oct4_tar.potent)))
oct4_potent_tar$gtable$grobs[[6]]$gp=grid::gpar(col=color, fontsize=rep(3,length(oct4_tar.potent)))
oct4_potent_tar
ggsave('figures/oct4_potent_tar.jpg',oct4_potent_tar, device='jpg', width = 15, height = 15, dpi=400)
rm(names)
rm(color)
rm(tar_length)


# permutation test -> most correlated genes of Oct4 ############################
spearman_corr <- function(sce_obj, clust_num, gene_name, ctrl_genes, clust_name, top_n){
  # gene_name is the name of the core gene used for correlation calculation
  # con_genes means control gene names
  sce <- sce_obj[,sce_obj$label %in% clust_num]
  
  # filter genes that are expressed only in less than 20% cells
  cell_num <- round(ncol(sce)/5,0)
  filtered_genes <- rownames(sce)[nexprs(sce, byrow=TRUE) > cell_num]
  print(paste0('Total number of genes: ',length(filtered_genes)))
  
  # calculate corr of Oct4 and all genes
  get_corr <- function(name1,name2,sce){
    if(name1==name2){
      same <- t(as.matrix(c(1,0,0)))
      colnames(same) <- c('rho','p.value','FDR')
      rownames(same) <- name1
      return(same)
    }
    else{
      r <- as.matrix(correlatePairs(sce,subset.row=c(name1,name2), 
                                    assay.type = "logcounts")[c('rho','p.value','FDR')])
      rownames(r) <- name1
      return(r)
    }
  }
  print(system.time(corr <- mclapply(filtered_genes, FUN = get_corr, 
                                     name2=gene_name, sce=sce, mc.cores=48)))
  corr <- do.call(rbind, corr)
  corr <- as.data.frame(corr)
  
  # filter genes with FDR > 0.01
  corr <- corr[rownames(corr)!=gene_name,]
  corr <- corr[corr$FDR < 0.01,]
  print(paste0('totol number of genes with FDR < 0.01: ', nrow(corr)))
  
  # select positive Rho genes
  corr_posi <- corr[corr$rho >0,]
  corr_posi <- corr_posi[order(corr_posi$FDR, decreasing = FALSE),]
  print(paste0('positive genes number: ', nrow(corr_posi)))
  
  # select negative Rho genes
  corr_nega <- corr[corr$rho <0,]
  corr_nega <- corr_nega[order(corr_nega$FDR, decreasing = FALSE),]
  print(paste0('negative genes number: ', nrow(corr_nega)))
  
  # check control genes' correlations
  ctrl_genes <- ctrl_genes[ctrl_genes != gene_name]
  filtered_ctrl_genes <- ctrl_genes[!ctrl_genes %in% filtered_genes]
  print(paste0('Express < 20% cells control genes : ',paste(filtered_ctrl_genes, collapse = ', ')))
  top_ctrl_genes <- ctrl_genes[ctrl_genes %in% c(rownames(corr_posi[1:top_n,]),rownames(corr_nega[1:top_n,]))]
  print(paste0('Top 30 correlated control genes: ',paste(top_ctrl_genes, collapse = ', ')))
  ctrl_genes <- ctrl_genes[!ctrl_genes %in% c(rownames(corr_posi[1:top_n,]),rownames(corr_nega[1:top_n,]))]
  
  # sort positive control genes
  posi_ctrl <- ctrl_genes[ctrl_genes %in% rownames(corr_posi)]
  posi_order <- order(corr_posi[rownames(corr_posi) %in% posi_ctrl,]$FDR, decreasing = FALSE)
  posi_ctrl <- posi_ctrl[posi_order]
  print(paste0('Positively correlated control genes with FDR <0.01: ',paste(posi_ctrl, collapse = ', ')))
  
  # sort negative control genes
  nega_ctrl <- ctrl_genes[ctrl_genes %in% rownames(corr_nega)]
  nega_order <- order(corr_nega[rownames(corr_nega) %in% nega_ctrl,]$FDR, decreasing = FALSE)
  nega_ctrl <- nega_ctrl[nega_order]
  print(paste0('Negatively correlated control genes with FDR <0.01: ',paste(nega_ctrl, collapse = ', ')))
  
  # sum up and plot gene logcounts
  top_corr_genes <- c(rownames(na.omit(corr_posi[1:top_n,])),
                        gene_name,posi_ctrl,nega_ctrl,
                        rownames(na.omit(corr_nega[1:top_n,])))
  gap1 <- nrow(na.omit(corr_posi[1:top_n,]))
  gap2 <- gap1+length(c(gene_name,posi_ctrl))
  gap3 <- gap2+length(nega_ctrl)
  ifelse(gap3==gap2, gap <- c(gap1,gap2), gap <- c(gap1,gap2,gap3))
  top_corr_plot <- plotHeatmap(sce, exprs_values = 'logcounts', 
                               features = top_corr_genes,
                               columns = names(sort(assays(sce)$logcounts[gene_name,])),
                               cluster_rows = FALSE, cluster_cols=FALSE, 
                               center = TRUE, zlim = c(-5,5),
                               color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(40),
                               main = paste0('Top ',top_n,' positively & negatively correlated genes in ', clust_name),
                               gaps_row = gap,
                               color_columns_by = 'label')
  
  # change the color of target gene
  fontsize = top_corr_plot$gtable$grobs[[3]]$gp$fontsize
  names = top_corr_plot$gtable$grobs[[3]]$label
  color = rep('black', length(top_corr_genes))
  color[grep(paste0('^',gene_name,'$'),names)] <- 'red'
  for(n in c(top_ctrl_genes,posi_ctrl,nega_ctrl)){
    color[grep(paste0('^',n,'$'),names)] <- 'blue'
  }
  top_corr_plot$gtable$grobs[[3]]$gp <- grid::gpar(col=color, fontsize=fontsize)
  
  return(list(corr_df=corr, heatmap=top_corr_plot))
}

fdr_lt_0.01_num <- function(sce_obj,clust_num,gene_name){
  sce <- sce_obj[,sce_obj$label %in% clust_num]
  cell_num <- round(ncol(sce)/5,0)
  filtered_genes <- rownames(sce)[nexprs(sce, byrow=TRUE) > cell_num]
  get_corr <- function(name1,name2,sce){
    if(name1==name2){
      same <- t(as.matrix(c(1,0,0)))
      colnames(same) <- c('rho','p.value','FDR')
      rownames(same) <- name1
      return(same)
    }
    else{
      r <- as.matrix(correlatePairs(sce,subset.row=c(name1,name2), assay.type = "logcounts")[c('rho','p.value','FDR')])
      rownames(r) <- name1
      return(r)
    }
  }
  corr <- mclapply(filtered_genes, FUN = get_corr, name2=gene_name, sce=sce, mc.cores=20)
  corr <- do.call(rbind, corr)
  corr <- as.data.frame(corr)
  corr <- corr[rownames(corr)!=gene_name,]
  corr <- corr[corr$FDR < 0.01,]
  posi_corr <- corr[corr$rho >0,]
  nega_corr <- corr[corr$rho <0,]
  count <- t(as.matrix(c(gene_name,nrow(corr), nrow(posi_corr), nrow(nega_corr))))
  colnames(count) <- c('gene_name','counts','posi', 'nega')
  rownames(count) <- gene_name
  return(count)
}

positive_control_genes <- c("Pou5f1",'Sox2','Nanog','Zfp42','Klf4','Klf2','Klf5',
                            'Esrrb','Eras','Nacc1','Utf1','Lefty1', 'Sall4',
                            'Smad3','Smad1','Gdf3', 'Lin28a','Tfcp2l1','Id3',
                            'Fgf4','Tcl1','Spp1','Upp1','Fbxo15','Dppa3','Trp53',
                            'Dppa5a','Nr0b1','Stat3','Cdyl', 'Mycbp', 'Tbx3', 
                            'Zfx', "Prdm14",'Foxd3', 'Gbx2','Zfp143','Otx2',
                            'Cdx2', "Gata3", "Eomes", 'Tcf15',"Tcf3",'Dnmt3a','Dnmt3b')

# check how many targets (FDR < 0.01) each gene has in serum
print(system.time(fdr_count <- mclapply(unique(c(positive_control_genes,oct4_tar.potent)), 
                                        FUN = fdr_lt_0.01_num, 
                                        sce_obj=sce_dev, clust_num=c(3,5), mc.cores=2)))
fdr_count <- do.call(rbind,fdr_count)
fdr_count <- as.data.frame(fdr_count)
fdr_count <- fdr_count[order(fdr_count$counts, decreasing=TRUE),]
dim(fdr_count)
head(fdr_count, 20)

# cluster3 corr ################################################################
serum_clust3_spear_corr$out <- capture.output(serum_clust3_spear_corr <- 
                                                spearman_corr(sce_dev, c(3),'Pou5f1',
                                                              positive_control_genes , 
                                                              'cluster 3',30))
serum_clust3_spear_corr$out
head(serum_clust3_spear_corr$corr_df)
dim(serum_clust3_spear_corr$corr_df)
ggsave('figures/serum_clust3.jpg',serum_clust3_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# cluster5 corr ################################################################
serum_clust5_spear_corr$out <- capture.output(serum_clust5_spear_corr <- 
                                                spearman_corr(sce_dev, c(5),'Pou5f1',
                                                              positive_control_genes , 
                                                              'cluster 5',30))
serum_clust5_spear_corr$out
head(serum_clust5_spear_corr$corr_df)
dim(serum_clust5_spear_corr$corr_df)
ggsave('figures/serum_clust5.jpg',serum_clust5_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# serum corr ###################################################################
serum_spear_corr$out <- capture.output(serum_spear_corr <- 
                                         spearman_corr(sce_dev, c(3,5),'Pou5f1',
                                                       positive_control_genes , 
                                                       'serum',30))
serum_spear_corr$out
head(serum_spear_corr$corr_df)
dim(serum_spear_corr$corr_df)
ggsave('figures/serum.jpg',serum_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# a2i corr #####################################################################
a2i_spear_corr$out <- capture.output(a2i_spear_corr <- 
                                     spearman_corr(sce_dev, c(2),'Pou5f1',
                                                       positive_control_genes , 
                                                       'a2i',30))
a2i_spear_corr$out
head(a2i_spear_corr$corr_df)
dim(a2i_spear_corr$corr_df)
ggsave('figures/a2i.jpg',a2i_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# cluster 1 corr ###############################################################
clust1_2i_spear_corr$out <- capture.output(clust1_2i_spear_corr <- 
                                           spearman_corr(sce_dev, c(1),'Pou5f1',
                                                       positive_control_genes , 
                                                       'cluster 1',30))
clust1_2i_spear_corr$out
head(clust1_2i_spear_corr$corr_df)
dim(clust1_2i_spear_corr$corr_df)
ggsave('figures/2i_clust1.jpg',clust1_2i_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# cluster 4 corr ###############################################################
clust4_2i_spear_corr$out <- capture.output(clust4_2i_spear_corr <- 
                                           spearman_corr(sce_dev, c(4),'Pou5f1',
                                                       positive_control_genes , 
                                                       'cluster 4',30))
clust4_2i_spear_corr$out
head(clust4_2i_spear_corr$corr_df)
dim(clust4_2i_spear_corr$corr_df)
ggsave('figures/2i_clust4.jpg',clust4_2i_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# 2i corr ######################################################################
i2_spear_corr$out <- capture.output(i2_spear_corr <- 
                                    spearman_corr(sce_dev, c(1,4),'Pou5f1',
                                                       positive_control_genes , 
                                                       '2i',30))
i2_spear_corr$out
head(i2_spear_corr$corr_df)
dim(i2_spear_corr$corr_df)
ggsave('figures/2i.jpg',i2_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# cluster by Oct4 candidate genes ##############################################
sce_seurat <- sce_origin
sce_seurat <- sce_seurat[, !cell_filter$discard]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Pou5f1$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[,assays(sce_seurat)$counts[grepl('^Sox2$',rownames(sce_seurat)),] != 0]
sce_seurat <- sce_seurat[!is.spike,]
sce_seurat <- sce_seurat[rowSums(assays(sce_seurat)$counts > 0) >= 3,]
sce_seurat
sce_seurat <- computeSumFactors(sce_seurat, cluster = quickCluster(sce_seurat))
sce_seurat <- logNormCounts(sce_seurat)
seurat <- as.Seurat(sce_seurat, counts = "counts", data = "logcounts")
seurat <- ScaleData(seurat, features = rownames(seurat))
seurat <- RunPCA(seurat, features = oct4_tar.potent, seed.use = seed)
ElbowPlot(seurat)
VizDimLoadings(seurat, dims = 1:4, reduction = "pca")
PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
PC_num.gml
PC_num.gml <- round(PC_num.gml$dim.est,0)
PC_num.gml
seurat <- RunTSNE(seurat, seed.use = seed, perplexity = 20, dims = 1:PC_num.gml)
DimPlot(seurat, reduction = "tsne" , group.by = 'Cell.Type')
sce_targets <- as.SingleCellExperiment(seurat)
plotReducedDim(sce_targets, "TSNE", colour_by="Cell.Type")  # make sure it's the same as seurat's tSNE
counts(sce_targets) <- as.matrix(counts(sce_targets))
assays(sce_targets)$norm_counts <- 2^(logcounts(sce_targets))-1
sce_targets
sweep_para <- function(seurat, seed, features, res, feature_num){
  para_res <- c()
  sil_score <- c()
  cluster_num <- c()
  ft_num <- c()
  PC_num <- c()
  for(j in feature_num){
    print(paste0('feature number = ', j))
    seurat <- RunPCA(seurat, features = features[1:j], seed.use = seed)
    PC_num.gml <- intrinsicDimension::maxLikGlobalDimEst(seurat[['pca']]@cell.embeddings, k = 10)
    PC_num.gml
    PC_sele <- round(PC_num.gml$dim.est,0)
    seurat <- FindNeighbors(seurat, dims = 1:PC_sele)
    for(i in res){
      print(paste0('res = ', i))
      para_res <- c(para_res,i)
      seurat <- FindClusters(seurat, resolution = i, algorithm = 4)
      cluster <- seurat@meta.data[[paste0('originalexp_snn_res.',i)]]
      sil <- cluster::silhouette(as.integer(cluster), dist(seurat[['pca']]@cell.embeddings[,1:PC_sele]))
      sil.data <- as.data.frame(sil)
      score <- mean(sil.data$sil_width)
      sil_score <- c(sil_score,score)
      cluster_num <- c(cluster_num, length(table(cluster)))
      ft_num <- c(ft_num, j)
      PC_num <- c(PC_num, PC_sele)
    }
    cat('\n')
  }
  r <- data.frame(res=para_res, sil_score=sil_score, clust_num=cluster_num,
                  feature_num=ft_num, PC_num=PC_num)
  return(r)
}
out <- sweep_para(seurat, seed, oct4_tar.potent, 
                  res=seq(0.1,2.0,0.1),
                  feature_num = c(100, 200, length(oct4_tar.potent)))
out
ggplot(out, aes(x=res,y=sil_score,label=paste0('(',res,', ',round(sil_score,3),')'))) +
  geom_line(aes(color = factor(feature_num))) +
  ylab('Silhouette Score') +
  xlab('Leiden Resolution') +
  guides(color=guide_legend(title="Feature Number")) +
  scale_color_brewer(palette="Paired")
ggsave('figures/targets_leiden_para.jpg',device='jpg', width = 8, height = 6)
res_sele <- 1
seurat <- FindNeighbors(seurat, dims = 1:PC_num.gml)
seurat <- FindClusters(seurat, resolution = res_sele, algorithm = 4)  # algorithm 4 is "Leiden"; 1 is "Louvain"
DimPlot(seurat, reduction = "tsne", label = TRUE, shape.by = 'Cell.Type')
ggsave('figures/targets_leiden_res_1.0.jpg',device='jpg', width = 8, height = 6)
colLabels(sce_targets) <- seurat@meta.data[[paste0('originalexp_snn_res.',res_sele)]]
group_marker <- FindAllMarkers(seurat, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
group_marker %>% dplyr::filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  dplyr::slice_min(n = 10, order_by = p_val_adj) -> top10
tar_leiden_markers <- plotHeatmap(sce_targets, exprs_values = "logcounts", 
                              order_columns_by=c("label", "Cell.Type"), 
                              features=top10$gene, cluster_rows = FALSE, center = TRUE,
                              gaps_col = cumsum(as.numeric(table(colLabels(sce_targets)))),
                              gaps_row = seq(0,50,10), zlim= c(-6,6),
                              color = colorRampPalette(rev(RColorBrewer::brewer.pal(10,"RdYlBu")))(30))
ggsave('figures/targets_leiden_markers.jpg', tar_leiden_markers, device='jpg', width = 10, height = 8)
plotReducedDim(sce_targets, "TSNE", shape_by="Cell.Type", colour_by = 'label') 
ggsave('figures/targets_leiden_tsne.jpg', device='jpg', width = 8, height = 7)
sce_targets$dev_clust <- colLabels(sce_dev)
plotReducedDim(sce_targets, "TSNE", shape_by="Cell.Type", colour_by = 'dev_clust')
ggsave('figures/targets_leiden_tsne_dev.jpg', device='jpg', width = 8, height = 7)
table(tar=colLabels(sce_targets), dev=colLabels(sce_dev))
rm(out)
rm(res_sele)
rm(group_marker)
rm(tar_leiden_markers)
rm(PC_num.gml)


# all clusters corr ############################################################
all_spear_corr$out <- capture.output(all_spear_corr <- 
                                      spearman_corr(sce_dev, c(1:5),'Pou5f1',
                                                    positive_control_genes , 
                                                    'all cells',30))
all_spear_corr$out
head(all_spear_corr$corr_df)
dim(all_spear_corr$corr_df)
ggsave('figures/all.jpg',all_spear_corr$heatmap, device='jpg', width = 12, height = 12)


# venn diagram summary ########################################################
ggVennDiagram(list(cluster3_serum = rownames(serum_clust3_spear_corr$corr_df), 
                   cluster5_serum = rownames(serum_clust5_spear_corr$corr_df), 
                   a2i = rownames(a2i_spear_corr$corr_df),
                   cluster1_2i = rownames(clust1_2i_spear_corr$corr_df), 
                   cluster4_2i = rownames(clust4_2i_spear_corr$corr_df)),
                               label = "count", label_size = 5) +
  theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
  ggtitle('intersections of correlated genes in each clusters (Oct4)') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/oct4_venn_1.pdf",device='pdf', width = 12, height = 12)

ggVennDiagram(list(cluster3_serum = rownames(serum_clust3_spear_corr$corr_df), 
                   cluster5_serum = rownames(serum_clust5_spear_corr$corr_df), 
                   a2i = rownames(a2i_spear_corr$corr_df),
                   cluster1_2i = rownames(clust1_2i_spear_corr$corr_df), 
                   cluster4_2i = rownames(clust4_2i_spear_corr$corr_df),
                   candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']),
              label = "count", label_size = 5) +
  theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
  ggtitle('intersections of correlated genes in each clusters (Oct4)') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/oct4_venn_2.pdf",device='pdf', width = 12, height = 12)

ggVennDiagram(list(serum = rownames(serum_spear_corr$corr_df), 
                   a2i = rownames(a2i_spear_corr$corr_df),
                   cluster1_2i = rownames(clust1_2i_spear_corr$corr_df), 
                   cluster4_2i = rownames(clust4_2i_spear_corr$corr_df),
                   candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']),
              label = "count", label_size = 5) +
  theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
  ggtitle('intersections of correlated genes in each clusters (Oct4)') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/oct4_venn_3.pdf",device='pdf', width = 12, height = 12)

# check intersection
find_inter <- function(li){
  r <- list()
  for(i in 2:length(li)){
    arr <- combn(1:length(li),i)
    names <- c()
    for(j in 1:ncol(arr)){
      n <- Reduce(intersect, li[arr[,j]])
      names <- c(names, n)
    }
    r[[i]] <- unique(names)
  }
  return(r)
}
find_inter(list(rownames(serum_spear_corr$corr_df),
                rownames(a2i_spear_corr$corr_df),
                rownames(clust1_2i_spear_corr$corr_df),
                rownames(clust4_2i_spear_corr$corr_df)))

find_inter(list(rownames(serum_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(serum_clust3_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(serum_clust5_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(a2i_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(clust1_2i_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(clust4_2i_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))


# serum Nanog corr #############################################################
# check Nanog variation
nanog_var <- plotExpression(sce_dev,features = c('Nanog','Pou5f1','Sox2',
                                                 'Klf4','Esrrb','Tcl1'), 
                            x='Cell.Type')
ggsave('figures/nanog_variation.jpg',nanog_var, device = 'jpg', width = 10, height = 10)

# plot corr heatmap
serum_nanog_spear_corr$out <- capture.output(serum_nanog_spear_corr <- 
                                       spearman_corr(sce_dev, c(3,5),'Nanog',
                                                     positive_control_genes , 
                                                     'serum (Nanog)',30))
serum_nanog_spear_corr$out
head(serum_nanog_spear_corr$corr_df)
dim(serum_nanog_spear_corr$corr_df)
ggsave('figures/serum_nanog.jpg',serum_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# Serum clust3 Nanog corr #############################################################
serum_clust3_nanog_spear_corr$out <- capture.output(serum_clust3_nanog_spear_corr <- 
                                             spearman_corr(sce_dev, c(3),'Nanog',
                                                           positive_control_genes , 
                                                           'cluster 3 (Nanog)',30))
serum_clust3_nanog_spear_corr$out
head(serum_clust3_nanog_spear_corr$corr_df)
dim(serum_clust3_nanog_spear_corr$corr_df)
ggsave('figures/serum_clust3_nanog.jpg',serum_clust3_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# Serum clust5 Nanog corr #############################################################
serum_clust5_nanog_spear_corr$out <- capture.output(serum_clust5_nanog_spear_corr <- 
                                                      spearman_corr(sce_dev, c(5),'Nanog',
                                                                    positive_control_genes , 
                                                                    'cluster 5 (Nanog)',30))
serum_clust5_nanog_spear_corr$out
head(serum_clust5_nanog_spear_corr$corr_df)
dim(serum_clust5_nanog_spear_corr$corr_df)
ggsave('figures/serum_clust5_nanog.jpg',serum_clust5_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# a2i Nanog corr #############################################################
a2i_nanog_spear_corr$out <- capture.output(a2i_nanog_spear_corr <- 
                                       spearman_corr(sce_dev, c(2),'Nanog',
                                                     positive_control_genes , 
                                                     'a2i (Nanog)',30))
a2i_nanog_spear_corr$out
head(a2i_nanog_spear_corr$corr_df)
dim(a2i_nanog_spear_corr$corr_df)
ggsave('figures/a2i_nanog.jpg',a2i_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# 2i_clust1 Nanog corr #############################################################
clust1_2i_nanog_spear_corr$out <- capture.output(clust1_2i_nanog_spear_corr <- 
                                             spearman_corr(sce_dev, c(1),'Nanog',
                                                           positive_control_genes , 
                                                           'cluster 1 (Nanog)',30))
clust1_2i_nanog_spear_corr$out
head(clust1_2i_nanog_spear_corr$corr_df)
dim(clust1_2i_nanog_spear_corr$corr_df)
ggsave('figures/2i_clust1_nanog.jpg',clust1_2i_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# 2i_clust4 Nanog corr #############################################################
clust4_2i_nanog_spear_corr$out <- capture.output(clust4_2i_nanog_spear_corr <- 
                                                   spearman_corr(sce_dev, c(4),'Nanog',
                                                                 positive_control_genes , 
                                                                 'cluster 4 (Nanog)',30))
clust4_2i_nanog_spear_corr$out
head(clust4_2i_nanog_spear_corr$corr_df)
dim(clust4_2i_nanog_spear_corr$corr_df)
ggsave('figures/2i_clust4_nanog.jpg',clust4_2i_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# 2i_Nanog corr #############################################################
i2_nanog_spear_corr$out <- capture.output(i2_nanog_spear_corr <- 
                                                   spearman_corr(sce_dev, c(1,4),'Nanog',
                                                                 positive_control_genes , 
                                                                 '2i (Nanog)',30))
i2_nanog_spear_corr$out
head(i2_nanog_spear_corr$corr_df)
dim(i2_nanog_spear_corr$corr_df)
ggsave('figures/2i_nanog.jpg',i2_nanog_spear_corr$heatmap, device='jpg', width = 12, height = 12)

# venn diagram summary ########################################################
ggVennDiagram(list(cluster3_serum = rownames(serum_clust3_nanog_spear_corr$corr_df), 
                   cluster5_serum = rownames(serum_clust5_nanog_spear_corr$corr_df),
                   cluster1_2i = rownames(clust1_2i_nanog_spear_corr$corr_df), 
                   cluster4_2i = rownames(clust4_2i_nanog_spear_corr$corr_df)),
              label = "count", label_size = 5) +
  theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
  ggtitle('intersections of correlated genes in each clusters (Nanog)') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/nanog_venn_1.pdf",device='pdf', width = 12, height = 12)

ggVennDiagram(list(cluster3_serum = rownames(serum_clust3_nanog_spear_corr$corr_df), 
                   cluster5_serum = rownames(serum_clust5_nanog_spear_corr$corr_df),
                   a2i = rownames(a2i_nanog_spear_corr$corr_df),
                   cluster1_2i = rownames(clust1_2i_nanog_spear_corr$corr_df), 
                   cluster4_2i = rownames(clust4_2i_nanog_spear_corr$corr_df),
                   candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']),
              label = "count", label_size = 5) +
  theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
  ggtitle('intersections of correlated genes in each clusters (Nanog)') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/nanog_venn_2.pdf",device='pdf', width = 12, height = 12)

ggVennDiagram(list(serum = rownames(serum_nanog_spear_corr$corr_df), 
                   a2i = rownames(a2i_nanog_spear_corr$corr_df),
                   cluster1_2i = rownames(clust1_2i_nanog_spear_corr$corr_df), 
                   cluster4_2i = rownames(clust4_2i_nanog_spear_corr$corr_df),
                   candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']),
              label = "count", label_size = 5) +
  theme(legend.position = "none", plot.margin=margin(1,1,1,1,'cm')) +
  ggtitle('intersections of correlated genes in each clusters (Nanog)') +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("figures/nanog_venn_3.pdf",device='pdf', width = 12, height = 12)

# check intersection
find_inter(list(serum = rownames(serum_nanog_spear_corr$corr_df), 
                a2i = rownames(a2i_nanog_spear_corr$corr_df), 
                cluster1_2i = rownames(clust1_2i_nanog_spear_corr$corr_df), 
                cluster4_2i = rownames(clust4_2i_nanog_spear_corr$corr_df)))

find_inter(list(rownames(serum_nanog_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(serum_clust3_nanog_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(serum_clust5_nanog_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(a2i_nanog_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(clust1_2i_nanog_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
find_inter(list(rownames(clust4_2i_nanog_spear_corr$corr_df),
                candidates = oct4_tar.potent[oct4_tar.potent != 'Pou5f1']))
