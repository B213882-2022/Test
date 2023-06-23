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
library(Seurat)
library(scry)
library(scDblFinder)
library(ggVennDiagram)
library(ggplot2)
library(ggpmisc)
library(ggrepel)
library(ggpubr)

# set working directory
setwd('/home/s2321661/Test')

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

# cell QC
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

# Table summary of QC
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

# Venn Diagram the summarizes the failures in QC
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

# Histogram summary of QC
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

# gene filtering
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

# Normalisation
sce <- computeSumFactors(sce, cluster = quickCluster(sce))
sce <- logNormCounts(sce)
sce  # assays(sce) has 'logcounts'

# check oct4 variation
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

# correlation (scatter plot)
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
hk_genes <- c('Actb','Gapdh','Pgk1','Ppia','Rpl38','Rpl13a')
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
rm(hk_genes)

# fearture selection
assays(sce)$norm_counts <- 2^(logcounts(sce))-1
# scran HVG (model gene variance by trend)
#some how modelGeneVarWithSpikes(sce,'spikes') fails...
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
top100_hvg.log <- plotExpression(sce,features_hvg,exprs_values = "logcounts") + ylim(0,20)
top100_hvg.log
ggsave('figures/top100_hvg_log.jpg',top100_hvg.log, device='jpg', width = 20, height = 6)
# normalised raw counts of HVG
top100_hvg.raw <- plotExpression(sce,features_hvg,exprs_values = "norm_counts") + ylim(0,2e+5)
top100_hvg.raw
ggsave('figures/top100_hvg_raw.jpg',top100_hvg.raw, device='jpg', width = 20, height = 6)

# Deviance (select genes that are both highly expressed and highly variable)
dev <- rowData(devianceFeatureSelection(sce,assay = "counts", fam='binomial'))$binomial_deviance
# log counts of deviance
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:100]
top100_dev.log <- plotExpression(sce,features_dev,exprs_values = "logcounts")+ ylim(0,20)
top100_dev.log
ggsave('figures/top100_dev_log.jpg',top100_dev.log, device='jpg', width = 20, height = 6)
# normalised raw counts of deviance
top100_dev.raw <- plotExpression(sce,features_dev,exprs_values = "norm_counts") + ylim(0,2e+5)
top100_dev.raw
ggsave('figures/top100_dev_raw.jpg',top100_dev.raw, device='jpg', width = 20, height = 6)

# check intersection between features selected by HVG and DEV
features_hvg <- getTopHVGs(var.out,n=2000)
features_dev <- names(dev[order(dev,decreasing=TRUE)])[1:2000]
features_venn <- ggVennDiagram(data.frame(features_dev=features_dev,features_hvg=features_hvg)) +
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

# Dimensionality reduction
set.seed(213882)
#set.seed(NULL)

# tSNE perplexity test
# HVG:1000-4000 + perplexity:10,15,20,25,30 (5行4列) 

# different number of features
# by hvg
reduce_dim_hvg.func <- function(sce,num=c(100,250,500,750,1000,1500,2000,3000,4000),ref){
  # num is an array of numbers
  p <- list()
  for(n in num){
    print(n)
    features_hvg <- getTopHVGs(ref,n=n)
    sce_hvg <- runPCA(sce,subset_row = features_hvg)
    sce_hvg <- runTSNE(sce_hvg,perplexity=10, subset_row = features_hvg)
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
results <- reduce_dim_hvg.func(sce,ref=var.out)
dim_reduce_hvg <- ggarrange(plotlist=results[['plots']], 
                            ncol=3, nrow=3, common.legend = TRUE, legend="right", 
                            legend.grob = results[['legend']]) +
                  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
dim_reduce_hvg
ggsave('figures/dim_reduce_hvg.jpg',dim_reduce_hvg, device='jpg', width = 15, height = 35)
rm(results)
# by dev
reduce_dim_dev.func <- function(sce,num=c(100,250,500,750,1000,1500,2000,3000,4000),ref){
  # num is an array of numbers
  p <- list()
  for(n in num){
    print(n)
    features_dev <- names(ref[order(ref,decreasing=TRUE)])[1:n]
    sce_dev <- runPCA(sce,subset_row = features_dev)
    sce_dev <- runTSNE(sce_dev,perplexity=10, subset_row = features_dev)
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
results <- reduce_dim_dev.func(sce,ref=dev)
dim_reduce_dev <- ggarrange(plotlist=results[['plots']], 
                            ncol=3, nrow=3, common.legend = TRUE, legend="right", 
                            legend.grob = results[['legend']]) +
                  theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm"))
dim_reduce_dev
ggsave('figures/dim_reduce_dev.jpg',dim_reduce_dev, device='jpg', width = 15, height = 35)
rm(results)

# Doublet test
dbl.dens <- computeDoubletDensity(sce_dev, subset.row=features_dev, 
                                  d=ncol(reducedDim(sce_dev)))
sce_dev$DoubletScore <- dbl.dens
doublet_tSNE <- plotTSNE(sce_dev, colour_by="DoubletScore")
ggsave('figures/doublet_tSNE.jpg',doublet_tSNE, device='jpg', width = 8, height = 7)
