library(SingleCellExperiment)
library(scater)
library(scran)
library(biomaRt)
library(pheatmap)
library(data.table)
library(gplots)

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

# gene annotation
ah <- AnnotationHub()
#display(ah)  # search 'GRCm39' and 'EnsDb' to find v109 ID
ens.mm.v109 <- AnnotationHub()[["AH109655"]]
columns(ens.mm.v109)
gene_anno <- AnnotationDbi::select(ens.mm.v109, keys=rownames(sce), 
                                   keytype='GENEID', column='GENENAME')
rowData(sce)$GENENAME <- make.names(gene_anno$GENENAME[match(rownames(sce),gene_anno$GENEID)],
                                    unique = TRUE)
rowData(sce)$ENSEMBL <- rownames(sce)
rownames(sce) <- rowData(sce)$GENENAME
count(grepl("^ERCC", rowData(sce)$ENSEMBL))
is.spike <- grepl("^ERCC", rowData(sce)$ENSEMBL)
spike_names <- rowData(sce)$ENSEMBL[grepl("^ERCC", rowData(sce)$ENSEMBL)]
rownames(sce)[is.spike] <- spike_names
sce
head(rowData(sce))
length(grep('^NA',rowData(sce)$GENENAME))  # the number of failures in gene name annotation
length(grep('^NA',rownames(sce)))

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

# get mitochondrial genes from annotation files got from Biomart Website
MT_genes <- read.csv("MT_genes.csv")
head(MT_genes)
MT_id <- paste(MT_genes$Gene.stable.ID, collapse='|')
head(MT_id)

# cell QC
count(grepl(MT_id, rowData(sce)$ENSEMBL))
is.mito <- grepl(MT_id, rowData(sce)$ENSEMBL)
stats <- perCellQCMetrics(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
# or we can store the stats directly to colData 
# by addPerCellQC(sce, subsets=list(ERCC=is.spike, Mt=is.mito))
head(stats)
sce$sum <- stats$sum/1e6
sce$detected <- stats$detected
sce$mito <- stats$subsets_Mt_percent
sce$ERCC <- stats$subsets_ERCC_percent

# check QC failed cells (filter cells)
reasons <- perCellQCFilters(stats,sub.fields=c("subsets_Mt_percent", "subsets_ERCC_percent"))
colSums(as.matrix(reasons))  # check how many cells were dropped 
summary(reasons$discard)
sce_origin <- sce
#sce <- sce[, !reasons$discard]  # drop 40 cells in total
dim(sce)

# Dot plots that summarize cell QC
sce_origin$discard <- reasons$discard
sce_origin$sum <- sce_origin$stats$sum/1e6
sce_origin$detected <- sce_origin$stats$detected
sce_origin$mito <- sce_origin$stats$subsets_Mt_percent
sce_origin$ERCC <- sce_origin$stats$subsets_ERCC_percent
thres_sum <- attributes(reasons$low_lib_size)$thresholds[1]/1e6
thres_detected <- attributes(reasons$low_n_features)$thresholds[1]
thres_mito <- attributes(reasons$high_subsets_Mt_percent)$thresholds[2]
thres_ERCC <- attributes(reasons$high_subsets_ERCC_percent)$thresholds[2]
QC_sum <- gridExtra::grid.arrange(
  plotColData(sce_origin, x="Cell.Type", y="sum", colour_by="discard") +
    ggtitle("Library sizes") + ylab('Total Reads Count') + xlab('Cell Type') + 
    geom_hline(yintercept = thres_sum,
               linetype = "dashed", color = "red") + 
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(labels = c(0,
                                  sapply(sort(seq(5,max(sce_origin$sum),5)), 
                                         function(x){latex2exp::TeX(paste0('$',x,'\\times 10^6$'))})),
                       breaks = sort(c(0,seq(5,max(sce_origin$sum),5)))) +
    annotate('text', x = 3.6, y = thres_sum,
              label = latex2exp::TeX(paste0('$',round(thres_sum,2),'\\times 10^6$')), 
              vjust = 1.2, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce_origin, x="Cell.Type", y="detected", colour_by="discard")+
    ggtitle("Detected genes") + ylab('Number of detected genes') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_detected,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce_origin$detected),2000)) +
    annotate('text',x = 3.5, y = thres_detected,
              label = round(thres_detected,2), 
              vjust = 1.7, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce_origin, x="Cell.Type", y="mito", colour_by="discard")+
    ggtitle("Mito percent") + ylab('Mitochondrial proportion (%)') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_mito,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce_origin$mito),10)) +
    annotate('text',x = 3.5, y = thres_mito,
              label = round(thres_mito,2), 
              vjust = -0.8, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  plotColData(sce_origin, x="Cell.Type", y="ERCC", colour_by="discard")+
    ggtitle("ERCC percent") + ylab('ERCC proportion (%)') +
    xlab('Cell Type') +
    geom_hline(yintercept = thres_ERCC,
               linetype = "dashed", color = "red") +
    theme(plot.title = element_text(size=12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x=element_text(size = 12),
          axis.title.y=element_text(size = 12),
          legend.text=element_text(size=10)) + 
    guides(color=guide_legend("Discard", title.theme = element_text(size = 12))) + 
    scale_y_continuous(breaks=seq(0,max(sce_origin$ERCC),10)) +
    annotate('text',x = 3.5, y = thres_ERCC,
              label = round(thres_ERCC,2), 
              vjust = -0.8, size=3.9, color='red') + 
    coord_cartesian(clip = 'off'),
  ncol=1
)
ggsave('figures/QC_summary.pdf',QC_sum, device='pdf', width = 8, height = 24)

# Table summary of QC
total.cells <- table(colData(sce_origin)$Cell.Type)
keep.cells <- table(colData(sce_origin)[!colData(sce_origin)$discard,]$Cell.Type)
discard.cells <- table(colData(sce_origin)[colData(sce_origin)$discard,]$Cell.Type)
summary.cells <- cbind(total.cells, discard.cells, keep.cells)
summary.cells <- rbind(summary.cells, colSums(summary.cells))
rownames(summary.cells) <- c(rownames(summary.cells)[1:3],'Total')
summary.cells <- as.data.frame(summary.cells)
cell.types <- rownames(summary.cells)
summary.cells <- cbind(cell.types,summary.cells)
colnames(summary.cells) <- c('Cell Types', 'Original', 'Discard', 'Keep')
summary.cells
write.csv(summary.cells ,'table_summary_QC.csv', row.names = FALSE)

# Venn Diagram the summarizes the failures in QC
library(ggVennDiagram)
library(ggplot2)
library(ggpmisc)
venn <- as.data.frame(reasons)
venn_diagram <- ggVennDiagram(apply(venn[1:4], 2, function(x) which(x == TRUE)),
                              label_alpha=0, 
                              set_color = c("deepskyblue2","darkolivegreen","darkorange","darkorchid"), 
                              category.names = c(paste("Library Size (",sum(venn[1]),")",sep = ''),
                                                 paste("Detected Genes (",sum(venn[2]),")",sep = ''),
                                                 paste("Mito Percent (",sum(venn[3]),")",sep = ''),
                                                 paste("ERCC percent (",sum(venn[4]),")",sep = ''))) + 
  scale_fill_gradient(low="white",high = "coral1") +
  ggtitle('Summary of failures in cell QC') + 
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_color_manual(values = c("deepskyblue2","darkolivegreen","darkorange","darkorchid")) +
  annotate('table',label = summary.cells,x = 0.61, y = 0.1,vjust = 0, hjust = -0.3) +
  theme(legend.position = c(0.95, 0.6))
ggsave('venn.jpg',venn_diagram, device='jpg', width = 8, height = 8)

# Histogram summary of QC
hist <- gridExtra::grid.arrange(
  ggplot(as.data.frame(sce$stats), aes(x=sum/1e6))+
    geom_histogram(color='grey80',bins=50) +
    geom_density(alpha=.2, fill="blue") +
    xlab('Library sizes (millions)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_sum), color = 'red', linetype="dashed") +
    annotate('text',x=thres_sum-2,y=100, label=round(thres_sum,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(sce$stats), aes(x=detected))+
    geom_histogram(color='grey80',bins=50) +
    xlab('Number of detected genes') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_detected), color = 'red', linetype="dashed") +
    annotate('text',x=thres_detected-1400,y=60, label=round(thres_detected,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(sce$stats), aes(x=subsets_Mt_percent))+
    geom_histogram(color='grey80',bins=50) +
    xlab('Mitochondrial proportion (%)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_mito), color = 'red', linetype="dashed") +
    annotate('text',x=thres_mito+6,y=200, label=round(thres_mito,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ggplot(as.data.frame(sce$stats), aes(x=subsets_ERCC_percent))+
    geom_histogram(color='grey80',bins=50) +
    xlab('ERCC proportion (%)') +
    ylab('Number of cells') +
    geom_vline(aes(xintercept=thres_ERCC), color = 'red', linetype="dashed") +
    annotate('text',x=thres_ERCC+5,y=300, label=round(thres_ERCC,2),color = 'red') +
    theme(plot.margin=margin(1,1,1,1,'cm')),
  ncol=2
)
ggsave('figures/hist.jpg',hist, device='jpg', width = 10, height = 8)
