library(DESeq2)
library(IHW)

#load read count table
df_te <- read.table('No_ChrE/ReadCounts_NoChrE.csv', sep = ',',header = T) # dir = path where df_te is stored
rownames(df_te) <- df_te[,1]
df_te[,1]<-NULL

###call design matrix
colData = read.table('/Users/akankshitadash/Desktop/DESEQ/de_te_design.txt',sep='\t',header=TRUE) # dir = path where de_te_design is stored

###deseq analysis
###use a design formula that models the assay different at t0, the diff over time,
###and any assay-specific differences over time (the interaction term assay:time)
design = ~ assay + time + assay:time
dds <- DESeqDataSetFromMatrix(countData = df_te, colData = colData,design = design)

# ###run analysis to produce genes which at one or more time points after time 0 showed a TE effect
# reduced = ~ assay + time
# dds <- DESeq(dds, test='LRT', reduced = reduced)

#run DESeq
dds <- DESeq(dds)

alpha <- 0.05
res <- results(dds,alpha = alpha)
#res <- results(dds,alpha = alpha, filterFun = IHW)
resultsNames(dds)

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = list(factors),test = "Wald",alpha = alpha)
  
  summary(res_df)
  res_df <- na.omit(res_df)
  # #
  # res_df$SYMBOL <- mapIds(org.Mm.eg.db,
  #                         keys = as.character(rownames(res_df)),
  #                         column = 'SYMBOL',
  #                         keytype = 'ACCNUM',
  #                         multiVals = 'first')
  # ###remove statistically insignificant genes
  # res_df <- res_df[res_df$padj<alpha,]
  # res_df <- data.frame(ACCNUM=row.names(data.frame(res_df)), res_df,row.names = NULL)
  # # #
  # ###select columns Accnum and Log2FC
  # res_df <- res_df[,c("ACCNUM","log2FoldChange")]
  # res_df <- data.frame(res_df, row.names=NULL)
  # # 
  # # ###select AccNum, log2FC and padj
  # # res_df <- data.frame(ACCNUM = row.names(data.frame(res_df)), res_df, row.names = NULL)
  # # res_df <- res_df[,c("ACCNUM","log2FoldChange","padj")]
  
  ###add symbol and select padj
  # res_df <- data.frame(res_df)
  # rownames(res_df) <- NULL
  # res_df <- res_df[,c("SYMBOL","log2FoldChange")]
}

###deseq at individual time points for RNA
res_1_rna <- filter_sig_res(dds, c("time_t1_vs_t0"))
res_2_rna <- filter_sig_res(dds, c("time_t2_vs_t0"))
res_3_rna <- filter_sig_res(dds, c("time_t3_vs_t0"))
res_4_rna <- filter_sig_res(dds, c("time_t4_vs_t0"))

###deseq at individual time points for RPF
res_1_rpf <- filter_sig_res(dds, c("time_t1_vs_t0","assayrpf.timet1"))
res_2_rpf <- filter_sig_res(dds, c("time_t2_vs_t0","assayrpf.timet2"))
res_3_rpf <- filter_sig_res(dds, c("time_t3_vs_t0","assayrpf.timet3"))
res_4_rpf <- filter_sig_res(dds, c("time_t4_vs_t0","assayrpf.timet4"))

###deseq at individual time points for TE, controlling for the baseline
res_1_te_0 <- filter_sig_res(dds, c("assayrpf.timet1"))
res_2_te_0 <- filter_sig_res(dds, c("assayrpf.timet2"))
res_3_te_0 <- filter_sig_res(dds, c("assayrpf.timet3"))
res_4_te_0 <- filter_sig_res(dds, c("assayrpf.timet4"))

###deseq at individual time points for TE, true is using LRT, not controlling for the baseline
# res_1_te <- filter_sig_res(dds, c("assay_rpf_vs_rna","assayrpf.timet1"))
# res_2_te <- filter_sig_res(dds, c("assay_rpf_vs_rna","assayrpf.timet2"))
# res_3_te <- filter_sig_res(dds, c("assay_rpf_vs_rna","assayrpf.timet3"))
# res_4_te <- filter_sig_res(dds, c("assay_rpf_vs_rna","assayrpf.timet4"))
# res_5_te <- filter_sig_res(dds, c("assay_rpf_vs_rna","assayrpf.timet5"))
# res_6_te <- filter_sig_res(dds, c("assay_rpf_vs_rna","assayrpf.timet6"))

list_all_rna <- list(res_1_rna, res_2_rna, res_3_rna, res_4_rna)
list_all_rpf <- list(res_1_rpf, res_2_rpf, res_3_rpf, res_4_rpf)
list_all_te_0 <- list(res_1_te_0, res_2_te_0, res_3_te_0, res_4_te_0)

###combined results
list_colnames <- c("ACCNUM","lfc_1","lfc_2","lfc_3","lfc_4")
# list_colnames <- c("ACCNUM","lfc_1",'padj_1',"lfc_2",'padj_2',"lfc_3",'padj_3',
#                    "lfc_4",'padj_4',"lfc_5",'padj_5',"lfc_6",'padj_6')
# combined_res <- function(list_res) {
#   res_all <- Reduce(function(x,y) merge(x,y, by = 'AccNum'),list_res)
#   res_all <- setNames(res_all, list_colnames)
# }

# res_all_rna <- combined_res(list_all_rna)
# res_all_rpf <- combined_res(list_all_rpf)
# res_all_te_0 <- combined_res(list_all_te_0)

sessionInfo()

####get normalized matrix
norm_count <- data.frame(counts(dds, normalized = T))
colnames(norm_count) <-c('cdReads0_RNA_1','cdReads1_RNA_1','cdReads2_RNA_1','cdReads3_RNA_1','cdReads4_RNA_1','cdReads0_RNA_2','cdReads1_RNA_2','cdReads2_RNA_2','cdReads3_RNA_2','cdReads4_RNA_2','cdReads0_RPF_1','cdReads1_RPF_1','cdReads2_RPF_1','cdReads3_RPF_1','cdReads4_RPF_1','cdReads0_RPF_2','cdReads1_RPF_2','cdReads2_RPF_2','cdReads3_RPF_2','cdReads4_RPF_2')