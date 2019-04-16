library(DESeq2)

###load read count table
df_te_enriched <- read.table('/Users/akankshitadash/Desktop/DESEQ/',sep = '\t',header = T) # dir = path where df_te_enriched is stored
rownames(df_te_enriched) <- df_te_enriched[,1]
df_te_enriched[,1] <- NULL

###call design matrix
colData = read.table('/Users/akankshitadash/Desktop/DESEQ/de_te_enriched_design.txt',sep='\t',header=TRUE) # dir = path where de_te_enriched_design is stored

###deseq analysis
###setup DESeq matrix
#design doesn't normalize to t0, if normalize to t0, use design = ~time + time:assay + time:cond + assay:cond + time:cond:assay 
design = ~time + time:assay + time:cond + time:cond:assay 
dds <- DESeqDataSetFromMatrix(countData = df_te_enriched, colData = colData,design = design)

###run DESeq
dds <- DESeq(dds)

alpha <- 0.1
res <- results(dds,alpha = alpha)
resultsNames(dds)

###function to compare different combinations of factors
filter_sig_res <- function(dds, factors) {
  res_df <- results(dds, contrast = list(factors),test = "Wald",alpha = alpha)
  
  summary(res_df)
  res_df <- na.omit(res_df)
  # # #
  # res_df$SYMBOL <- mapIds(org.Mm.eg.db,
  #                         keys = as.character(rownames(res_df)),
  #                         column = 'SYMBOL',
  #                         keytype = 'ACCNUM',
  #                         multiVals = 'first')
  # res_df <- data.frame(res_df)
  # res_df <- res_df[,c("SYMBOL",'log2FoldChange')]
  
  # ###remove statistically insignificant genes
  # res_df <- res_df[res_df$padj < alpha,]
  # res_df <- data.frame(ACCNUM=row.names(data.frame(res_df)), res_df,row.names = NULL)
  # #
  # ###select columns Accnum and Log2FC
  # res_df <- res_df[,c("ACCNUM","log2FoldChange")]
  # res_df <- data.frame(res_df, row.names=NULL)
  # 
  # ###select AccNum, log2FC and padj
  # res_df <- data.frame(ACCNUM = row.names(data.frame(res_df)), res_df, row.names = NULL)
  # res_df <- res_df[,c("ACCNUM","log2FoldChange","padj")]
  
  ###add symbol and select padj
  # res_df <- data.frame(res_df)
  # rownames(res_df) <- NULL
  # res_df <- res_df[,c("SYMBOL","log2FoldChange","padj")]

}

###DESeq at individual time points for RNA Cyto
# res_1_rna_cyt <- filter_sig_res(dds, c("time_t1_vs_t0"))
# res_2_rna_cyt <- filter_sig_res(dds, c("time_t2_vs_t0"))
# res_3_rna_cyt <- filter_sig_res(dds, c("time_t3_vs_t0"))
# res_4_rna_cyt <- filter_sig_res(dds, c("time_t4_vs_t0"))
# res_5_rna_cyt <- filter_sig_res(dds, c("time_t5_vs_t0"))
# res_6_rna_cyt <- filter_sig_res(dds, c("time_t6_vs_t0"))

###DESeq at individual time points for RNA Crude
# res_1_rna_cru <- filter_sig_res(dds, c("time_t1_vs_t0", 'timet1.condb'))
# res_2_rna_cru <- filter_sig_res(dds, c("time_t2_vs_t0", 'timet2.condb'))
# res_3_rna_cru <- filter_sig_res(dds, c("time_t3_vs_t0", 'timet3.condb'))
# res_4_rna_cru <- filter_sig_res(dds, c("time_t4_vs_t0", 'timet4.condb'))
# res_5_rna_cru <- filter_sig_res(dds, c("time_t5_vs_t0", 'timet5.condb'))
# res_6_rna_cru <- filter_sig_res(dds, c("time_t6_vs_t0", 'timet6.condb'))

###DESeq at individual time points for RNA enrichment
res_0_rna_enrich <- filter_sig_res(dds, c('timet0.condb'))
res_1_rna_enrich <- filter_sig_res(dds, c("timet1.condb"))
res_2_rna_enrich <- filter_sig_res(dds, c("timet2.condb"))
res_3_rna_enrich <- filter_sig_res(dds, c("timet3.condb"))
res_4_rna_enrich <- filter_sig_res(dds, c("timet4.condb"))
res_5_rna_enrich <- filter_sig_res(dds, c("timet5.condb"))
res_6_rna_enrich <- filter_sig_res(dds, c("timet6.condb"))

###DESeq at individual time points for RPF Cyto
# res_1_rpf_cyt <- filter_sig_res(dds, c("time_t1_vs_t0", "timet1.assayrpf"))
# res_2_rpf_cyt <- filter_sig_res(dds, c("time_t2_vs_t0", "timet2.assayrpf"))
# res_3_rpf_cyt <- filter_sig_res(dds, c("time_t3_vs_t0", "timet3.assayrpf"))
# res_4_rpf_cyt <- filter_sig_res(dds, c("time_t4_vs_t0", "timet4.assayrpf"))
# res_5_rpf_cyt <- filter_sig_res(dds, c("time_t5_vs_t0", "timet5.assayrpf"))
# res_6_rpf_cyt <- filter_sig_res(dds, c("time_t6_vs_t0", "timet6.assayrpf"))

###DESeq at individual time points for RPF Crude
# res_1_rpf_cru <- filter_sig_res(dds, c("time_t1_vs_t0", "timet1.assayrpf", "timet1.condb","timet1.assayrpf.condb"))
# res_2_rpf_cru <- filter_sig_res(dds, c("time_t2_vs_t0", "timet2.assayrpf", "timet2.condb","timet2.assayrpf.condb"))
# res_3_rpf_cru <- filter_sig_res(dds, c("time_t3_vs_t0", "timet3.assayrpf", "timet3.condb","timet3.assayrpf.condb"))
# res_4_rpf_cru <- filter_sig_res(dds, c("time_t4_vs_t0", "timet4.assayrpf", "timet4.condb","timet4.assayrpf.condb"))
# res_5_rpf_cru <- filter_sig_res(dds, c("time_t5_vs_t0", "timet5.assayrpf", "timet5.condb","timet5.assayrpf.condb"))
# res_6_rpf_cru <- filter_sig_res(dds, c("time_t6_vs_t0", "timet6.assayrpf", "timet6.condb","timet6.assayrpf.condb"))

###DESeq at individual time points for RPF enrichment
res_0_rpf_enrich <- filter_sig_res(dds, c("timet0.condb","timet0.assayrpf.condb"))
res_1_rpf_enrich <- filter_sig_res(dds, c("timet1.condb","timet1.assayrpf.condb"))
res_2_rpf_enrich <- filter_sig_res(dds, c("timet2.condb","timet2.assayrpf.condb"))
res_3_rpf_enrich <- filter_sig_res(dds, c("timet3.condb","timet3.assayrpf.condb"))
res_4_rpf_enrich <- filter_sig_res(dds, c("timet4.condb","timet4.assayrpf.condb"))
res_5_rpf_enrich <- filter_sig_res(dds, c("timet5.condb","timet5.assayrpf.condb"))
res_6_rpf_enrich <- filter_sig_res(dds, c("timet6.condb","timet6.assayrpf.condb"))

###DESeq at individual time points for TE Cyto
# res_1_te_cyt <- filter_sig_res(dds, c("timet1.assayrpf"))
# res_2_te_cyt <- filter_sig_res(dds, c("timet2.assayrpf"))
# res_3_te_cyt <- filter_sig_res(dds, c("timet3.assayrpf"))
# res_4_te_cyt <- filter_sig_res(dds, c("timet4.assayrpf"))
# res_5_te_cyt <- filter_sig_res(dds, c("timet5.assayrpf"))
# res_6_te_cyt <- filter_sig_res(dds, c("timet6.assayrpf"))

###DESeq at individual time points for TE Crude
# res_1_te_cru <- filter_sig_res(dds, c("timet1.assayrpf", "timet1.assayrpf.condb"))
# res_2_te_cru <- filter_sig_res(dds, c("timet2.assayrpf", "timet2.assayrpf.condb"))
# res_3_te_cru <- filter_sig_res(dds, c("timet3.assayrpf", "timet3.assayrpf.condb"))
# res_4_te_cru <- filter_sig_res(dds, c("timet4.assayrpf", "timet4.assayrpf.condb"))
# res_5_te_cru <- filter_sig_res(dds, c("timet5.assayrpf", "timet5.assayrpf.condb"))
# res_6_te_cru <- filter_sig_res(dds, c("timet6.assayrpf", "timet6.assayrpf.condb"))

###DESeq at individual time points for TE enrichment
res_0_te_enrich <- filter_sig_res(dds, c('timet0.assayrpf.condb'))
res_1_te_enrich <- filter_sig_res(dds, c("timet1.assayrpf.condb"))
res_2_te_enrich <- filter_sig_res(dds, c("timet2.assayrpf.condb"))
res_3_te_enrich <- filter_sig_res(dds, c("timet3.assayrpf.condb"))
res_4_te_enrich <- filter_sig_res(dds, c("timet4.assayrpf.condb"))
res_5_te_enrich <- filter_sig_res(dds, c("timet5.assayrpf.condb"))
res_6_te_enrich <- filter_sig_res(dds, c("timet6.assayrpf.condb"))

list_all_rna_enrich <- list(res_0_rna_enrich,res_1_rna_enrich, res_2_rna_enrich, res_3_rna_enrich,
                            res_4_rna_enrich, res_5_rna_enrich, res_6_rna_enrich)

list_all_rpf_enrich <- list(res_0_rpf_enrich,res_1_rpf_enrich, res_2_rpf_enrich, res_3_rpf_enrich,
                            res_4_rpf_enrich, res_5_rpf_enrich, res_6_rpf_enrich)

list_all_te_enrich <- list(res_0_te_enrich,res_1_te_enrich, res_2_te_enrich, res_3_te_enrich,
                           res_4_te_enrich, res_5_te_enrich, res_6_te_enrich)

sessionInfo()