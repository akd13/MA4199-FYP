###Volcano Plot
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)
library(AnnotationDbi)

###res_tableDE: object from running results in DESeq
###type: rnabulk/rnacyt/..etc 
###list_of_genes: particular group of genes to be highlighted, can leave it as c() if none highlighted
###list_of_gtpasec: particular group of genes to be highlighted, can leave it as c() if none highlighted
###label: label of particular group of genes to be highlighted, can leave it as '' if none highlighted
###folder_label: name of folder to save plots
###xlims: c(x_min, x_max); ylims: c(y_min, y_max)
plot_volcano <- function(res_tableDE, type, day,list_of_genes,label,folder_label,xlims,ylims){
  
  res_tableDE <- data.frame(res_tableDE)
  res_tableDE$SYMBOL <- mapIds(org.Hs.eg.db,
                               keys = as.character(rownames(res_tableDE)),
                               column = 'SYMBOL',
                               keytype ='ACCNUM',
                               multiVals = 'first')
  threshold_DE <- res_tableDE$padj < 0.1
  #res_tableDE[1,]$SYMBOL<-'chrE'
  res_tableDE$threshold <- threshold_DE
  #res_tableDE$threshold[res_tableDE$SYMBOL %in% list_of_genes] <- 'highlight'
  list_of_gtpase <-c('ITPR3','ZFHX3','TNRC18','USP49','LRFN4',
                  'NBEAL2','ZNF704','PLEKHM2','ZNF462','PPP6R2',
                  'MEGF8','KIAA0100','ABL1','FAM193A','CSNK1G2',
                  'TRAK1','MYO10','CACNA1G','MFHAS1','RGMB','AXIN1',
                  'DOT1L','KMT2C','USP36','PPP1R16A','TRIM71','TFDP1',
                  'MED13L','PPP1R37','ATP10A','TBC1D9B','PITX1','SGSM2',
                  'KAT6A','PCDH17','NSD1','TRAPPC10','KCNC1','B4GALNT4',
                  'SLX4','PLXNA3','CIT','RNF10','GBA2','PIGH','BRD1','AGAP2',
                  'ADAMTS5','GPR153','FBF1','INPP5E','COL27A1','SKI','NXPE3',
                  'TET3','CHD2','ATP8B2','PKD1','DLG5','TRIM8','NSMF','SEL1L3','BTAF1',
                  'MBTPS1','MAST2','ATXN1','PTPRS','AADAT','SHB','MAN1A2','MYO18B','POMGNT2',
                  'SH2B3','GLI3','MAPKAPK5','ATP11A','TNRC6C','CAPN15','SMO','PCNX3','CHD6','DAPK1',
                  'ALKBH5','SLC35B2','RAPGEF1','KIFC2','SSH1','ZC3H11A','SPRED1','ZNF324B','UTS2R',
                  'DENND4B','KIF26A','ADCY3','KIF1C','ARHGAP39','SBF1','ATXN7L3','FOXN3','LRP8','MAP3K11','PER3','GABBR1','EXT1','RALGDS','HS6ST1','PLXNA1','UBN1','AGBL5','ZFAT','ARHGEF12','TRIO','TTBK2','RBBP8','KIF21B','ARHGAP23','PITPNM2','FOXM1','ZNF71','EFL1','KDM4B','ITGB8','HIPK2','HACE1','ZNF496','ABTB2','TBC1D9','TAB3','CACNA1H','ERCC8','HIRA','DAB2IP','LRRC8B','HIVEP2','IL11','NF2','MARK3','MLLT6','MPRIP','ZNF236','IRS1','ZNF835','PPFIA4','PAPD7','CREBBP','SLC40A1','CDC42BPB','PTPN21','PDE8A','CIC','B4GALNT3','UBE3C','CAMSAP1','EPB41L1','STK11','ST3GAL2','DZIP1','IPO13','JARID2','VPS54','HIVEP3','BRSK1','SLF2','TIAM1','PTPRG','REXO1','B4GALT2','TNK2','TOP3B','RAB11FIP3','GRIN3A','RC3H2','LRP3','GNAZ','CELSR3','FBXW7','CCDC88C','USP31','TMC6','NDST1','PCF11','TULP4','SH3RF3','LATS2','SLC43A1','CUEDC1','TRAK2','PACS2','HOMER2','FBXO41','SHANK3','AZIN1','MED13','FAM20C','SETD1A','ZDHHC6','BRD4','GPBP1L1','SETD5','MTCL1','SIPA1L2','STARD9','SHD','TWNK','SEMA6B','SPTBN2','EP300','CRAMP1','PREX1','PHF2','ADGRL1','SOBP','TRAF3','ARHGEF11',
                  'DOCK5','SCN8A','GLI2','NIPA2','MYO9B','TNRC6A','AGAP1','KMT2B')
  
  res_tableDE$genelabels <- ""
  res_tableDE$genelabels[res_tableDE$SYMBOL %in% list_of_genes] <- label

  ggplot(res_tableDE,aes(x = log2FoldChange, y = -log10(padj))) +
    
    ylim(ylims[[1]],ylims[[2]]) +
    xlim(xlims[[1]], xlims[[2]]) +
    
    #threshold lines
    geom_vline(xintercept = -1, linetype = 'dashed', color = 'brown',size = 0.7) +
    geom_vline(xintercept = 1, linetype = 'dashed', color = 'brown', size = 0.7) +
    geom_hline(yintercept = 1, linetype = 'dashed', color = 'brown', size = 0.7) +

  #points
  geom_point(aes(colour = threshold),alpha = 0.5) +
    geom_point(data = subset(res_tableDE, SYMBOL %in% list_of_gtpase),colour = 'purple')+
    geom_point(data = subset(res_tableDE, SYMBOL %in% list_of_genes),
               aes(fill = genelabels),colour = 'blue')+
    #scale_fill_manual(name = 'Legend',
     #                 values = cols)+
    guides(colour = F)+
    
    ###add gene labels for particular group from list_of_genes
    geom_text_repel(data = subset(res_tableDE, SYMBOL %in% list_of_genes),# & log2FoldChange > 1),
                    #data = subset(res_tableDE, SYMBOL %in% c('Htra2','Tomm7','Sqstm1')), #for rna autophagy
                    aes(label = SYMBOL), size = 3,
                    box.padding = unit(0.4, 'lines'),
                    point.padding = unit(0.4, 'lines'),
                    segment.size = 0.2, segment.colour = 'grey50')+
    
    theme_bw() + 
    
    ggtitle(paste(type,day,sep = '')) +
    xlab(bquote(~log[2]~ "FC")) +
    ylab(bquote(~-log[10]~italic(p-adj)))+
    
    theme(legend.position = c(0.15,0.9),
          #legend.position = 'none',
          legend.title = element_blank(),
          legend.background = element_rect(color = 'black',size = 0.5,linetype= 'solid'),
          legend.text = element_text(size = 12),
          plot.title = element_text(size = rel(1.5),hjust = 0.5),
          axis.title = element_text(size = rel(1.25)),
          axis.text.y = element_text(size = 12),
          axis.text.x = element_text(size=14))
  save_dir <- '/Users/akankshitadash/Desktop/Programs/MA4199/DESEQ/No_ChrE/VolcanoPlots/'
  filepath <- paste(save_dir,folder_label,'/',type,toString(day),'.png',sep = '')
  ggsave(filepath,width = 9, height = 7)
}
