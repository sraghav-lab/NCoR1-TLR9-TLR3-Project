load("results/Control_H3K27ac.RData")
load("results/NCoR1_ChIP-Seq.RData")
load("results/Control_RNASeq.RData")

# Extract enhacners associated with CpG, pIC specific and common genes
merge(diff_h3K27ac_2k_summit_peak.bed.annotation.df,H3K27ac_cluster_lrt_df,by.x="PeakID",by.y="genes") %>% 
       filter(SYMBOL %in% CpG_specific_genes) %>% 
       dplyr::select(c(2,3,4,1)) %>% 
       write.table("results/CpG_specific_H3K27ac_peaks.bed",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
merge(diff_h3K27ac_2k_summit_peak.bed.annotation.df,H3K27ac_cluster_lrt_df,by.x="PeakID",by.y="genes") %>% 
       filter(SYMBOL %in% CpG_pIC_common_genes) %>% 
       dplyr::select(c(2,3,4,1)) %>% 
       write.table("results/CpG_pIC_commmon_H3K27ac_peaks.bed",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
merge(diff_h3K27ac_2k_summit_peak.bed.annotation.df,H3K27ac_cluster_lrt_df,by.x="PeakID",by.y="genes") %>% 
       filter(SYMBOL %in% pIC_specific_genes) %>%
       dplyr::select(c(2,3,4,1)) %>% 
       write.table("results/pIC_specific_H3K27ac_peaks.bed",sep = "\t",row.names = FALSE,col.names = FALSE,quote = FALSE)
  
# Gene overlap of H3K27ac cluster and RNAseq cluster

gs.RNASeq = 25000 # Total no. of genes associated with diff H3K27ac peaks
h3k27ac_rnaseq_cluster = newGOM(diff_h3k27ac_cluster_genes.list, 
                                 CpG_pIC_specific_genes,
                                 genome.size=gs.RNASeq)

h3k27ac_rnaseq_cluster_pval = getMatrix(h3k27ac_rnaseq_cluster,name="pval") %>% 
                              as.data.frame() %>% .[c(4,3,2,1),] %>%
                              -log10(.) %>%
                              sapply(.,round,2) %>% 
                              as.data.frame() %>%
                              sapply(.,as.character)

p=drawHeatmap(h3k27ac_rnaseq_cluster,ncolused=5, grid.col="Blues", note.col="black",adj.p = TRUE)
p=t(p$carpet) %>% as.data.frame() 
p=p[c(4,3,2,1),]
colnames(p)= c("Cluster1(pIC)","Cluster2(CpG)")
pdf("plots/H3K27ac_cluster_gene_RNASeq_cluster_overlap.pdf",height = 4,width = 3)
draw(Heatmap(p,name = "Odds Ratio",
        col= colorRamp2(c(0,10,20),c("blue","yellow","red")),
        heatmap_legend_param=list(at=c(0,5,10,15,20),color_bar="continuous", 
                                  legend_direction="horizontal", legend_width=unit(5,"cm"),
                                  title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
        cluster_columns = FALSE,
        cluster_rows = FALSE,
        rect_gp = gpar(col = "grey"),
        column_names_side = "top",
        column_names_rot = 30,
        column_names_gp = gpar(fontsize = 15),
        row_names_gp = gpar(fontsize = 15),
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf(h3k27ac_rnaseq_cluster_pval[i, j]), x, y, gp = gpar(fontsize = 15))}),
        heatmap_legend_side="bottom")
dev.off()

#######################################################################################################################
# H3K27ac density on pIC spefic genes
my_comparisons = list(c("Emp.2hr.CpG","Emp.2hr.pIC"),c("Emp.6hr.CpG","Emp.6hr.pIC"))

H3k27ac_density_at_CpG_pIC_specific_genes = function(df,plotTitle){
  df.melt = melt(df)
  df.melt$variable = factor(df.melt$variable, levels = c("Emp.0hr","Emp.2hr.CpG","Emp.6hr.CpG","Emp.2hr.pIC",    
                                                         "Emp.6hr.pIC","Emp.2hr.CpG.pIC","Emp.6hr.CpG.pIC"))
  df.melt$Stimulation = gsub("Emp.[0-9]hr.|Emp.","",df.melt$variable)
  df.melt$Stimulation = factor(df.melt$Stimulation, levels = c("0hr","CpG","pIC","CpG.pIC"))
  ggplot(df.melt,aes(x=variable,y=value,color=Stimulation))+
    geom_boxplot(lwd=0.8,width=0.5)+
    gg_theme+
    theme(axis.text.x = element_blank())+
    xlab("")+
    ylab("Scaled vst (DESeq2)")+
    scale_color_manual(values = c("#f8766d", "#3eca9a", "#12b5f6", "#ea7ef4"))+
    ggtitle(plotTitle)+
    stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",
                       y_position = 14,method = "wilcox",paired = FALSE,size=4) 
  
}
CpG_specific_genes.peaks= merge(Emp.H3K27ac_vst.df.scaled,H3K27ac_cluster_annotation,,by.y="genes",by.x=0) %>% 
                                    filter(SYMBOL %in% CpG_specific_genes) %>% 
                                    dplyr::select(c(26,2:8))
pIC_specific_genes.peaks= merge(Emp.H3K27ac_vst.df.scaled,H3K27ac_cluster_annotation,,by.y="genes",by.x=0) %>% 
                                    filter(SYMBOL %in% pIC_specific_genes) %>% 
                                    dplyr::select(c(2:8))
CpG_pIC_common_genes.peaks= merge(Emp.H3K27ac_vst.df.scaled,H3K27ac_cluster_annotation,,by.y="genes",by.x=0) %>% 
                                    filter(SYMBOL %in% CpG_pIC_common_genes) %>% 
                                    dplyr::select(c(26,2:8))

p1_CpG= H3k27ac_density_at_CpG_pIC_specific_genes(CpG_specific_genes.peaks,"CpG_specific")
p1_pIC= H3k27ac_density_at_CpG_pIC_specific_genes(pIC_specific_genes.peaks,"pIC_Specific")
p1_CpG_pIC= H3k27ac_density_at_CpG_pIC_specific_genes(CpG_pIC_common_genes.peaks,"Common CpG and pIC")

pdf("plots/H3K27a_at_CpG_pIC_boxplot.pdf",height = 5,width = 20)
grid.arrange(p1_CpG,p1_CpG_pIC,p1_pIC,ncol=3)
dev.off()

###############################################################################################################################
# Correlation of Enhancer activity and gene expression
H3K27ac_RNAseq_cor = merge(H3K27ac_cluster_vst_gene_level,Emp_vst_rep_merged.scaled,by.x="SYMBOL",by.y=0) %>% 
                              filter(SYMBOL %in% unlist(CpG_pIC_specific_genes))
H3K27ac_RNAseq_cor = H3K27ac_RNAseq_cor %>% 
                              mutate(Genes =case_when(SYMBOL %in% CpG_specific_genes ~ "CpG",
                                                      SYMBOL %in% pIC_specific_genes ~ "pIC",
                                                      SYMBOL %in% CpG_pIC_common_genes ~ "CpG_pIC"))

labRow = c("C1qa","C1qb","C1qc","C2","Ccl3","Ccl4","Ccl6","Ccl9","Ccr4","H2Q-10","Il12a",
           "Il1b","Socs3","Il10", # CpG response genes
           "Ccr7","Il12b","Il6","Ebi3","Cd40","Cd44","Il1a", # Common genes
           "Irf7","Il27","Ifit1","Ifit3","Ifit2","Ifnb1","Cxcl9", # pIC specific
           "Cxcl10","Ccl5","Ddx58","Il15","Isg15","Isg20","Stat1")
pdf("plots/H3K27ac_RNASeq_corrplot.pdf",height = 5,width = 19)
grid.arrange(
plot_RNaseq_h3k27ac_scatter(H3K27ac_RNAseq_cor,"Emp_6hr_CpG","Emp.6hr.CpG"),
plot_RNaseq_h3k27ac_scatter(H3K27ac_RNAseq_cor,"Emp_6hr_pIC","Emp.6hr.pIC"),
plot_RNaseq_h3k27ac_scatter(H3K27ac_RNAseq_cor,"Emp_6hr_all3","Emp.6hr.CpG.pIC"),ncol=3)
dev.off()
################################################################################################################################
# H3K27ac enrichment at Synergy genes
my_comparisons = list(c("Emp.2hr.pIC","Emp.6hr.CpG.pIC"),c("Emp.6hr.CpG","Emp.6hr.CpG.pIC"))
synergy_genes.peaks= merge(Emp.H3K27ac_vst.df.scaled,H3K27ac_cluster_annotation,,by.y="genes",by.x=0) %>% 
                          filter(SYMBOL %in% synergy_genes) %>% 
                          dplyr::select(c(26,2:8))
p1_synergy =     H3k27ac_density_at_CpG_pIC_specific_genes(synergy_genes.peaks,"Synergy")

# H3K27ac enrichment at CpG antagonist genes
my_comparisons = list(c("Emp.6hr.CpG","Emp.6hr.CpG.pIC"))
CpG_antagonist_genes.peaks= merge(Emp.H3K27ac_vst.df.scaled,H3K27ac_cluster_annotation,,by.y="genes",by.x=0) %>% 
                                filter(SYMBOL %in% CpG_antagonoist_genes) %>% 
                                dplyr::select(c(26,2:8))

p1_CpG_antagonist = H3k27ac_density_at_CpG_pIC_specific_genes(CpG_antagonist_genes.peaks,"CpG_Antagonist")

# H3K27ac enrichment at pIC antagonist genes
my_comparisons = list(c("Emp.2hr.pIC","Emp.6hr.CpG.pIC"))

pIC_antagonist_genes.peaks= merge(Emp.H3K27ac_vst.df.scaled,
                                  H3K27ac_cluster_annotation,,by.y="genes",by.x=0) %>% 
                              filter(SYMBOL %in% pIC_antagonoist_genes) %>% 
                              dplyr::select(c(26,2:8))
p1_pIC_antagonist =     H3k27ac_density_at_CpG_pIC_specific_genes(pIC_antagonist_genes.peaks,"pIC_antagonist")

pdf('plots/H3K27ac_synergy_antagonist_boxplot.pdf',,height = 5,width = 19)
grid.arrange(p1_synergy,p1_CpG_antagonist,p1_pIC_antagonist,ncol=3)
dev.off()
###########################################################################################################################################
#NCoR1 binding acluster association with CpG and pIC specific genes
diff_NCoR1.annotation.df.list =list()
for (i in unique(diff_NCoR1.annotation.df$Clusters)){
  
  cluster = diff_NCoR1.annotation.df %>%
              filter(Clusters ==i & SYMBOL %in% clusters$df$genes) %>% 
              dplyr::select("SYMBOL") %>% 
              unique(.) %>% 
              pull(.,"SYMBOL")
  diff_NCoR1.annotation.df.list[[i]]= cluster
}

# Gene overlap of NCoR1 cluster and RNAseq cluster
gs.RNASeq = 25000 
NCoR1_Emp_rnaseq_cluster = newGOM(diff_NCoR1.annotation.df.list,
                                  CpG_pIC_specific_genes,
                                  gs.RNASeq)
inter.nl = getNestedList(NCoR1_Emp_rnaseq_cluster,"intersection")

NCoR1_Emp_rnaseq_cluster_pval = getMatrix(NCoR1_Emp_rnaseq_cluster,name="pval") %>% 
                                  as.data.frame() %>%
                                  -log10(.) %>%
                                  sapply(.,round,2) %>% 
                                  as.data.frame() %>% 
                                  sapply(.,as.character)
NCoR1_Emp_rnaseq_cluster_pval = NCoR1_Emp_rnaseq_cluster_pval[nrow(NCoR1_Emp_rnaseq_cluster_pval):1, ]
p=drawHeatmap(NCoR1_Emp_rnaseq_cluster,ncolused=5, grid.col="Blues", note.col="black",adj.p = TRUE)
p=t(p$carpet) %>% as.data.frame() 
#p=p[order(row.names(p)), ]
#colnames(p)= c("Cluster1(pIC)","Cluster2(CpG)")
pdf("plots/NCoR1_diffPeak_CpG_pIC_gene_overlap.pdf", height = 6,width = 5)
draw(Heatmap(p[c(5,4,3,1,2),],name = "Odds Ratio",
             col= colorRamp2(c(0,3,6.5),c("blue","yellow","red")),
             heatmap_legend_param=list(at=c(0,3,6.5),color_bar="continuous", 
                                       legend_direction="horizontal", legend_width=unit(5,"cm"),
                                       title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             rect_gp = gpar(col = "black"),
             column_names_side = "top",
             column_names_rot = 30,
             column_names_gp = gpar(fontsize = 15),
             row_names_gp = gpar(fontsize = 15),
             cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf(NCoR1_Emp_rnaseq_cluster_pval[i, j]), x, y, gp = gpar(fontsize = 15))}),
              heatmap_legend_side="bottom")
dev.off()
############################################################################################################################
# NFKB, IRF3 and PU.1 motif density @ Diff NCoR1 binding sites.
motif.density  = list.files(path = "Diff_NCoR1_associated_genes/",pattern = "H3K27ac_at_2fold_*")
motif.density = motif.density[c(1,5,3,2,4)]
motif.density.plot = list()
for (i in motif.density) {
  ptitle = gsub(".motif_density.txt|H3K27ac_at_","",i)
  motif.den = read.csv(paste0("Diff_NCoR1_associated_genes/",i),sep = "\t",header = T) %>% dplyr::select(c(1,2,5,8)) 
  colnames(motif.den) = c("Distance","IRF3","NFkB","Spib(PU.1)")
  motif.den.melt = melt(motif.den,id.vars = "Distance",variable.name = "TF Motifs")
  motif.density.plot[[i]] = motif.den.melt %>% ggplot(aes(x=Distance,y=value,color=`TF Motifs`))+geom_line(size=1)+gg_theme+
                            #ggtitle(ptitle)+
                            labs(
                            x="Distance From NCoR1 Center", y = "Motifs per peak per base")

}
do.call(grid.arrange, c(motif.density.plot,ncol=1))


#########################################################################################################################################################

