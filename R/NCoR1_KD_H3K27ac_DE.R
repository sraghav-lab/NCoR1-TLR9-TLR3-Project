
# Perform DESeq2 to identify significant differentially acetylated regions between NCoR1 KD vs Control condition
consensusToCount = import.bed("../../NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered.bed")
H3K27ac = read.csv("../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/Emp_NCoR1_all_condition_merged_summits_2kb.annotation",
                   sep = "\t",header = T)[,c(1,2,3,4,5)] 
colnames(H3K27ac)[1] = "GeneID"
H3K27ac = H3K27ac %>% filter(GeneID %in% Total_diff_peaks_merged_filtered$PeakID)
bamsToCount_Rep1 <- dir("../../NCOR1/NCoR1_H3K27ac_analysis/bowtie_out/", full.names = TRUE, pattern = "*.25M.bam")
bamsToCount_Rep2 <- dir("../../NCOR1/NCoR1_H3K27ac_Rep2_analysis/bowtie_out/", full.names = TRUE, pattern = "*.dup.filtered.srt.q10.bam$")
#bamsToCount <- subset(bamsToCount,!grepl("*_2h*", bamsToCount))
fcResults <- featureCounts(c(bamsToCount_Rep1,bamsToCount_Rep2), 
                            annot.ext = H3K27ac, isPairedEnd = FALSE,
                            countMultiMappingReads = FALSE, ,nthreads = 15)
myCounts <- fcResults$counts

colnames(myCounts)[2] ="Emp.2hr.CpG.pIC.25M.bam"
colnames(myCounts)[5] = "Emp.6hr.CpG.pIC.25M.bam"
colnames(myCounts)[9] ="NCoR1.2hr.CpG.pIC.25M.bam"
colnames(myCounts)[12] = "NCoR1.6hr.CpG.pIC.25M.bam"

colnames(myCounts) = gsub(".25M.bam|.Rep2.dup.filtered.srt.q10.bam","",colnames(myCounts))
colnames(myCounts)[c(15,24)] = c("NCoR1.2hr.CpG","Emp.0hr")
metaData <- data.frame(Group=colnames(myCounts))
metaData$Rep = factor(c(rep(1,14),rep(2,14)))
metaData$Stimulation = factor(gsub("Emp.|NCoR1.","",metaData$Group))
metaData$Genotype = factor(gsub("\\.[0-9].*","",metaData$Group))
metaData$Sample = paste0(metaData$Group,"_",metaData$Rep)
rownames(metaData)= metaData$Sample
colnames(myCounts) =metaData$Sample

metaData = metaData[order(metaData$Group),]

myCounts = myCounts[,metaData$Sample]

H3K27ac_DDS <- DESeqDataSetFromMatrix(myCounts, colData = metaData, design=~Group, rowRanges = consensusToCount)
H3K27ac_DDS <- DESeq(H3K27ac_DDS)

H3K27ac_Rlog <- vst(H3K27ac_DDS)
H3K27ac_Rlog.df = as.data.frame(assay(H3K27ac_Rlog))
H3K27ac_Rlog_rep_merged.df = H3K27ac_Rlog.df %>% merge_rep(.)
colnames(H3K27ac_Rlog_rep_merged.df) = gsub("_1","",colnames(H3K27ac_Rlog_rep_merged.df))

#plotPCA(H3K27ac_Rlog, intgroup = c("Genotype","Stimulation"))
distsRL <- dist(t(assay(H3K27ac_Rlog[,c(1:ncol(H3K27ac_Rlog))])))
H3K27ac_mat<- as.matrix(distsRL)
rownames(H3K27ac_mat) <- colnames(H3K27ac_mat) <- paste(metaData[c(1:28),]$Group,":", metaData[c(1:28),]$Rep)
                                       paste(Group,Rep, sep=' : ')))
hc <- hclust(distsRL)

pdf("plots/H3K27ac_NCoR1_KD_sample_hclust.pdf",height = 6,width = 7.5)
Heatmap(H3K27ac_mat,col = my_palette,name = "Euclidean Distance")
dev.off()


#######################################################################################################################
# 
h3k27ac_exp_plotlist =list()
# NCoR1 2hr CpG vs Emp 6hr CpG
NCoR1.2hrCpG_vs_Emp.2hrCpG = results(H3K27ac_DDS,contrast=c("Group","NCoR1.2hr.CpG","Emp.2hr.CpG")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.2hrCpG_vs_Emp.2hrCpG,"NCoR1.2hr.CpG vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.2hrCpG.png",dpi=1200,height = 5,width = 5.7)

# NCoR1 6hr CpG vs Emp 6hr CpG
NCoR1.6hrCpG_vs_Emp.6hrCpG = results(H3K27ac_DDS,contrast=c("Group","NCoR1.6hr.CpG","Emp.6hr.CpG")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.6hrCpG_vs_Emp.6hrCpG,"NCoR1.6hr.CpG vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.6hrCpG.png",dpi=1200,height = 5,width = 5.7)

# NCoR1 2hr pIC vs Emp 2hr pIC
NCoR1.2hrpIC_vs_Emp.2hrpIC = results(H3K27ac_DDS,contrast=c("Group","NCoR1.2hr.pIC","Emp.2hr.pIC")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.2hrpIC_vs_Emp.2hrpIC,"NCoR1.2hr.pIC vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.2hrpIC.png",dpi=1200,height = 5,width = 5.7)

# NCoR1 6hr pIC vs Emp 6hr pIC
NCoR1.6hrpIC_vs_Emp.6hrpIC = results(H3K27ac_DDS,contrast=c("Group","NCoR1.6hr.pIC","Emp.6hr.pIC")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.6hrpIC_vs_Emp.6hrpIC,"NCoR1.6hr.pIC vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.6hrpIC.png",dpi=1200,height = 5,width = 5.7)

# NCoR1 6hr CpG+pIC vs Emp 6hr CpG+pIC
NCoR1.2hrCpG.pIC_vs_Emp.2hrCpG.pIC = results(H3K27ac_DDS,contrast=c("Group","NCoR1.2hr.CpG.pIC","Emp.2hr.CpG.pIC")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.2hrCpG.pIC_vs_Emp.2hrCpG.pIC,"NCoR1.2hr.CpG.pIC vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.2hrCpGpIC.png",1200=800,height = 5,width = 5.7)

# NCoR1 6hr CpG+pIC vs Emp 6hr CpG+pIC
NCoR1.6hrCpG.pIC_vs_Emp.6hrCpG.pIC = results(H3K27ac_DDS,contrast=c("Group","NCoR1.6hr.CpG.pIC","Emp.6hr.CpG.pIC")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.6hrCpG.pIC_vs_Emp.6hrCpG.pIC,"NCoR1.6hr.CpG.pIC vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.6hrCpGpIC.png",dpi=1200,height = 5,width = 5.7)

NCoR1.0hr_vs_Emp.0hr = results(H3K27ac_DDS,contrast=c("Group","NCoR1.0hr","Emp.0hr")) %>% as.data.frame() 
KD.H3K27ac.Volcano.plot(NCoR1.0hr_vs_Emp.0hr,"NCoR1.0hr vs Emp")
ggsave("plots/Diff_H3K27ac_scatter_plot_NCoR1.0hr.png",dpi=1200,height = 5,width = 5.7)

#pdf("plots/Diff_H3K27ac_scatter_plot_Emp_NCoR1.pdf",height = 3,width = 5)
#do.call("grid.arrange", c(h3k27ac_exp_plotlist,ncol=2,nrow=2))
#dev.off()
#####################################################################################
H3K27ac_NCoR1_KD_DE_PeakID = lapply(list(NCoR1.0hr_vs_Emp.0hr,
                                         NCoR1.2hrCpG_vs_Emp.2hrCpG,
                                         NCoR1.6hrCpG_vs_Emp.6hrCpG,
                                         NCoR1.2hrpIC_vs_Emp.2hrpIC,
                                         NCoR1.6hrpIC_vs_Emp.6hrpIC,
                                         NCoR1.2hrCpG.pIC_vs_Emp.2hrCpG.pIC,
                                         NCoR1.6hrCpG.pIC_vs_Emp.6hrCpG.pIC),  
                              function(x) { x %>% 
                                            filter(abs(log2FoldChange) >=1 & padj  <=0.05) %>% 
                                            rownames_to_column("PeakID") %>% pull(PeakID) }) %>% 
                                            unique() %>% unlist() %>% unique()