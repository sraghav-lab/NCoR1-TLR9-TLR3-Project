# Read H3K27ac peaks summit files
h3k27ac_peak_summit_files = list.files(path="../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/",pattern = "*_summits.bed")
h3k27ac_peak_summit.bed = lapply(h3k27ac_peak_summit_files, import.bed)
names(h3k27ac_peak_summit.bed) = h3k27ac_peak_summit_files

# Extend 1kb up and downstream to peak Summit
resizeRanges <- lapply(h3k27ac_peak_summit.bed, resize ,width = 1000,fix = 'center')

# Write 1kb extended peak file
write_bed= function(x){
  bed_file = as.data.frame(x)
  out_filename = unique(gsub("_peak_[0-9]*","",bed_file$name))
  write.table(bed_file[,c(1:3,6,7)],
              paste0("../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/",out_filename,"_summits_2kb.bed"),
              sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
lapply(resizeRanges, write_bed)

# merge 2kb  peak summit peaks
# bedops -m *_summits_2kb.bed >Emp_NCoR1_all_condition_merged_summits_2kb.bed
# add peak identifier 
# cat peak_file |parallel --verbose 'awk "{print(\$1,\$2,\$3,"Peak"NR)}" OFS="\t" {}_summits_2kb_nfr.bed >{}_summits_2kb_nfr.bed.bed'
# cat peak_file | parallel --verbose 'mv {}_summits_2kb_nfr.bed.bed {}_summits_2kb_nfr.bed'

# get differential H3K27ac merged peaks from all the control and NCoR1 KD condition
# working directory : ../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks 

###################################################### script to get differential peaks ################################################################################################
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_CpG/ ../tag_dir/Emp_0hr/ -F 2  >../diff_H3K27ac_peaks/Emp_6h_CpG_vs_Emp_0h_2fold_up_summits_2kb.txt
# Using fixed size peaks 
# Differential Peaks: 10009 of 67560 (14.81% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_CpG/ ../tag_dir/Emp_0hr/ -F 2 -rev >../diff_H3K27ac_peaks/Emp_6h_CpG_vs_Emp_0h_2fold_down_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 9736 of 67560 (14.41% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_CpG/ ../tag_dir/Emp_0hr/ -F 2 -rev >../diff_H3K27ac_peaks/Emp_2h_CpG_vs_Emp_0h_2fold_down_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 8746 of 67560 (12.95% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_CpG/ ../tag_dir/Emp_0hr/ -F 2  >../diff_H3K27ac_peaks/Emp_2h_CpG_vs_Emp_0h_2fold_up_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 8213 of 67560 (12.16% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_pIC/ ../tag_dir/Emp_0hr/ -F 2  >../diff_H3K27ac_peaks/Emp_6h_pIC_vs_Emp_0h_2fold_up_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 9693 of 67560 (14.35% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_pIC/ ../tag_dir/Emp_0hr/ -F 2 -rev  >../diff_H3K27ac_peaks/Emp_6h_pIC_vs_Emp_0h_2fold_down_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 5045 of 67560 (7.47% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_pIC/ ../tag_dir/Emp_0hr/ -F 2  >../diff_H3K27ac_peaks/Emp_2h_pIC_vs_Emp_0h_2fold_up_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 9754 of 67560 (14.44% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_pIC/ ../tag_dir/Emp_0hr/ -F 2 -rev >../diff_H3K27ac_peaks/Emp_2h_pIC_vs_Emp_0h_2fold_down_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 11918 of 67560 (17.64% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2h_CpG_pIC/ ../tag_dir/Emp_0hr/ -F 2  >../diff_H3K27ac_peaks/Emp_2h_CpG_pIC_vs_Emp_0h_2fold_up_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 10912 of 67560 (16.15% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2h_CpG_pIC/ ../tag_dir/Emp_0hr/ -F 2 -rev >../diff_H3K27ac_peaks/Emp_2h_CpG_pIC_vs_Emp_0h_2fold_down_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 12803 of 67560 (18.95% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6h_CpG_pIC/ ../tag_dir/Emp_0hr/ -F 2  >../diff_H3K27ac_peaks/Emp_6h_CpG_pIC_vs_Emp_0h_2fold_up_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 13181 of 67560 (19.51% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6h_CpG_pIC/ ../tag_dir/Emp_0hr/ -F 2 -rev >../diff_H3K27ac_peaks/Emp_6h_CpG_pIC_vs_Emp_0h_2fold_down_summits_2kb.txt
# Using fixed size peaks
# Differential Peaks: 12373 of 67560 (18.31% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_CpG/ ../tag_dir/Emp_2hr_pIC/ -F 2  >../diff_H3K27ac_peaks/Emp_2h_CpG_vs_Emp_2h_pIC_2fold_up_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 6906 of 67560 (10.22% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_CpG/ ../tag_dir/Emp_2hr_pIC/ -F 2 -rev >../diff_H3K27ac_peaks/Emp_2h_CpG_vs_Emp_2h_pIC_2fold_down_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 6772 of 67560 (10.02% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_CpG/ ../tag_dir/Emp_6hr_pIC/ -F 2  >../diff_H3K27ac_peaks/Emp_6h_CpG_vs_Emp_6h_pIC_2fold_up_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 5130 of 67560 (7.59% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_CpG/ ../tag_dir/Emp_6hr_pIC/ -F 2  -rev >../diff_H3K27ac_peaks/Emp_6h_CpG_vs_Emp_6h_pIC_2fold_down_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 8970 of 67560 (13.28% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_2hr_CpG/ ../tag_dir/Emp_2hr_CpG/ -F 2  -rev >../diff_H3K27ac_peaks/Emp_2h_CpG_vs_Emp_6h_CpG_2fold_up_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 0 of 67560 (0.00% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_CpG/ ../tag_dir/Emp_2hr_CpG/ -F 2   >../diff_H3K27ac_peaks/Emp_6h_CpG_vs_Emp_2h_CpG_2fold_up_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 3934 of 67560 (5.82% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_CpG/ ../tag_dir/Emp_2hr_CpG/ -F 2 -rev  >../diff_H3K27ac_peaks/Emp_6h_CpG_vs_Emp_2h_CpG_2fold_down_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 2799 of 67560 (4.14% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_pIC/ ../tag_dir/Emp_2hr_pIC/ -F 2   >../diff_H3K27ac_peaks/Emp_6h_pIC_vs_Emp_2h_pIC_2fold_up_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 7393 of 67560 (10.94% passed)
# getDifferentialPeaks Emp_NCoR1_all_condition_merged_summits_2kb.bed ../tag_dir/Emp_6hr_pIC/ ../tag_dir/Emp_2hr_pIC/ -F 2 -rev  >../diff_H3K27ac_peaks/Emp_6h_pIC_vs_Emp_2h_pIC_2fold_down_summits_2kb.txt 
# Using fixed size peaks
# Differential Peaks: 3067 of 67560 (4.54% passed)
############################################################################################################################################################################################################
# Run DESeq2 to generate normalized H3K27ac enrichment value across the replicates of all the samples
consensusToCount = import.bed("../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/Emp_NCoR1_all_condition_merged_summits_2kb.bed")
H3K27ac = read.csv("../..//NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/Emp_NCoR1_all_condition_merged_summits_2kb.annotation",
                   sep = "\t",header = T)
H3K27ac = H3K27ac[,c(16,1,2,3,4,5)]
colnames(H3K27ac)[c(1,2)] = c("Gene","GeneID")
bamsToCount <- dir("../../NCOR1/NCoR1_H3K27ac_analysis/bowtie_out/", full.names = TRUE, pattern = "*.25M.bam")
#bamsToCount <- subset(bamsToCount,!grepl("*_2h*", bamsToCount))
fcResults <- featureCounts(bamsToCount, 
                           annot.ext = H3K27ac[,c(2:6)],
                           isPairedEnd = FALSE,
                           countMultiMappingReads = FALSE,
                           maxFragLength = 100,nthreads = 15)

# Extract raw count                            
myCounts <- fcResults$counts
colnames(myCounts) = gsub(".25M.bam","",colnames(myCounts))

# Metadata
metaData <- data.frame(Group=colnames(myCounts), row.names = colnames(myCounts))

# Run DESEq2
Emp.H3K27ac_DDS_Rep1 <- DESeqDataSetFromMatrix(myCounts, metaData, ~Group, rowRanges = consensusToCount)
Emp.H3K27ac_DDS_Rep1 <- DESeq(Emp.H3K27ac_DDS)

# Extract variance stabilized transformed values
H3K27ac_vst_Rep1 <- vst(Emp.H3K27ac_DDS)
# Plot PCA
plotPCA(H3K27ac_vst_Rep1, intgroup = "Group", ntop = nrow(Emp.H3K27ac_DDS))

H3K27ac_vst_Rep1_df = as.data.frame(assay(H3K27ac_vst_Rep1))

H3K27ac_vst_Rep1_df = merge(H3K27ac,H3K27ac_vst_Rep1_df,by.x="GeneID",by.y=0,all.x="TRUE")
head(H3K27ac_vst_Rep1_df)
saveRDS(H3K27ac_vst_Rep1_df,file="../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/H3K27ac_Rlog.rds")

H3K27ac_Rlog_df =readRDS("../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/Emp_NCoR1_all_condition_merged_summits_2kb_Rlog.rds")

# differential H3K27ac peaks identified using HOMER function getDifferentialPeaks
down_h3k27ac_peak_list = list.files(path = "../diff_H3K27ac_peaks/",
                                    pattern = "_2fold_down_summits_2kb.txt",full.names = TRUE)
down_h3k27ac_peak_name = gsub("../diff_H3K27ac_peaks//","",down_h3k27ac_peak_list)
down_h3k27ac_peak_name = gsub("_summits_2kb.txt","",down_h3k27ac_peak_name)
down_h3k27ac_peak = lapply(down_h3k27ac_peak_list, read.csv,skip=18,
                           header=TRUE, stringsAsFactors=FALSE,sep="\t")
down_h3k27ac_peak = lapply(down_h3k27ac_peak, "[", c(1:4,10))
down_h3k27ac_peak = lapply(down_h3k27ac_peak, transform, log2FC = log2(Fold.Change.vs..Background))
names(down_h3k27ac_peak) = down_h3k27ac_peak_name

up_h3k27ac_peak_list = list.files(path = "../diff_H3K27ac_peaks/",
                                  pattern = "_2fold_up_summits_2kb.txt",full.names = TRUE)
up_h3k27ac_peak_name = gsub("../diff_H3K27ac_peaks//","",up_h3k27ac_peak_list)
up_h3k27ac_peak_name = gsub("_summits_2kb.txt","",up_h3k27ac_peak_name)
up_h3k27ac_peak = lapply(up_h3k27ac_peak_list, read.csv,skip=18,
                         header=TRUE, stringsAsFactors=FALSE,sep="\t")
up_h3k27ac_peak = lapply(up_h3k27ac_peak, "[", c(1:4,10))
up_h3k27ac_peak = lapply(up_h3k27ac_peak, transform, log2FC = log2(Fold.Change.vs..Background))
names(up_h3k27ac_peak) = up_h3k27ac_peak_name

Total_up_diff_peaks = unique(do.call(rbind,up_h3k27ac_peak)$X.PeakID)
Total_down_diff_peaks = unique(do.call(rbind,down_h3k27ac_peak)$X.PeakID)

Total_diff_peaks = as.data.frame(unique(c(Total_up_diff_peaks,Total_down_diff_peaks))) 
colnames(Total_diff_peaks)[1] = "PeakID"

Total_diff_peaks_list = c(up_h3k27ac_peak,down_h3k27ac_peak)
Total_diff_peaks = Total_diff_peaks_list %>% 
                            lapply(., "[", c(1,6)) %>% 
                            do.call(rbind,.) %>% rownames_to_column() 
Total_diff_peaks$rowname = gsub("\\.[0-9]*","",Total_diff_peaks$rowname)
Total_diff_peaks = reshape2::dcast(Total_diff_peaks,X.PeakID~rowname)
Total_diff_peaks[is.na(Total_diff_peaks)] <- 0
merge_H3K27ac_FC= function(col1){
  N <- ncol(col1)-1
  name <- colnames(col1[,-1])
  obj <- vector("list",ncol(col1)/2)
  k=1
  
  for(i in 1:N) {
    
    if(i%%2 ==1 && i <= N)    {
      
      #print(i)
      ID <-rowSums(col1[,c(i+1,i+2)]) 
      obj[[k]] <- ID
      nam <- colnames(col1)[i+1]
      nam <- str_replace(nam,"_2fold_down","")
      names(obj)[k] <- nam
      names(obj[[k]]) <- col1[,1]
      #print(k)
      k=k+1
    }
  }
  mat_merged <- as.data.frame(t(do.call(rbind, obj)))
  colnames(mat_merged) = names(obj)
  return(mat_merged)
}

Total_diff_peaks = merge_H3K27ac_FC(Total_diff_peaks) 


write.table(Total_diff_peaks,"../diff_H3K27ac_peaks/Total_diff_peaks.txt",
            sep = "\t",quote = FALSE,row.names = FALSE)

# Merge up and down group in excel readxl::read_xlsx("../diff_H3K27ac_peaks/Total_diff_peaks.xlsx")
# read excel file 
# readxl::read_xlsx("../diff_H3K27ac_peaks/Total_diff_peaks.xlsx")
Total_diff_peaks_merged = Total_diff_peaks 
colnames(Total_diff_peaks_merged)
h3K27ac_2k_summit_peak = read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/Emp_NCoR1_all_condition_merged_summits_2kb.bed",
                                  sep="\t",header = F)

Total_diff_peaks_merged_filtered = Total_diff_peaks_merged %>% 
                                        rownames_to_column("PeakID") %>%
                                        filter(PeakID %in% rownames(H3K27ac_vst_Rep1_df[rowSums(H3K27ac_vst_Rep1_df) >=100,]))
dim(Total_diff_peaks_merged_filtered)

diff_h3K27ac_2k_summit_peak = h3K27ac_2k_summit_peak  %>% filter(V4 %in% Total_diff_peaks_merged_filtered$PeakID)
head(diff_h3K27ac_2k_summit_peak)
write.table(diff_h3K27ac_2k_summit_peak,
              file = paste0("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered.bed"),
              quote = FALSE,sep = "\t",row.names = FALSE,
              col.names = FALSE)

# annotate Total diff_h3K27ac_2k_summit_peak
diff_h3K27ac_2k_summit_peak.bed = diff_h3K27ac_2k_summit_peak
colnames(diff_h3K27ac_2k_summit_peak.bed) =c("Chr","Start","End","PeakID") 
diff_h3K27ac_2k_summit_peak.bed= GRanges(diff_h3K27ac_2k_summit_peak.bed) 
diff_h3K27ac_2k_summit_peak.bed.annotation =  annotatePeak(diff_h3K27ac_2k_summit_peak.bed, 
                                                           TxDb=txdb,tssRegion=c(-1000, 1000),
                                                           annoDb="org.Mm.eg.db",
                                                           genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic")) 

plotAnnoBar(diff_h3K27ac_2k_summit_peak.bed.annotation) +gg_theme 


diff_h3K27ac_2k_summit_peak.bed.annotation.df = diff_h3K27ac_2k_summit_peak.bed.annotation %>% as.data.frame() 
############################################################################################################################################################################################################
Total_diff_peaks_merged_filtered = read.csv("../../NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered.bed",
                                            sep = "\t",header = F)
colnames(Total_diff_peaks_merged_filtered) = c("chr","Start","End","PeakID","Score")
consensusToCount = import.bed("../../NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered.bed")
H3K27ac = read.csv("../../NCOR1/NCoR1_H3K27ac_analysis/macs_narrow_peaks/Emp_NCoR1_all_condition_merged_summits_2kb.annotation",
                   sep = "\t",header = T)[,c(1,2,3,4,5)] 
colnames(H3K27ac)[1] = "GeneID"
H3K27ac = H3K27ac %>% filter(GeneID %in% Total_diff_peaks_merged_filtered$PeakID)
bams_Rep1 <- dir("../../NCOR1/NCoR1_H3K27ac_analysis/bowtie_out/", full.names = TRUE, pattern = "*.25M.bam")
bams_Rep2 <- dir("../../NCOR1/NCoR1_H3K27ac_Rep2_analysis/bowtie_out/", full.names = TRUE, pattern = "*.dup.filtered.srt.q10.bam$")
#bamsToCount <- subset(bamsToCount,!grepl("*_2h*", bamsToCount))
fcResults <- featureCounts(c(bams_Rep1,bams_Rep2), 
                           annot.ext = H3K27ac, isPairedEnd = FALSE,
                           countMultiMappingReads = FALSE, 
                           maxFragLength = 100,
                           nthreads = 15)
myCounts <- fcResults$counts 
#myCounts = myCounts[, -grep("2h", colnames(myC))]
colnames(myCounts)[2] ="Emp.2hr.CpG.pIC.25M.bam"
colnames(myCounts)[5] = "Emp.6hr.CpG.pIC.25M.bam"
colnames(myCounts)[9] ="NCoR1.2hr.CpG.pIC.25M.bam"
colnames(myCounts)[12] = "NCoR1.6hr.CpG.pIC.25M.bam"

colnames(myCounts) = gsub(".25M.bam|.Rep2.dup.filtered.srt.q10.bam","",colnames(myCounts))
colnames(myCounts)[c(15,24)] = c("NCoR1.2hr.CpG","Emp.0hr")
myCounts = myCounts[, -grep("NCoR1", colnames(myCounts))]
#myCounts = myCounts[,c(1,5,6,7,8,12,13,14,19,20,21,22,24,26,27,28)]

metaData <- data.frame(Group=colnames(myCounts))
metaData$Rep = factor(c(rep(1,7),rep(2,7)))
metaData$Stimulation = factor(gsub("Emp.|NCoR1.","",metaData$Group))
metaData$Genotype = factor(gsub("\\.[0-9].*","",metaData$Group))
metaData$Sample = paste0(metaData$Group,"_",metaData$Rep)
rownames(metaData)= metaData$Sample

colnames(myCounts) = metaData$Sample
metaData = metaData %>% arrange(Group)
rownames(metaData) = metaData$Sample

myCounts = myCounts[,metaData$Sample]

Emp.H3K27ac_DDS <- DESeqDataSetFromMatrix(myCounts, colData = metaData, design=~Group, rowRanges = consensusToCount)

Emp.H3K27ac_DDS <- DESeq(Emp.H3K27ac_DDS)
Emp.H3K27ac_vst <- vst(Emp.H3K27ac_DDS)
Emp.H3K27ac_vst

############################################################################################################################################################################################################
Emp.H3K27ac_vst_df = as.data.frame(assay(Emp.H3K27ac_vst)) %>% merge_rep(.) 
colnames(Emp.H3K27ac_vst_df)  = gsub("_1","",colnames(Emp.H3K27ac_vst_df))
Emp.H3K27ac_vst.df.scaled = Emp.H3K27ac_vst_df %>% apply(., 1, scale) %>% t() %>% as.data.frame()
colnames(Emp.H3K27ac_vst.df.scaled) = colnames(Emp.H3K27ac_vst_df)
############################################################################################################################################################################################################
# Calculated euclidena distance between H3K27ac ChIP-Seq samples
distsRL <- dist(t(assay(Emp.H3K27ac_vst)))
H3K27ac_mat<- as.matrix(distsRL)
rownames(H3K27ac_mat) <- colnames(H3K27ac_mat) <- with(colData(Emp.H3K27ac_DDS),
                                                       paste(Group,Rep, sep=' : '))
hc <- hclust(distsRL,method = "complete")
fviz_dend(hc, cex = 0.5)
my_palette <- colorRampPalette(c("red", "white", "blue"))(n = 50)
pdf("plots/H3K27ac_Emp_sample_hclust.pdf",height = 8,width = 10)
Heatmap(H3K27ac_mat,col = my_palette,name = "Euclidean Distance")
dev.off()
############################################################################################################################################################################################################
# Emp 2hr CpG vs Emp 0hr
Emp2hrCpG_vs_Emp0hr = results(Emp.H3K27ac_DDS,contrast=c("Group","Emp.2hr.CpG","Emp.0hr")) %>% as.data.frame()
#Emp2hrCpG_vs_Emp0hr.plt = 
H3K27ac.Volcano.plot(Emp2hrCpG_vs_Emp0hr,"2hr.CpG vs 0hr")
ggsave("plots/Diff_H3K27ac_scatter_plot_Emp2hrCpG.png",dpi=800,height = 5,width = 5.7)

# Emp 6hr CpG vs Emp 0hr
Emp6hrCpG_vs_Emp0hr = results(Emp.H3K27ac_DDS,contrast=c("Group","Emp.6hr.CpG","Emp.0hr")) %>% as.data.frame()
#Emp6hrCpG_vs_Emp0hr.plt =
H3K27ac.Volcano.plot(Emp6hrCpG_vs_Emp0hr,"6hr.CpG vs 0hr")
ggsave("plots/Diff_H3K27ac_scatter_plot_Emp6hrCpG.png",dpi=800,height = 5,width = 5.7)

# Emp 2hr pIC vs Emp 0hr
Emp2hrpIC_vs_Emp0hr = results(Emp.H3K27ac_DDS,contrast=c("Group","Emp.2hr.pIC","Emp.0hr")) %>% as.data.frame()
#Emp2hrpIC_vs_Emp0hr.plt = 
H3K27ac.Volcano.plot(Emp2hrpIC_vs_Emp0hr,"2hr.pIC vs 0hr")
ggsave("plots/Diff_H3K27ac_scatter_plot_Emp2hrpIC.png",dpi=800,height = 5,width = 5.7)

# Emp 6hr pIC vs Emp 0hr 
Emp6hrpIC_vs_Emp0hr = results(Emp.H3K27ac_DDS,contrast=c("Group","Emp.6hr.pIC","Emp.0hr")) %>% as.data.frame()
#Emp6hrpIC_vs_Emp0hr.plt = 
H3K27ac.Volcano.plot(Emp6hrpIC_vs_Emp0hr,"6hr.pIC vs 0hr")
ggsave("plots/Diff_H3K27ac_scatter_plot_Emp6hrpIC.png",dpi=800,height = 5,width = 5.7)

# Emp 2hr CpG+pIC vs Emp 0hr
Emp2hrCpGpIC_vs_Emp0hr = results(Emp.H3K27ac_DDS,contrast=c("Group","Emp.2hr.CpG.pIC","Emp.0hr")) %>% as.data.frame()
#Emp2hrCpGpIC_vs_Emp0hr.plt =
H3K27ac.Volcano.plot(Emp2hrCpGpIC_vs_Emp0hr,"2hr.CpG.pIC vs 0hr")

# Emp 6hr CpG+pIC vs Emp 0hr 
Emp6hrCpGpIC_vs_Emp0hr = results(Emp.H3K27ac_DDS,contrast=c("Group","Emp.6hr.CpG.pIC","Emp.0hr")) %>% as.data.frame()
#Emp6hrCpGpIC_vs_Emp0hr.plt =
H3K27ac.Volcano.plot(Emp6hrCpGpIC_vs_Emp0hr,"6hr.CpG.pIC vs 0hr")
dev.off()
pdf("plots/Diff_H3K27ac_scatter_plot_Emp.pdf",height = 3,width = 4)
plots = grid.arrange(Emp2hrCpG_vs_Emp0hr.plt,Emp6hrCpG_vs_Emp0hr.plt,
             Emp2hrpIC_vs_Emp0hr.plt,Emp6hrpIC_vs_Emp0hr.plt,
             Emp2hrCpGpIC_vs_Emp0hr.plt,Emp6hrCpGpIC_vs_Emp0hr.plt,
             ncol=2,nrow=3)
############################################################################################################################################################################################################
############################################################################################################################################################################################################
# Average H3K27ac enrichment at gene level
H3K27ac_vst.df.scaled_gene_level = merge(diff_h3K27ac_2k_summit_peak.bed.annotation.df,Emp.H3K27ac_vst.df.scaled,
                                         by.x="PeakID",by.y=0) %>% 
                                        dplyr::select(c(17:ncol(.))) %>% 
                                        group_by(SYMBOL) %>% 
                                        summarise_all(list(mean)) %>% as.data.frame()

############################################################################################################################################################################################################
# # Likelihood ratio test
Emp.H3K27ac_DDS_lrt <- DESeq(Emp.H3K27ac_DDS, test="LRT",reduced = ~1)
# # Extract results
H3K27ac_res_LRT <- results(Emp.H3K27ac_DDS_lrt)
# # Subset the LRT results to return genes with padj < 0.05 and log2FoldChange >0.58
H3K27ac_sig_res_LRT <- H3K27ac_res_LRT %>%
                        data.frame() %>%
                        rownames_to_column(var="Peak") %>% 
                        as_tibble() %>% 
                        filter(padj < 0.05 & abs(log2FoldChange) >0.58 )
# 
# # Get sig gene lists
H3K27ac_sigLRT_genes <- H3K27ac_sig_res_LRT %>% pull(Peak)
# 
length(H3K27ac_sigLRT_genes)
# 
# 
H3K27ac_clustering_sig_genes <- H3K27ac_sig_res_LRT %>%
                                  arrange(padj) %>%
                                  filter(abs(H3K27ac_sig_res_LRT$log2FoldChange) > 0.58)

# # Obtain vst values for those significant Peaks
H3K27ac_cluster_vst <- as.data.frame(assay(Emp.H3K27ac_vst))[H3K27ac_clustering_sig_genes$Peak, ]

# Clustering 
H3K27ac_clusters <- degPatterns(H3K27ac_cluster_vst, metadata = metaData,time = "Group",col = NULL)

table = H3K27ac_clusters[["normalized"]]
time = metaData$Group
color = NULL
min_genes = 10
process = FALSE
points = TRUE
boxes = TRUE
smooth = TRUE
lines = TRUE
facet = TRUE
cluster_column = "cluster"
prefix_title = "Group: "


if (cluster_column  %in% colnames(table)){
  table[["cluster"]] = table[[cluster_column]]
}
if (process){
  table <- .process(table, time, color)
}

if ("cluster"  %in% colnames(table)){
  counts <- table(distinct(table, genes, cluster)[["cluster"]])
  counts <- counts[counts>=min_genes]
  if (length(counts)==0)
    stop("No clusters with min_genes > ", min_genes)
  table <- inner_join(table,
                      data.frame(cluster = as.integer(names(counts)),
                                 title = paste(prefix_title,
                                               names(counts),
                                               "- genes:" ,
                                               counts),
                                 stringsAsFactors = FALSE),
                      by = "cluster")
}
table$title = factor(table$title,levels = str_sort(unique(table$title),numeric = TRUE))
if (is.null(color)){
  color = "dummy"
  table[[color]] = ""
  lines = FALSE
}
table[["line_group"]] = paste(table[["genes"]],
                              table[[color]])


splan <- length(unique(table$Sample)) - 1L

table$Sample = factor(table$Sample,levels = c("Emp.0hr_1","Emp.2hr.CpG_1",
                                              "Emp.6hr.CpG_1","Emp.2hr.pIC_1",
                                              "Emp.6hr.pIC_1","Emp.2hr.CpG.pIC_1",
                                              "Emp.6hr.CpG.pIC_1"))

############################################################################################################################################################################################################
H3K27ac_cluster_lrt_df = H3K27ac_clusters$df

H3K27ac_cluster_lrt_df = H3K27ac_cluster_lrt_df %>%
                               mutate(Condition = case_when(cluster %in% c(2,5,10,25,28,29) ~ "pIC",
                                                            cluster %in% c(1,6,15) ~ "CpG",
                                                            cluster %in% c(3,9,12,14,19,20,26,32,34,35) ~ "CpG_pIC",
                                                            cluster %in% c(4,7,8,11,13,16,17,18,21,22,23,24,27,33) ~ "Down"))


table = table %>% 
            mutate(Condition = case_when(cluster %in% c(2,5,10,25,28,29) ~ "pIC",
                                         cluster %in% c(1,6,15) ~ "CpG",
                                         cluster %in% c(3,9,12,14,19,20,26,32,34,35) ~ "CpG_pIC",
                                         cluster %in% c(4,7,8,11,13,16,17,18,21,22,23,24,27,33) ~ "Down")
                                         )

# Write bed file of H3K27ac cluster
for ( i in unique(H3K27ac_cluster_lrt_df$Condition)) {
    diff_h3K27ac_2k_summit_peak %>% 
    filter(V4 %in% H3K27ac_cluster_lrt_df[which(H3K27ac_cluster_lrt_df$Condition ==i),]$genes) %>% 
    dplyr::select(c(1,2,3,4)) %>%
    write.table(.,file= paste0("../../NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/","H3K27ac_",i,".bed"),
                sep = "\t",quote = FALSE,col.names = FALSE,row.names = FALSE)
  
}

H3K27ac_cluster_vst_scaled= H3K27ac_vst.df.scaled[H3K27ac_cluster_lrt_df$genes,]

# Heatmap
# split cluster
split <- factor(H3K27ac_cluster_lrt_df$Condition)
H3K27ac_cluster_hmap <- Heatmap(H3K27ac_cluster_vst_scaled,km = 1,
                              name="Z-score",
                              col= colorRamp2(c(-1.5,0,1.5),c("blue","white","red")),
                              heatmap_legend_param=list(at=c(-1.5,0,1.5),color_bar="continuous", 
                                                        legend_direction="horizontal", legend_width=unit(5,"cm"),
                                                        title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),

                              #Split heatmap by Cluster
                              split = split,
                              border=TRUE,
                              # box colour
                              #rect_gp = gpar(col = "black"),

                              #Row annotation configurations
                              cluster_rows=TRUE,
                              cluster_row_slices = FALSE,
                              show_row_dend=TRUE,
                              #row_title="", 
                              row_title_side="left",
                              row_title_gp=gpar(fontsize=8),
                              show_row_names=FALSE,
                              row_names_side="right",
                              #row_title_rot=30,

                              #Column annotation configuratiions
                              cluster_columns=FALSE,
                              show_column_dend=TRUE,
                              #column_title="Samples",
                              column_title_side="top",
                              column_names_rot = 45, 
                              column_names_side = "top",
                              column_title_gp=gpar(fontsize=10, fontface="bold"),
                              #column_title_rot=45,
                              column_names_gp = gpar(fontsize = 10, fontface="bold"),
                              show_column_names=TRUE,

                              #Dendrogram configurations: columns
                              clustering_distance_columns="euclidean",
                              clustering_method_columns="single",
                              column_dend_height=unit(10,"mm"),

                              #Dendrogram configurations: rows
                              clustering_distance_rows="euclidean",
                              clustering_method_rows="single",
                              row_dend_width=unit(10,"mm"))

pdf("plots/H3K27ac_diif_peaks_Emp_heatmap.pdf",height = 7,width = 5)
draw(H3K27ac_cluster_hmap ,heatmap_legend_side="right")
dev.off()
############################################################################################################################################################################################################

# Generate boxplot of each H3K27ac Cluster 
my_comparisons = list(c("Emp.2hr.CpG","Emp.2hr.pIC"),
                      c("Emp.6hr.CpG","Emp.6hr.pIC"),
                      c("Emp.0hr","Emp.2hr.pIC"),
                      c("Emp.0hr","Emp.6hr.CpG"))

H3K27ac_cluster_length = lengths(split(H3K27ac_cluster_lrt_df$genes,H3K27ac_cluster_lrt_df$Condition))
H3K27ac_cluster_facet_label = c(paste0("CpG(",H3K27ac_cluster_length[1],")"),
                                paste0("CpG & pIC(",H3K27ac_cluster_length[2],")"),
                                paste0("Down(",H3K27ac_cluster_length[3],")"),
                                paste0("pIC(",H3K27ac_cluster_length[4],")"))
names(H3K27ac_cluster_facet_label) =c("CpG","CpG_pIC","Down","pIC")

pdf("plots/H3K27ac_diif_peaks_Emp_boxplot.pdf",height = 12,width = 5)
merge(H3K27ac_cluster_lrt_df,H3K27ac_vst_df,by.x="genes",by.y=0) %>% 
  dplyr::select(c(3:10)) %>%
  melt(.,id.vars="Condition") %>% 
  mutate(Condition = factor(Condition,levels = c("CpG","CpG_pIC","pIC","Down")),
         variable = factor(variable, levels = c("Emp.0hr", "Emp.2hr.CpG", "Emp.6hr.CpG",
                                                "Emp.2hr.pIC","Emp.6hr.pIC",
                                                "Emp.2hr.CpG.pIC", "Emp.6hr.CpG.pIC"))) %>%
  ggplot(.,aes(x=variable,y=value,color=variable))+
  geom_boxplot(lwd=0.2,outlier.shape=NA)+
  facet_wrap(~Condition,ncol = 1,labeller = labeller(Condition = H3K27ac_cluster_facet_label))+
  ylab("Noramlized vst (DESeq2)")+
  ylim(0,16)+
  stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",y_position = 14,method = "wilcox",paired = FALSE,size=4) +
  gg_theme_xrot
dev.off()
############################################################################################################################################################################################################
# Add gene annotation to H3K27ac Cluster
H3K27ac_cluster_annotation =  merge(H3K27ac_cluster_lrt_df,diff_h3K27ac_2k_summit_peak.bed.annotation.df,by.x="genes",by.y="PeakID") 
H3K27ac_cluster_annotation %>% 
  dplyr::select(c(4,5,6,1)) %>%
  write.table(.,"../../NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/H3K27ac_diff_peaks_LRT.bed",
              quote = FALSE,row.names = FALSE,col.names = FALSE,sep = "\t")

H3K27ac_cluster_vst_gene_level =  merge(H3K27ac_cluster_annotation,Emp.H3K27ac_vst.df.scaled,by.x="genes",by.y=0) %>%
                                        dplyr::select(c(19,21:ncol(.))) %>% #head()
                                        group_by(SYMBOL) %>%
                                        summarise_all(list(mean)) %>%
                                        as.data.frame()

diff_h3k27ac_cluster_genes.list =list()
for (i in 1:4){
    c = unique(H3K27ac_cluster_annotation$Condition)[i]
    cluster= H3K27ac_cluster_annotation %>% 
                filter(Condition ==c) %>% 
                dplyr::select("SYMBOL") %>% 
                unique(.) %>% 
                pull(.,"SYMBOL")
  diff_h3k27ac_cluster_genes.list[[i]]= cluster
  names(diff_h3k27ac_cluster_genes.list)[i] = unique(H3K27ac_cluster_annotation$Condition)[i]
}
diff_h3k27ac_cluster_genes.list = diff_h3k27ac_cluster_genes.list[lapply(diff_h3k27ac_cluster_genes.list,length)>0]
gene_df_toupper = lapply(diff_h3k27ac_cluster_genes.list, toupper) 
term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
compEnricher =compareCluster(gene_df_toupper, fun="enricher",TERM2GENE=term2gene)

# Pathway barplot
pdf("plots/H3K27ac_diif_peaks_Emp_pathway.pdf",height = 13,width = 9)
compEnricher %>%
  as.data.frame() %>% #head()
  group_by(Cluster) %>% #View()
  slice(1:6) %>%
  ungroup %>% 
  mutate(Description = gsub("REACTOME|_"," ",Description),
         Description = reorder_within(Description, -pvalue, Cluster)) %>% 
  
  as.data.frame() %>% 
  #filter(Cluster %in% lrt_cluster) %>%
  mutate(Cluster = factor(Cluster,levels = c("CpG","CpG_pIC","pIC","Down"))) %>%
  #.[c(-4,-5,-7,-8,-10),] %>%
  #Description = reorder(Description, pvalue))
  ggplot(aes(x=reorder(Description,-pvalue), y=-log(pvalue), fill = p.adjust)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_wrap(~Cluster, scales = "free",ncol = 1) +
  #scale_x_reordered() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  scale_fill_gradient(low="red",high="blue")+
  scale_y_continuous(expand = c(0,0)) +
  gg_theme
dev.off()

############################################################################################################################################################################################################
save(diff_h3K27ac_2k_summit_peak.bed.annotation.df, 
     H3K27ac_cluster_lrt_df,
     Emp.H3K27ac_vst.df.scaled,
     H3K27ac_vst.df.scaled_gene_level,
     H3K27ac_cluster_vst_gene_level,
     file = "results/Control_H3K27ac.RData")



