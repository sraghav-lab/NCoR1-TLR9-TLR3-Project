############################################################################################################################################################################################################
# Overlap of Diff peaks with RelA and IRF3
# Directory : /home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks
# parallel --verbose 'echo {} $(bedtools intersect -a {} -b /media/Hard_disk1/Analysis_backup/Analysis_result_of_ChIPseq_Lajos/IRF3_peaks/IRF3_pIC_90min_Rep_merged_200bp.bed -wa |bedops -m - | wc -l)' ::: Total_diff_peaks_rlog_filtered_cluster_[0-9]_nfr.bed  >Total_diff_peaks_rlog_filtered_cluster_IRF3_90min_overlap.txt
# parallel --verbose 'echo {} $(bedtools intersect -a {} -b /home/imgsb/Gyan/NCOR1/Analysis_result_of_NCoR1_TF/Rela_0hr_2hr_6hr_CpG_peaks/Rela_0hr_CpG_peaks.narrowPeak -wa |bedops -m - | wc -l)' ::: Total_diff_peaks_rlog_filtered_cluster_[0-9]_nfr.bed >Total_diff_peaks_rlog_filtered_cluster_Rela_2hr_overlap.txt

Total_cluster_nfr = read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered_cluster_1-9_nfr.txt",
                             sep="\t",header = F)[,c(2,1)]
colnames(Total_cluster_nfr)= c("V1","V2")
Total_cluster_nfr$V1 = gsub("Total_diff_peaks_rlog_filtered_|_nfr.bed","",Total_cluster_nfr$V1)
Rela = read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered_cluster_Rela_2hr_overlap.txt",
                sep=" ",header = F)
Rela$V1 = gsub("Total_diff_peaks_rlog_filtered_|_nfr.bed","",Rela$V1)
Irf3 = read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered_cluster_IRF3_90min_overlap.txt",
                sep=" ",header = F)
Irf3$V1 = gsub("Total_diff_peaks_rlog_filtered_|_nfr.bed","",Irf3$V1)

Rela_Irf3 = cbind(Total_cluster_nfr,Rela,Irf3)[,c(1,2,4,6)]
colnames(Rela_Irf3) = c("Cluster","Total_nfr","Rela","Irf3")
Rela_Irf3$Rela_overlap = Rela_Irf3$Rela/Rela_Irf3$Total_nfr*100
Rela_Irf3$Irf3_overlap = Rela_Irf3$Irf3/Rela_Irf3$Total_nfr*100
Rela_Irf3_melt = melt(Rela_Irf3[,c(1,5,6)])
Rela_Irf3_melt$Factor= gsub("_overlap","",Rela_Irf3_melt$variable)
Rela_Irf3_melt$Cluster = factor(Rela_Irf3_melt$Cluster,levels = c("cluster_3","cluster_5","cluster_4","cluster_6","cluster_7","cluster_8",
                                                                  "cluster_1","cluster_2"))
ggplot(Rela_Irf3_melt,aes(x=Factor,y=value,fill=Factor)) + 
  facet_wrap(~Cluster,ncol = 1)+ 
  geom_bar(stat = "identity") +
  gg_theme +
  ggsave("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered_cluster_Rela_Irf3_overlap.png",height = 25,width = 5)

############################################################################################################################################################################################################

# IRF3 and RelA density and motif score
# cat Total_diff_peaks_rlog_filtered_cluster_[0-9].bed >Total_diff_peaks_rlog_filtered_cluster.bed

# bedtools intersect -a Total_diff_peaks_rlog_filtered_cluster.bed -b ../../Analysis_result_of_NCoR1_TF/Rela_0hr_2hr_6hr_CpG_peaks/Rela_0hr_CpG_peaks.narrowPeak -wa >Total_diff_peaks_rlog_filtered_cluster_Rela.bed
# annotatePeaks.pl Total_diff_peaks_rlog_filtered_cluster_Rela.bed mm10 -d ../../Analysis_result_of_NCoR1_TF/tagdir/Rela_2hr/ /media/Hard_disk1/Analysis_backup/Analysis_result_of_ChIPseq_Lajos/IRF3_pIC_90min_Rep1/ -m ../motif_analysis/Total_diff_peaks_rlog_filtered_cluster_3_nfr/homerResults/motif1.motif  -log >Total_diff_peaks_rlog_filtered_cluster_RelA.density.txt

# bedtools intersect -a Total_diff_peaks_rlog_filtered_cluster.bed -b /media/Hard_disk1/Analysis_backup/Analysis_result_of_ChIPseq_Lajos/IRF3_pIC_90min_Rep1/ -wa >Total_diff_peaks_rlog_filtered_cluster_IRF3.bed
# annotatePeaks.pl Total_diff_peaks_rlog_filtered_cluster_IRF3.bed mm10 -d ../../Analysis_result_of_NCoR1_TF/tagdir/Rela_2hr/ /media/Hard_disk1/Analysis_backup/Analysis_result_of_ChIPseq_Lajos/IRF3_pIC_90min_Rep1/ -m ../motif_analysis/Total_diff_peaks_rlog_filtered_cluster_3_nfr/homerResults/motif1.motif  -log >Total_diff_peaks_rlog_filtered_cluster_IRF3.density.txt


diff_h3k27ac = read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered_cluster_RelA.density.txt",
                        sep="\t",header = T)[,c(1,22:24)]
head(diff_h3k27ac)
colnames(diff_h3k27ac) = c("peakid","RelA","Irf3","NFkB-p65")
rownames(diff_h3k27ac) =diff_h3k27ac$peakid

diff_h3k27ac =merge(diff_h3k27ac,x,by=0) 

diff_h3k27ac_filtered = diff_h3k27ac %>%
  filter( str_detect(`NFkB-p65`,"TC")) 

diff_h3k27ac_filtered = diff_h3k27ac_filtered %>% 
  mutate(Group = case_when(str_detect(`NFkB-p65`,"TC") ~ "RelA"))

ggplot() +
  # draw the original data series with grey
  geom_point(aes(Emp.2hr.CpG, Emp.2hr.pIC), data = x, colour = alpha("grey", 0.5),size =0.5) +
  # colourise only the filtered data
  geom_point(aes(Emp.2hr.CpG, Emp.2hr.pIC, colour = Group), data = diff_h3k27ac_filtered,size=0.5) +
  #scale_color_gradient(low="blue",high = "red")+
  scale_color_manual(values = "blue")+
  gg_theme+
  #scale_color_viridis(option = "D")+
  guides(colour = guide_legend(override.aes = list(size=10)))+
  theme(text = element_text(size = 25))



# Tag density at  Overlap of Diff peaks with RelA and IRF3
# parallel --verbose 'bedtools intersect -b {} -a /media/Hard_disk1/Analysis_backup/Analysis_result_of_ChIPseq_Lajos/IRF3_peaks/IRF3_pIC_90min_Rep_merged_200bp.bed -wo |cut -f 1,2,3,11' ::: Total_diff_peaks_rlog_filtered_cluster_[0-9].bed >Total_diff_peaks_rlog_filtered_cluster_IRF3_90min.bed
# annotatePeaks.pl Total_diff_peaks_rlog_filtered_cluster_IRF3_90min.bed mm10 -d /media/Hard_disk1/Analysis_backup/Analysis_result_of_ChIPseq_Lajos/IRF3_pIC_90min_Rep1/ >Total_diff_peaks_rlog_filtered_cluster_IRF3_90min_density.txt

# parallel --verbose 'bedtools intersect -b {} -a ../../Analysis_result_of_NCoR1_TF/Rela_0hr_2hr_6hr_CpG_peaks/Rela_0hr_CpG_peaks.narrowPeak -wo |cut -f 1,2,3,14' ::: Total_diff_peaks_rlog_filtered_cluster_[0-9].bed >Total_diff_peaks_rlog_filtered_cluster_RelA_2hr_CpG.txt
# 


diff_h3k27ac_Irf3 = read.csv("/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/Total_diff_peaks_rlog_filtered_cluster_IRF3_90min_density.txt",
                             sep="\t",header = T)[,c(1,20)]
colnames(diff_h3k27ac_Irf3) = c("peakid","Irf3")
diff_h3k27ac_Irf3$peakid =gsub("-[0-9]*","",diff_h3k27ac_Irf3$peakid)
diff_h3k27ac_Irf3 = diff_h3k27ac_Irf3 %>%
  group_by(peakid) %>%
  summarize(Irf3=mean(Irf3)) %>% as.data.frame(.) 
rownames(diff_h3k27ac_Irf3) = diff_h3k27ac_Irf3$peakid
my_comparisons <- list( c("CpG", "pIC"), c("CpG and pIC", "CpG"),c("CpG and pIC","pIC"))
p2 = merge(diff_h3k27ac_Irf3,x,by=0) %>% 
                    filter(subgroup == "Cluster4" |subgroup == "Cluster6" |subgroup == "Cluster3" |subgroup == "Cluster5" ) %>%
  
                    mutate(Group = case_when(subgroup == "Cluster3" ~ "CpG",
                           subgroup == "Cluster4" | subgroup == "Cluster6" ~ "pIC",
                           subgroup == "Cluster5" ~ "CpG and pIC")) %>%
ggplot(aes(x=Group,y=Irf3,fill=Group))+
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",
                     y_position = 9,method = "wilcox",paired = FALSE,)+
  
  scale_fill_manual(values = c("blue", "green","red"))+
  gg_theme+
  theme(text = element_text(size = 25),
        axis.text.x = element_blank())+
  labs(x="",y="Log ( Normalized Tag Count)",title = "IRF3 90min CpG")
  
grid.arrange(p1,p2,ncol=2)

############################################################################################################################################################################################################
diff_h3k27ac_cluster = list.files("../NCoR1_H3K27ac_analysis/diff_H3K27ac_peaks/",pattern = "Total_diff_peaks_rlog_filtered_cluster_[0-9].bed")
