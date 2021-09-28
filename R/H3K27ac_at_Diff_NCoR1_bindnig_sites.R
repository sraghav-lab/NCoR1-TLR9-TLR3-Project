# ##########################################################################################
# Ologram output plot
overlap = list.files(path = "/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/H3k27ac_at_Diff_NCoR1_sites/ologram_output/",
                     pattern = "*.tsv",full.names = TRUE)
overlap.names = lapply(overlap, basename) %>% unlist() %>% gsub("Total_diff_peaks_rlog_filtered_|.tsv","",.)
overlap.files = lapply(overlap,read.csv, sep="\t",header=T)
overlap.files = lapply(overlap.files, "[", c(1,2,3,6,7))
names(overlap.files) = overlap.names
for (i in 1:length(overlap.files)) {
       overlap.files[[i]]$feature_type = gsub("[Query + | + ... ]","",overlap.files[[i]]$feature_type)
       overlap.files[[i]]$feature_type = gsub("\\[|\\]","",overlap.files[[i]]$feature_type)
       
}
overlap.df= bind_rows(overlap.files, .id = "column_label") %>% melt(,measure.vars = c("nb_intersections_expectation_shuffled","nb_intersections_true"))
overlap.df$nb_intersections_variance_shuffled[overlap.df$variable == "nb_intersections_true"] <- 0  
overlap.df$nb_intersections_pvalue[overlap.df$variable == "nb_intersections_expectation_shuffled"] <- ""
overlap.df$variable = ifelse(overlap.df$variable == "nb_intersections_true","True","Shuffled")
overlap.df$sd = sqrt(overlap.df$nb_intersections_variance_shuffled)
overlap.df$feature_type = factor(overlap.df$feature_type,
                                 levels = c("2fold_p_CpG_pIC_all3_vs_Uns","2fold_p_all3_vs_CpG_pIC",
                                            "2fold_p_pIC_vs_CpG_all3","2fold_no_chang_ns_CpG_pIC_all3",
                                            "2fold_p_Uns_vs_CpG_pIC_all3"),
                                labels = c("Cluster I","Cluster II","Cluster III","Cluster IV","Cluster V"))

pdf("plots/H3K27ac_NCoR1_overlap_binding_diff_peaks.pdf", height = 6, width = 14)
overlap.df %>% #filter(column_label == "H3K27ac_Down") %>%
ggplot(.,aes(x=feature_type,y=value,fill=variable))+
       geom_bar(stat="identity",position=position_dodge())+
       geom_errorbar(aes(ymin = value-sd, ymax = value+sd), position =position_dodge(.9), width = 0.25)+
       geom_text(aes(label=nb_intersections_pvalue),size=3,hjust=0.5,vjust=0.1, position=position_dodge(width=0.9), angle=0)+
       facet_wrap(~column_label,ncol = 8,scales = "free")+
       gg_theme+
       theme(axis.text.x = element_text(angle=40,vjust = 1,hjust = 1))+
       xlab("NCoR1 Binding")+
       ylab("No. of NCoR1 peaks Overlap with H3K27ac")
dev.off()

overlap.df.perc = overlap.df %>% filter(variable == "True") %>% dplyr::select(c(1,2,6))
overlap.df.perc = overlap.df.perc %>% dplyr::group_by(column_label) %>% dplyr::mutate(percent = value/sum(value) *100 )
overlap.df.perc$column_label = factor(overlap.df.perc$column_label,levels = c("H3K27ac_CpG","H3K27ac_CpG_pIC","H3K27ac_pIC","H3K27ac_Down"))
pdf("plots/H3K27ac_NCoR1_overlap_binding_diff_peaks_percent.pdf", height = 7, width = 7)
ggplot(overlap.df.perc,aes(x=column_label,y=percent,fill=feature_type))+
        geom_bar(stat="identity")+
        gg_theme_xrot+
                geom_text(aes(y=percent,label=paste0(round(percent,2),"%")),position=position_stack(vjust=0.5))+
        labs(y= "percent of H3K27ac bound by NCoR1",fill = "NCoR1")
dev.off()
#####################################################################################################################
files = list.files(path = "/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/H3k27ac_at_Diff_NCoR1_sites/",
                   pattern = "H3K27ac_at.*.bed$",full.names = TRUE)
#files =  c("CpG_pIC_all3.bed","CpG_all3.bed","Uns_pIC.bed","Uns_vs_CpG_pIC_all3.bed","no-change_uns_CpG_pIC_all3.bed")
names(files) = c("No_Change","CpG_all3","CpG_pIC_all3","Uns_pIC","Uns")
#files = list.files("NCoR1_Uns_Unique.bed","NCoR1_CpG_Unique.bed","NCoR1_Uns_CpG_common.bed","NCoR1_SMRT_Uns_CpG_6hr_common_peaks.bed","SMRT_Uns_CpG_common.bed","SMRT_Uns_Unique.bed")
peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000),
                       annoDb= "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))

peak_annotation = lapply(peakAnnoList, annotation)
peak_annotation = do.call(cbind,peak_annotation)[c(1,3,4,5),c(1,2,4,6,8,10)]
colnames(peak_annotation) = gsub("\\.Freq","",colnames(peak_annotation))
colnames(peak_annotation)[1] = "Features"
peak_annotation_melt = melt(peak_annotation)
#peak_annotation_melt$variable = factor(peak_annotation_melt$variable,levels = c("SMRT_Uns","SMRT_Uns_CpG","SMRT_CpG","NCoR1_SMRT_CpG","NCoR1_SMRT_Uns_CpG","NCoR1_CpG","NCoR1_Uns"))
ggplot(peak_annotation_melt,aes(x=variable,y=value,fill=Features))+
       geom_bar(stat = "identity")+
       coord_flip() +
       gg_theme

diff_NCoR1.annotation.df = lapply(peakAnnoList, function(i) as.data.frame(i)) %>%
       do.call(rbind,.) %>%
       mutate(Clusters= gsub("\\..*","",rownames(.),""))

diff_NCoR1.annotation.df_10kb = diff_NCoR1.annotation.df %>% filter(abs(distanceToTSS) <= 10000)
diff_NCoR1.annotation.df.list =list()
for (i in unique(diff_NCoR1.annotation.df$Clusters)){
       
       cluster = diff_NCoR1.annotation.df %>% 
                        filter(Clusters ==i ) %>% 
                        dplyr::select("SYMBOL") %>%
                        unique(.) %>% 
                        pull(.,"SYMBOL")
       diff_NCoR1.annotation.df.list[[i]]= cluster
}




######################################################################################################################
# Density plot of H3K27ac at Diffeerntial NCoR1 binding sites.
dat = list.files(pattern = "*genes.peaks.density.txt",
                 path = "/home/imgsb/Gyan/NCOR1/NCoR1_H3K27ac_analysis/H3k27ac_at_Diff_NCoR1_sites/",
                 full.names = TRUE)
head(dat)
dat = dat[c(3,4,2,1,5)]
#CpG.dat <- read.csv("CpG_NCoR1_pu1.tag_density", sep = "\t",header =T,na.strings = "NA")
#colnames
plotlist =list()

for(i in 1:5){
c=2      
dat1 <- read.csv(dat[i],sep = "\t",header =T,na.strings = "NA")
#df1 <- dat1[,c(1,2,5,8,11,14,17,20,23,26,29,32,35,38,41)]
df1 <- dat1[,c(1,2,5,8,11,14,17,20,23,26,29,32,35,38,41)]
colname <- gsub("...tag_dir.","",colnames(df1))
colname <- gsub("..Coverage","",colname)
colname[1] <- "Dist"
colnames(df1) <- colname
        ptitle = gsub("_CpG_pIC_genes.peaks.density.txt","",basename(dat[i]))
        df2 <- melt(df1,id="Dist",variable.name = "Stimulation")
        head(df2)
        df2$condition  = gsub("Emp_|NCoR1_","",df2$Stimulation)
        df2$condition  = factor(df2$condition,levels = c("0hr","2hr_CpG","6hr_CpG","2hr_pIC","6hr_pIC","2h_CpG_pIC","6h_CpG_pIC"))
        df2$genotype   = gsub("_.*","",df2$Stimulation)
        x = min(df2$value)
        y = max(df2$value)
        #pname <- paste0("Plot-12041_no_change",i)
        den1 <- df2 %>% #dplyr::filter((grepl("2h",condition))) %>% 
                ggplot(aes(x=Dist  ,y=value,color=genotype))+
                geom_line(size=1) +
                gg_theme+
                facet_wrap(~condition,ncol = 7)+
                theme_bw()+
                ylim(x,y)+
                scale_color_manual(values = c("blue","red"))+
                theme(axis.line.x = element_line(color="black", size = 0.5),
                        plot.title = element_text(face="italic"),
                        axis.line.y = element_line(color="black", size = 0.5),
                        axis.text.y=element_text(size=10,colour = "black"),
                        axis.text.x=element_text(size=10,colour = "black"),
                        axis.title.y = element_text(size=10,colour = "black"),
                        axis.title.x = element_text(size=10,colour = "black"),
                        strip.text.x = element_text(size = 15, colour = "black"),
                        legend.text=element_text(size=10))+
                        #legend.position = "none")+
                ggtitle(ptitle)+
                labs(x="Distance from NCoR1 peaks summit (bp)",y="Normalized Tag Count)")
        plotlist[[i]] = den1
       # if( c < 15){
               #c = c +2 
        #}
}
pdf("plots/H3K27ac_at_all_diff_NCoR1_sites.pdf",height = 15,width = 14)
do.call(grid.arrange,  c(plotlist,ncol=1))
dev.off()
################################################################################################################################

dat = read.csv("2fold_up_all3_vs_CpG_pIC.density.txt",header  = T,sep="\t",row.names = 1)[,c(19:36)]
colnames(dat)
#colnames(dat)[1] = "PeakID"
colnames(dat) = gsub(".*tag_dir\\.","",colnames(dat))
colnames(dat) = gsub("\\.*Tag.Count.*","",colnames(dat))
colnames(dat) = gsub("......Analysis_result_of_NCoR1_ChIPSeq.NCOR1_mm10_alignment.Uniquley_mapped_reads.subsample_alignment.tag_directory.","",colnames(dat))
head(dat)
res= cor(dat[,c(1,2,3,4,5,9,13,17)])
Heatmap(res)

