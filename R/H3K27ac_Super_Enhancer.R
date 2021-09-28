# H3K27ac Super Enhancer analysis

SE = list.files(path="results/H3K27ac_SE/",
                recursive=TRUE,full.names=TRUE,pattern = "*summits_2kb_Gateway_SuperEnhancers.bed")
names(SE) = gsub("_summits_2kb_Gateway_SuperEnhancers.bed","",basename(SE))
peakAnnoList <- lapply(SE, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000),
                       annoDb= "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))
H3K27ac_SE.annotation.df = lapply(peakAnnoList, function(i) as.data.frame(i)) %>%
                              do.call(rbind,.) %>%
                              mutate(Clusters= gsub("\\..*","",rownames(.),""))%>%
                              remove_rownames(.)

H3K27ac_SE.annotation.df = H3K27ac_SE.annotation.df[!duplicated(H3K27ac_SE.annotation.df),]

# Super enhancer argument
SE.plot.list = list()

enhancerFile = "results/H3K27ac_SE/Emp_2hr_CpG/Emp_2hr_CpG_summits_2kb_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
enhancerName = "Emp_2hr_CpG"
Annotation = H3K27ac_SE.annotation.df %>% filter(Clusters == "Emp_2hr_CpG")
gene_label = labRow
SE.plot.list[[1]] = plot_SE(enhancerFile,enhancerName,Annotation,gene_label)

enhancerFile = "results/H3K27ac_SE/Emp_6hr_CpG/Emp_6hr_CpG_summits_2kb_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
enhancerName = "Emp_6hr_CpG"
Annotation = H3K27ac_SE.annotation.df %>% filter(Clusters == "Emp_6hr_CpG")
SE.plot.list[[2]] = plot_SE(enhancerFile,enhancerName,Annotation,gene_label)

enhancerFile = "results/H3K27ac_SE/Emp_2hr_pIC/Emp_2hr_pIC_summits_2kb_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
enhancerName = "Emp_2hr_pIC"
Annotation = H3K27ac_SE.annotation.df %>% filter(Clusters == "Emp_2hr_pIC")
SE.plot.list[[3]] = plot_SE(enhancerFile,enhancerName,Annotation,gene_label)

enhancerFile = "results/H3K27ac_SE/Emp_6hr_pIC/Emp_6hr_pIC_summits_2kb_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
enhancerName = "Emp_6hr_pIC"
Annotation = H3K27ac_SE.annotation.df %>% filter(Clusters == "Emp_6hr_pIC")
SE.plot.list[[4]] = plot_SE(enhancerFile,enhancerName,Annotation,gene_label)

enhancerFile = "results/H3K27ac_SE/Emp_2h_CpG_pIC/Emp_2h_CpG_pIC_summits_2kb_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
enhancerName = "Emp_2h_CpG_pIC"
Annotation = H3K27ac_SE.annotation.df %>% filter(Clusters == "Emp_2h_CpG_pIC")
SE.plot.list[[5]] = plot_SE(enhancerFile,enhancerName,Annotation,gene_label)


enhancerFile = "results/H3K27ac_SE/Emp_6h_CpG_pIC/Emp_6h_CpG_pIC_summits_2kb_12KB_STITCHED_ENHANCER_REGION_MAP.txt"
enhancerName = "Emp_6h_CpG_pIC"
Annotation = H3K27ac_SE.annotation.df %>% filter(Clusters == "Emp_6h_CpG_pIC")
SE.plot.list[[6]] = plot_SE(enhancerFile,enhancerName,Annotation,gene_label)

pdf("plots/H3K27ac_SE_plot.pdf",height = 8,width = 8)
do.call(grid.arrange, c(SE.plot.list,ncol=2,nrow=2))
dev.off()
