# Differential expression analysis of NCoR1 KD 2hr pIC vs Control 2hr pIC

count <- readRDS("results/RNASeq_RAWCount_count.rds")
head(count)

# select genes whose sum of expression across all the condition is greater than zero
countdata <- count[rowSums(count)>1,]
head(countdata)
dim(countdata)


replicate <- factor(c(rep(1,1),rep(2,1)))
condition <- gsub("_R[0-9]","",colnames(countdata))
coldata <- data.frame(row.names=colnames(countdata),condition)

# select only unstimulated and 2hr and 6hr pIC samples (Control and NCoR1 KD)
countdata = countdata[,c(1,2,11,12,7,8,17,18,23,24,27,28)]
coldata = coldata[c(1,2,11,12,7,8,17,18,23,24,27,28),]
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)
norm_count <- as.data.frame(counts(dds, normalized=TRUE))

vst <- varianceStabilizingTransformation(dds, blind=FALSE)
# PCA plot
plotPCA(vst , intgroup=c('condition'))
  
############################################################################

# Extract DEGs of comparison between 2hr NCoR1 KD pIC and 2hr Control pIC samples
NP2_vs_EP2 <- results(dds, contrast=c("condition","NCoR1_2hr_pIC","Emp_2hr_pIC"))%>% as.data.frame()
NP2_vs_EP2 <- merge(NP2_vs_EP2, as.data.frame(counts(dds, normalized=TRUE)[,c(1,2,7,8)]), by="row.names", sort=FALSE)
names(NP2_vs_EP2)[1] <- "Gene"

twohr_pIC.Volcano.plot =function(df,gene_label,Title){
  df <- df[which(df$log2FoldChange != "NA" & df$padj != "NA"),]
  df <- df %>%
    mutate(reg = case_when(
      df$log2FoldChange >= 1 & df$padj <= 0.05 ~ "Up",
      df$log2FoldChange <= -1 & df$padj <= 0.05 ~ "Down",
      abs(df$log2FoldChange) <= 1 & df$padj >= 0.05 ~ "No Change",
      abs(df$log2FoldChange) <= 1 & df$padj <= 0.05 ~ "No Change",
      abs(df$log2FoldChange) > 1 & df$padj >0.05 ~ "No Change"
    )) %>%
    mutate(reg = factor(reg, 
                        levels = c("Up", "No Change","Down")))
  
  label <- gene_label
  up_label = dim(df[which(df$reg =="Up"),])[1]
  down_label = dim(df[which(df$reg =="Down"),])[1]
  
  data= subset(df, Gene %in% label)
  data= data[which(abs(data$log2FoldChange) >=1),]
  data= data[which(abs(data$padj) <=0.05),]
  df.plt <- ggplot(df,aes(x=log2FoldChange,y=-log10(padj),label = Gene))+
    geom_point(aes(color=reg),size=0.2)+
    scale_color_manual(name = "Differential \n regulation",
                       values = c("Down" = "blue",
                                  "No Change" = "grey",
                                  "Up" = "red"))+
    theme_bw()+
    #ylim(0,max(-log10(df$padj)))+ 
    xlim(min(df$log2FoldChange),max(df$log2FoldChange))+
    xlab("log2 (Fold Change)")+
    ylab("-log10(adj p-value)")+
    gg_theme+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    geom_text_repel(
      data          = data,
      size          = 3.5,
      direction    = "y",
      angle        = 0,
      hjust=1,
      vjust        = 0.5,
      #segment.size = 0.2,
      #fill = data$stat1,
      #arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
      nudge_x = ifelse(data$log2FoldChange >1 , data$log2FoldChange,  -2+data$log2FoldChange) ,
      nudge_y = 20 + -log10(data$padj),
    segment.size  = 0.5,
    segment.color = "grey",
    max.overlaps = 1000    # max.overlaps = 1000
    ) +
  geom_vline(xintercept=c(-1,1), linetype="dashed",size=0.5)+
    geom_hline(yintercept=c(1.3), linetype="dashed",size=0.5)
  df.plt <- df.plt +annotate("text", x = -3, y = 78, label = down_label,color="blue",size=6)
  df.plt + annotate("text", x = 3, y = 78, label = up_label,color="red",size=6) +
    ggtitle(Title)
}

# Generate volcano plot of 2hr pIC comparison (NCoR1 KD vs Control)
pdf("plots/NCoR1_2hr_pI.volcano_plot.pdf",height = 7,width = 8)
twohr_pIC.Volcano.plot(NP2_vs_EP2,labRow,"NCoR1 KD vs Control\n(2hr pIC)")
dev.off()

#################################################################
# gene expression box plot in uns 2hr and 6hr pIC (NCoR1 KD and Control) for pIC specific genes)
twohr_sixhr_pIC_vst = merge_rep(vst) %>% rownames_to_column("Gene") %>% melt()
twohr_sixhr_pIC_vst$Stimulation = gsub("Emp_|NCoR1_","",twohr_sixhr_pIC_vst$variable)
twohr_sixhr_pIC_vst$Group = gsub("_.*","",twohr_sixhr_pIC_vst$variable)
twohr_sixhr_pIC_vst$Stimulation = factor(twohr_sixhr_pIC_vst$Stimulation,levels = c("Uns","2hr_pIC","6hr_pIC"))
pdf("plots/2hr_pIC_expressiom_boxplot.pdf",height = 4.5,width = 5)
twohr_sixhr_pIC_vst %>% 
        filter(Gene %in% pIC_specific_genes) %>%
        ggplot(.,aes(x=Stimulation,y=value,colour=Group))+
        geom_boxplot()+
        scale_colour_manual(values = c("blue","red"),labels = c("Control","NCoR1 KD"))+
        ylab ("vst (DESeq2)")+
        gg_theme_xrot+
        stat_compare_means(aes(group = Group), label = "p.signif",  label.x = 1.5,method = "wilcox",paired = FALSE)
dev.off()
###############################################################################################################


