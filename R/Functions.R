# H3K27ac density plot (Control)
H3k27ac_density = function(df,column,Title){
  dat1 <- read.csv(df,sep = "\t",header =T,na.strings = "NA")
  df1 <- dat1[,column]
  colname <- gsub("...NCoR1_H3K27ac_analysis.tag_dir.|...NCoR1_H3K27ac_Rep2_analysis.tag_dir.","",colnames(df1))

  colname <- gsub("..Coverage","",colname)
  colname[1] <- "Dist"
  colnames(df1) <- colname
  colnames(df1)[8] = "Emp_6hr_CpG_pIC"
  colnames(df1)[3] = "Emp_0hr_Rep2"
  df2 <- reshape2::melt(df1,id="Dist",variable.name = "Stimulation")
  df2 <- reshape2::melt(df1[,c(1,2,4,6,8)],id="Dist",variable.name = "Stimulation")
  df2$condition   =gsub("Emp_|NCoR1_|_Rep2","",df2$Stimulation)
  tgc <- summarySE(df2, measurevar="value", groupvars=c("Dist","condition"))
  tgc$condition = factor(tgc$condition,levels = c("0hr","6hr_CpG","6hr_pIC","6hr_CpG_pIC"),labels = c("Uns","CpG","pIC","CpG+pIC"))
  df2$condition = factor(df2$condition,levels = c("0hr","6hr_CpG","6hr_pIC","6hr_CpG_pIC"),labels = c("Uns","CpG","pIC","CpG+pIC"))
  
  den <- tgc %>%  ggplot(aes(x=Dist  ,y=value,color=condition))+
          #geom_boxplot()+
          geom_line(size=1) +
          gg_theme+
          #facet_wrap(~condition,ncol = 4)+
          theme_bw()+
          geom_ribbon(aes(ymax = value + sd, ymin = value - sd),
            alpha = 0.5,
            fill = "grey70",
            colour=NA)+

          #ylim(0,8)+
          #scale_color_manual(values = c("#f7766d","#38bf7d","#3daff5","#e789f0"))+
          gg_theme +
          labs(x="Distance from NCoR1 peaks center (bp)",y="Normalized Tag Count")
          ggtitle(Title)
  return(den1)
}


# H3K27ac (Control and NCoR1)
H3k27ac_density = function(df,column,Title){
  dat1 <- read.csv(df,sep = "\t",header =T,na.strings = "NA")
  df1 <- dat1[,column]
  colname <- gsub("...NCoR1_H3K27ac_analysis.tag_dir.","",colnames(df1))
  colname <- gsub("..Coverage","",colname)
  colname[1] <- "Dist"
  colnames(df1) <- colname
  
  df2 <- reshape2::melt(df1,id="Dist",variable.name = "Stimulation")
  df2$condition   =gsub("Emp_|NCoR1_","",df2$Stimulation)
  df2$Genotype = gsub("_[0-9].*","",df2$Stimulation)
  df2$condition = factor(df2$condition,levels = c("0hr","6hr_CpG","6hr_pIC","6h_CpG_pIC"),labels = c("Uns","CpG","pIC","CpG+pIC"))
  
  den1 <- df2 %>%  ggplot(aes(x=Dist  ,y=value,color=Genotype))+
    #geom_boxplot()+
    geom_line(size=1) +
    gg_theme+
    facet_wrap(~condition,ncol = 4)+
    theme_bw()+
    #ylim(0,8)+
    scale_color_manual(values = c("blue","red"))+
    gg_theme +
    labs(x="Distance from NCoR1 peaks center (bp)",y="Normalized Tag Count")+
    ggtitle(Title)
  return(den1)
}

# extract peak count from peakAnnolist
annotation = function(x){
  as.data.frame(x) %>% 
    dplyr::select(annotation) %>%
    separate(annotation, c("A","B"), " \\(",extra = "drop",fill = "right") %>% 
    dplyr::select(A) %>% table()  %>% as.data.frame() 
}


# Volcano plot 
library(ggrepel)
Volcano.plot =function(df,gene_label,Title){
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
            geom_point(aes(color=reg),size=0.5)+
            scale_color_manual(name = "Differential \n regulation",
                               values = c("Down" = "blue",
                                          "No Change" = "grey",
                                          "Up" = "red"))+
            theme_bw()+
            #ylim(0,max(-log10(df$padj)))+ 
            xlim(-6,11)+
            xlab("log2 (Fold Change)")+
            ylab("-log10(adj p-value)")+
            gg_theme+
            guides(colour = guide_legend(override.aes = list(size=3)))+
            #geom_text_repel(
            #  data          = data,
            #  size          = 5,
            #  #nudge_y      = 30,
            #  direction    = "x",
            #  angle        = 0,
            #  vjust        = 0,
            #  #segment.size = 0.2,
            #  # #angle         = 45,
            #  # #fill = data$stat1,
            #  #arrow = arrow(length = unit(0.01, "npc"), type = "closed", ends = "first"),
            #  nudge_y       = 5 + data$log2FoldChange ,
            #  segment.size  = 0.5,
            #  segment.color = "black"
            #  # direction     = "x"
            #) +
            geom_vline(xintercept=c(-1,1), linetype="dashed",size=0.5)+
            geom_hline(yintercept=c(1.3), linetype="dashed",size=0.5)
  df.plt <- df.plt +annotate("text", x = -4, y = max(-log10(df$padj)), label = down_label,color="blue",size=6)
  df.plt + annotate("text", x = 8, y = max(-log10(df$padj)), label = up_label,color="red",size=6) +
           ggtitle(Title)
}


plot_h3k27ac_Emp = function(df,x,y){
  plt = merge(df,H3K27ac_vst_df,by=0) %>% 
    #filter(log2FoldChange >=1 | log2FoldChange <= -1) %>%
    ggplot(.,aes_string(x=x,y=y))+
    #geom_hex(bins = 70) +
    geom_point(size=0.0001,aes(color=cut(log2FoldChange, c(-6, -1, 1, 6)))) +
    scale_color_manual(values = c("blue","grey","red"),label=c("Down","No Change","Up")) +
    theme_scatter #+guides(colour = guide_legend(override.aes = list(size=10)))
  plt = plt + annotate("text", x = 6, y = 10.5, label = dim(df[which(df$log2FoldChange >=1),])[1],color="red",size=6)
  plt = plt + annotate("text", x = 10.5, y = 6, label = dim(df[which(df$log2FoldChange <= -1),])[1],color="blue",size=6) 
  return(plt)
}


H3K27ac.Volcano.plot =function(df,Title){
  df <- df[which(df$log2FoldChange != "NA" & df$padj != "NA"),]
  df <- df %>%
    mutate(reg = case_when(
      df$log2FoldChange >= 1 & df$padj < 0.05 ~ "Up",
      df$log2FoldChange <= -1 & df$padj < 0.05 ~ "Down",
      abs(df$log2FoldChange) < 1 & df$padj >= 0.05 ~ "No Change",
      abs(df$log2FoldChange) < 1 & df$padj <= 0.05 ~ "No Change",
      abs(df$log2FoldChange) > 1 & df$padj >0.05 ~ "No Change"
    )) %>%
    mutate(reg = factor(reg, 
                        levels = c("Up", "No Change","Down"))) #%>% 
    #filter(reg != "No Change")
  
  up_label = dim(df[which(df$reg =="Up"),])[1]
  down_label = dim(df[which(df$reg =="Down"),])[1]
  
  df.plt <- ggplot(df,aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(aes(color=reg),size=0.08)+
    scale_color_manual(name = "Differential \n regulation",
                       values = c("Down" = "blue",
                                  "No Change" = "grey",
                                  "Up" = "red"))+
    theme_bw()+
    #ylim(0,max(-log10(df$padj)))+ 
    ylim(0,70)+ 
    #xlim(min(df$log2FoldChange),max(df$log2FoldChange))+
    xlim(-5,5)+
    xlab("log2 (Fold Change)")+
    ylab("-log10(adj p-value)")+
    gg_theme+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    geom_vline(xintercept=c(-1,1), linetype="dashed",size=0.5)+
    geom_hline(yintercept=c(1.3), linetype="dashed",size=0.5)
  df.plt <- df.plt +annotate("text", x = -4.5, y = 65, label = down_label,color="blue",size=6)
  df.plt + annotate("text", x = 4.5, y = 65, label = up_label,color="red",size=6) +
    ggtitle(Title)
}

KD.H3K27ac.Volcano.plot =function(df,Title){
  df <- df[which(df$log2FoldChange != "NA" & df$padj != "NA"),]
  df <- df %>%
    mutate(reg = case_when(
      df$log2FoldChange >= 1 & df$padj < 0.05 ~ "Up",
      df$log2FoldChange <= -1 & df$padj < 0.05 ~ "Down",
      abs(df$log2FoldChange) < 1 & df$padj >= 0.05 ~ "No Change",
      abs(df$log2FoldChange) < 1 & df$padj <= 0.05 ~ "No Change",
      abs(df$log2FoldChange) > 1 & df$padj >0.05 ~ "No Change"
    )) %>%
    mutate(reg = factor(reg, 
                        levels = c("Up", "No Change","Down"))) #%>% 
  #filter(reg != "No Change")
  
  up_label = dim(df[which(df$reg =="Up"),])[1]
  down_label = dim(df[which(df$reg =="Down"),])[1]
  
  df.plt <- ggplot(df,aes(x=log2FoldChange,y=-log10(padj)))+
    geom_point(aes(color=reg),size=0.08)+
    scale_color_manual(name = "Differential \n regulation",
                       values = c("Down" = "blue",
                                  "No Change" = "grey",
                                  "Up" = "red"))+
    theme_bw()+
    #ylim(0,max(-log10(df$padj)))+ 
    ylim(0,15)+ 
    #xlim(min(df$log2FoldChange),max(df$log2FoldChange))+
    xlim(-3,5)+
    xlab("log2 (Fold Change)")+
    ylab("-log10(adj p-value)")+
    gg_theme+
    guides(colour = guide_legend(override.aes = list(size=3)))+
    geom_vline(xintercept=c(-1,1), linetype="dashed",size=0.5)+
    geom_hline(yintercept=c(1.3), linetype="dashed",size=0.5)
  df.plt <- df.plt +annotate("text", x = -2.5, y =14, label = down_label,color="blue",size=6)
  df.plt + annotate("text", x = 2.5, y = 14, label = up_label,color="red",size=6) +
    ggtitle(Title)
}


###############################################################################
# Pathway enrichemnt analysis (fgsea)
library(fgsea)
DE_fgsea = function(df,rows){
  df <- df[which(abs(df$log2FoldChange) >=1 & df$padj <= 0.05 ),]
  df$fcsign <- sign(df$log2FoldChange)
  df$logP=-log10(df$padj)
  df$metric= df$logP/df$fcsign
  if("Gene" %in% colnames(df)){ 
    df$Gene <- toupper(df$Gene)
    }
  df <-df[,c("Gene", "metric")]
  ranks <- deframe(df)
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- ranks[which(ranks != "Inf")]
  # head(ranks, 20)
  pathways.hallmark <- gmtPathways("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
  #pathways.hallmark %>%   head() %>%   lapply(head)
  register(SerialParam())
  fgseaRes <- fgsea(pathways=pathways.hallmark, stats=ranks, nperm=10000)
  fgseaResTidy <- fgseaRes %>%   as_tibble() %>%   arrange(desc(NES)) %>% as.data.frame() %>% filter(pval<0.05)
  rownames(fgseaResTidy) = seq(length=nrow(fgseaResTidy))
  # Show in a nice table:
  #fgseaResTidy %>% 
  #  dplyr::select(-ES, -nMoreExtreme) %>% 
  #  filter(pval <=0.05) %>% View()
  #  arrange(padj) %>% 
  #  DT::datatable()
  fgseaResTidy %>% #View()
    slice(rows) %>%
    mutate(pathway = gsub("REACTOME|_"," ",pathway),
           reg=case_when(NES >0 ~ "Up",
                         NES < 0 ~ "Down")) %>% 
    ggplot(aes(reorder(pathway, NES), NES,fill=reg)) +
    geom_bar(stat="identity") +
    coord_flip() +
    ylab("") +
    ylab("Normalized Enrichment Score")+
    scale_fill_manual(values = c("blue","red"))+
    scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
    scale_y_continuous(expand = c(0,0)) +
    gg_theme
}


# pathway enrichment from clusterpofileR
library(tidytext)
Gene_Set_Enrichment = function(gene_set){
  term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
  gen_set =lapply(gene_set, toupper)
  compareCluster(diff_NCoR1.annotation.df.list_toupper, fun="enricher",TERM2GENE=term2gene) %>% as.data.frame() %>%
  #write.table("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/NCoR1_diff")
  group_by(Cluster) %>% 
  slice(1:5) %>%
  ungroup %>% 
  mutate(Cluster = factor(Cluster,levels = sort(unique(Cluster))),
         Description = gsub("REACTOME|_"," ",Description)) %>%
  #Description = reorder(Description, pvalue)) 
  #as.data.frame(.) %>% 
  ggplot(aes(x=reorder(Description,-pvalue), y=-log(pvalue), fill = p.adjust)) +
  geom_bar(stat="identity") +
  #geom_text(aes(label =str_wrap(geneID, width = 130), y = 0), color = "black", hjust = 0,size=4)+
  facet_wrap(~Cluster, scales = "free",ncol = 1) +
  coord_flip() +
  
  #scale_x_reordered() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 40))+
  scale_fill_gradient(low="red",high="blue")+
  scale_y_continuous(expand = c(0,0)) +
  gg_theme
}


# plot RNAseq and H3K27ac scatter plot
plot_RNaseq_h3k27ac_scatter = function(df,x,y,SYMBOL){
  ggplot(df,aes_string(x=x,y=y))+
    geom_point(aes(color=Genes),size=0.2) +
    stat_cor()+
    geom_smooth(method="lm",se = TRUE,color="grey")+
    scale_color_manual(values = c("blue","darkgreen","red"))+
    xlim(-2,2)+ylim(-2,2)+
    geom_vline(xintercept = 0,linetype = "dashed",size = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed",size = 0.5) +
    gg_theme+
    geom_label_repel(data = subset(df,SYMBOL %in% labRow),
                    aes(label = SYMBOL,fill=factor(Genes),segment.colour = factor(Genes)),,
                    color = "white" ,fontface="italic",box.padding = 0.5,max.overlaps = 1000,
                    size = 3.5,segment.size      = 0.2)+
    scale_fill_manual(values = c("blue","darkgreen","red"),aesthetics = c("fill", "segment.color"))
    #geom_label_repel(data = subset(df,SYMBOL %in% labRow),
    #                 aes(label = SYMBOL,fill=factor(Genes),segment.colour = factor(Genes)),
    #                 color = "white" ,box.padding = 0.5,max.overlaps = 40,
    #                 size = 3.5)+
    #scale_fill_manual(values = c("blue","darkgreen","red"),aesthetics = c("color", "segment.color"))    # scale_fill_manual(values = c("blue","darkgreen","red"),aesthetics = c("color", "segment.color"))
}
  
# plot H3K27ac differential enhancer scatter plot
theme_scatter = theme_bw(20) + 
  theme(axis.text.x=element_text(size=15,color="black"),
        axis.title.x=element_text(size=15),
        axis.text.y=element_text(size=15,color="black"),
        axis.title.y=element_text(size=15),
        legend.justification=c(0,1),
        legend.title = element_blank(),
        plot.title = element_text(hjust = 0.5,size=15))


plot_h3k27ac = function(df,x,y){
  plt = merge(df,H3K27ac_Rlog.df,by=0) %>% 
    ggplot(.,aes_string(x=x,y=y))+
    geom_point(size=0.1,aes(color=cut(log2FoldChange, c(-6, -1, 1, 6)))) +
    scale_color_manual(values = c("blue","grey","red"),label=c("Down","No Change","Up")) +
    theme_scatter +guides(colour = guide_legend(override.aes = list(size=10)))
  plt = plt + annotate("text", x = 6, y = 10.5, label = dim(df[which(df$log2FoldChange >=1),])[1],color="red",size=6)
  plt = plt + annotate("text", x = 10.5, y = 6, label = dim(df[which(df$log2FoldChange <= -1),])[1],color="blue",size=6) 
  return(plt)
}


#intersection
intersection = function(gene.df){
  ncor_cluster_gene_exp.list = list()
  for (j in 1:3){
    
    #print(j)
    if( j  ==1){
      common = Reduce(intersect, list(gene.df[[j]],gene.df[[j+1]],gene.df[[j+2]]))
      diff = setdiff(gene.df[[j]],gene.df[[j+1]])
      diff = setdiff(diff,gene.df[[j+2]])
      com = Reduce(intersect, list(gene.df[[j]],gene.df[[j+1]]))
      CpG_pIC_all3_vs_CpG_all3 = setdiff(com,common)
      ncor_cluster_gene_exp.list[["CpG_pIC_all3_vs_CpG_all3"]] = CpG_pIC_all3_vs_CpG_all3
      
    }
    #print(length(diff))
    if( j ==2){
      diff = setdiff(gene.df[[j]],gene.df[[j+1]])
      diff = setdiff(diff,gene.df[[1]])
      com = Reduce(intersect, list(gene.df[[j]],gene.df[[j+1]]))
      CpG_all3_vs_Uns_pIC = setdiff(com,common)
      ncor_cluster_gene_exp.list[["CpG_all3_vs_Uns_pIC"]] = CpG_all3_vs_Uns_pIC
    }
    
    if( j ==3){
      diff = setdiff(gene.df[[j]],gene.df[[j-1]])
      diff = setdiff(diff,gene.df[[j-2]])
      com = Reduce(intersect, list(gene.df[[j]],gene.df[[j-2]]))
      Uns_pIC_vs_CpG_pIC_all3 = setdiff(com,common)
      ncor_cluster_gene_exp.list[["Uns_pIC_vs_CpG_pIC_all3"]] = Uns_pIC_vs_CpG_pIC_all3
      
    }
    
    
    ncor_cluster_gene_exp.list[[names(gene.df[j])]] = diff
    ncor_cluster_gene_exp.list[["Common"]] = common
    
    
  }
  return(ncor_cluster_gene_exp.list)
}


# Generate expression bar plot of any genes based on only gene name
norm_count_bar_plot = function(gene){
  gene_exp = Emp_NCoR1_Normalized_count[gene,] %>% reshape2::melt()
  gene_exp$variable = gsub('_R2|_R1',"",gene_exp$variable)
  gene_exp$condition = gsub('_.*',"",gene_exp$variable)
  gene_exp$stimulation = gsub('Emp_|NCoR1_',"",gene_exp$variable)
  gene_exp$condition <- factor(gene_exp$condition,levels = c("Emp","NCoR1"))
  gene_exp$stimulation <- factor(gene_exp$stimulation,levels = c("Uns","IFNg","6hr_CpG","6hr_pIC","6hr_all3"))
  
  tgc <- summarySE(gene_exp, measurevar="value", groupvars=c("condition","stimulation"))
  print(ggplot(tgc,aes(x=stimulation,y=value,fill=condition)) +
                   geom_bar(stat="identity",position=position_dodge()) +
                   facet_wrap(~stimulation,scales = "free_x",ncol=5)+
                   geom_errorbar( aes(ymin=value-sd, ymax=value+sd), width=.2,position=position_dodge(.9))+
                   scale_fill_manual(values=c("blue","red"))+
                   #scale_color_hue()+
                   theme_bw()+
                   theme(axis.text.x=element_blank(),
                         axis.text.y=element_text(size=10,colour = "black",face= "bold"),
                         axis.title.y=element_text(size=10,colour = "black",face= "bold"),
                         plot.title = element_text(size = 20, face = "italic",hjust = 0.5),
                         legend.text = element_text(size=10,colour = "black"),
                         #legend.position = "none",
                         strip.text = element_text(size = 15,face="bold")) + 
                   labs(x = "",y="Normalized count",label =FALSE,title = "") +
                   ggtitle(gene))
}


merge_rep= function(col1){
  N <- ncol(col1)
  name <- colnames(col1)
  obj <- vector("list",ncol(col1)/2)
  k=1
  
  for(i in 1:N) {
    
    if(i%%2 ==1 && i <= N)    {
      
      #print(i)
      ID <-rowMeans(col1[,c(i,i+1)]) 
      obj[[k]] <- ID
      nam <- colnames(col1)[i]
      nam <- str_replace(nam,"_R[0-9]","")
      names(obj)[k] <- nam
      names(obj[[k]]) <- rownames(col1)
      #print(k)
      k=k+1
    }
  }
  mat_merged <- as.data.frame(t(do.call(rbind, obj)))
  colnames(mat_merged) = names(obj)
  return(mat_merged)
}


H3K27ac_Rep_scatter_plot= function(col1){
  N <- ncol(col1)
  name <- colnames(col1)
  name = gsub("_1|_2","",name)
  name = unique(name)
  k=1
  plotList = list()
  for(i in 1:N) {
    
    if(i%%2 ==1 && i <= N)    {
      
      #print(i)
      tmp = col1[,c(i,i+1)]
      x= colnames(tmp)[1]
      y= colnames(tmp)[2]
      nam <- str_replace(colnames(tmp[,1]),"_1","")
      tmp$density =  get_density(tmp[,1], tmp[,2], n = 100)
      p = ggplot(tmp) +
                        geom_point(aes_string(x=x, y=y, color = "density"),size =0.1) + 
                        scale_color_viridis()+
                        ggtitle(name[k]) +
                        xlab("Replicate 1") + ylab("Replicate 2") +
                        xlim(4,12) +ylim(4,12)
        plotList[[k]] = p + stat_cor(aes_string(x=x, y=y),method = "pearson", label.x = 4.5, label.y = 11.5)
    print(k)
     k=k+1
    }
  }
  #mat_merged <- as.data.frame(t(do.call(rbind, obj)))
  #colnames(mat_merged) = names(obj)
  return(plotList)
}
