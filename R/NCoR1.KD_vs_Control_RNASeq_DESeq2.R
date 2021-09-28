load("results/Control_RNASeq.RData")

# Load count data

count <- readRDS("results/RNASeq_RAWCount_count.rds")
head(count)

################################################################################################

# Differential expression analysis between Control and NCoR1 KD sample
# select genes whose sum of expression across all the condition is greater than zero
countdata <- count[rowSums(count)>1,][,c(1:20)]
head(countdata)
dim(countdata)

# Prepare metadata of the Samples based on column name
replicate <- factor(c(rep(1,1),rep(2,1)))
condition <- gsub("_R[0-9]","",colnames(countdata))
coldata <- data.frame(row.names=colnames(countdata),condition)

# Run DESeq2 
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

################################################################################################
# Create list of object to store the results and plots

# DEG list
DE.list = list()
# List of Genes that d significantly up or downregulated 
NCoR1_KD_No.change.list = list()
# List to contain volcano plot of DEGs
Volcano.plot.list =list()
# List to store output of Gene set enrichment analysis performed using fgsea
fgsea.list = list()
################################################################################################
label <- c("Il10","Il27","Ido1","Ido2","Cd274","Ctla4","Il6","Il12b",
           "Lag3", "Pparg", "Ifnb1", "Myd88", "Akt3", "Lag3","Cd83","Il12a","Ncor1")

# Extract DEGs for CpG stimulated condition
res = results(dds, contrast=c("condition","NCoR1_6hr_CpG","Emp_6hr_CpG")) %>% as.data.frame()
NC_vs_EC <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(5,6,15,16)]), by="row.names", sort=FALSE)
names(NC_vs_EC)[1] <- "Gene"
DE.list[[1]] = NC_vs_EC %>% 
                  filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                  dplyr::select("Gene") %>% pull(.,"Gene")
names(DE.list)[1]= "NC_vs_EC_up"
DE.list[[2]] = NC_vs_EC %>% 
                  filter(log2FoldChange <=-1 & padj <= 0.05) %>% 
                  dplyr::select("Gene") %>% pull(.,"Gene")
names(DE.list)[2]= "NC_vs_EC_down"
NCoR1_KD_No.change.list[[1]] = NC_vs_EC %>% 
                                  filter(abs(log2FoldChange) < 0.58 & padj > 0.05) %>% 
                                  dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_KD_No.change.list)[1] = "NC_vs_EC_noChange"
Volcano.plot.list[[1]] = Volcano.plot(NC_vs_EC,label,"NCoR1 KD vs Control\n(6hr CpG)")
# pathway rows to plot 
rows= c(6,9,20,26,112,113,114,115)
fgsea.list[[1]] = DE_fgsea(NC_vs_EC,rows)

# Extract DEGs for pIC stimulated condition
res = results(dds, contrast=c("condition","NCoR1_6hr_pIC","Emp_6hr_pIC")) %>% as.data.frame()
NP_vs_EP <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(7,8,17,18)]), by="row.names", sort=FALSE)
names(NP_vs_EP)[1] <- "Gene"
DE.list[[3]] = NP_vs_EP %>% 
                  filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                  dplyr::select("Gene") %>% pull(.,"Gene")
names(DE.list)[3]= "NP_vs_EP_up"
DE.list[[4]] = NP_vs_EP %>% 
                  filter(log2FoldChange <=-1 & padj <= 0.05) %>% 
                  dplyr::select("Gene") %>% pull(.,"Gene")
names(DE.list)[4]= "NP_vs_EP_down"
NCoR1_KD_No.change.list[[2]] = NP_vs_EP %>% 
                                  filter(abs(log2FoldChange) < 0.58 & padj > 0.05) %>% 
                                  dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_KD_No.change.list)[2] = "NP_vs_EP_noChange"
Volcano.plot.list[[2]] = Volcano.plot(NP_vs_EP,label,"NCoR1 KD vs Control\n(6hr pIC)")
#pathway rows to plot 
rows= c(1,2,3,4,8,120,121,122,123,124)
fgsea.list[[2]] = DE_fgsea(NP_vs_EP,rows)

# Extract DEGs for pIC+CpG+IFNg stimulated condition
res = results(dds, contrast=c("condition","NCoR1_6hr_all3","Emp_6hr_all3")) %>% as.data.frame()
NA_vs_EA <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(9,10,19,20)]), by="row.names", sort=FALSE)
names(NA_vs_EA)[1] <- "Gene"
DE.list[[5]] = NA_vs_EA %>% 
                  filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                  dplyr::select("Gene") %>% pull(.,"Gene")
names(DE.list)[5]= "NA_vs_EA_up"
DE.list[[6]] = NA_vs_EA %>%
                  filter(log2FoldChange <=-1 & padj <= 0.05) %>% 
                  dplyr::select("Gene") %>% pull(.,"Gene")
names(DE.list)[6]= "NA_vs_EA_down"
NCoR1_KD_No.change.list[[3]] = NA_vs_EA %>% 
                                  filter(abs(log2FoldChange) < 0.58 & padj > 0.05) %>% 
                                  dplyr::select("Gene") %>% pull(.,"Gene")
names(NCoR1_KD_No.change.list)[3] = "NA_vs_EA_noChange"
Volcano.plot.list[[3]] = Volcano.plot(NA_vs_EA,label,"NCoR1 KD vs Control\n(6hr pIC+CpG+IFNg)")
rows= c(8,17,34,44,73,75,76,78)
fgsea.list[[3]] = DE_fgsea(NA_vs_EA,rows)

# Generate pdf and png of volcano plot of NCoR1KD  vs Control for all the comparison 
grid.arrange(Volcano.plot.list)
png("plots/Emp_DE_genes_volcano_plot.png",width = 36,height = 12,units = "cm",res = 500)
do.call(grid.arrange, c(Volcano.plot.list, ncol=3))
dev.off()

pdf("plots/Emp_DE_genes_FGSEA_plot.pdf",width = 30,height = 6)
do.call(grid.arrange, c(fgsea.list, ncol=3))
dev.off()

saveRDS(DE.list,file = "results/DEG_NCoR1.KD_vs_Control.rds")
################################################################################################
# Compare Number of DEGs  obtained in Control stimulated vs unstimulated condition  and 
# NCOR1 KD vs Control activated condition
Comp = cbind(as.data.frame(lengths(DE.list)),as.data.frame(lengths(Emp_DE.list[c(-1,-2)])))
Comp$Comparison = gsub("_up|_down","",rownames(Comp))
Comp$Reg = gsub(".*_vs_.*_","",rownames(Comp))
colnames(Comp) = c("NCoR1 KD","Control","Comparison","Reg")
Comp.melt = melt(Comp,id.vars = c("Comparison","Reg"))
Comp.melt$Comparison = factor(Comp.melt$Comparison, 
                              levels = c("NC_vs_EC","NP_vs_EP","NA_vs_EA"),
                              labels = c("CpG","pIC","CpG+pIC"))
Comp.melt$variable = factor(Comp.melt$variable,levels = c("Control","NCoR1 KD"))
Comp.melt$Reg = factor(Comp.melt$Reg,
                       levels = c("up","down"),
                       labels = c("Up-regulated","Down-regulated"))

# Generate bar plot of the list of DEGs 
pdf("plots/Emp_NCoR1_DE_genes_comparison.pdf",height = 5,width = 7)
Comp.melt %>% #filter(Reg != "Down-regulated") %>% 
    ggplot(.,aes(x=Comparison,y=value,fill=variable))+
          scale_fill_manual(values = c("blue","red")) +
          geom_bar(stat="identity",position = "dodge")+
          facet_grid(~Reg)+
          gg_theme_xrot+
          ylab("Number of genes")
dev.off()

################################################################################################

# Generate nomralized count for all the condition 
normalized_count = as.data.frame(counts(dds, normalized=TRUE))
saveRDS(normalized_count,file = "results/Emp_NCoR1_Normalized_count.rds")
Emp_NCoR1_Normalized_count = readRDS("results/Emp_NCoR1_Normalized_count.rds")

# Genrate vst(variance stabilized transformation) values for all the samples
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
vst$condition = factor(vst$condition,
                       levels = c("Emp_Uns","Emp_IFNg","Emp_6hr_CpG","Emp_6hr_pIC","Emp_6hr_all3"))
# PCA plot
plotPCA(vst , 
        intgroup=c('condition'),ntop=500) +
        gg_theme+
        #ylim(-25,25)+
        geom_point(size=6)+
        scale_color_discrete(name = "", labels = c("Uns","IFNg","CpG","pIC","CpG+pIC+IFNg"))

Emp_NCoR1_vst =assay(vst)
saveRDS(Emp_NCoR1_vst,file = "results/Emp_NCoR1_vst.rds")

# Function to take average of replicates across all the condition
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

# Average of variance stabilized transformation values 
Emp_NCoR1_vst_rep_merged = merge_rep(Emp_NCoR1_vst)
# Scale values across the samples
Emp_NCoR1_vst_rep_merged.scaled <- t(apply(Emp_NCoR1_vst_rep_merged, 1, scale)) %>% 
                                        as.data.frame()  %>% 
                                        filter(!is.na(.))
colnames(Emp_NCoR1_vst_rep_merged.scaled) = colnames(Emp_NCoR1_vst_rep_merged)

# Average of normalized count
normalized_count_rep_merged <- merge_rep(normalized_count)
################################################################################################

# Likelihood ratio test to get genes showing significant vriablity across the condition
# Below analysis were carried out based on tutorial given at this link 
# https://hbctraining.github.io/DGE_workshop/lessons/08_DGE_LRT.html

dds_lrt <- DESeq(dds, test="LRT", reduced = ~ 1)
# Extract results
res_LRT <- results(dds_lrt)
# Subset the LRT results to return genes with padj < 0.05 and log2FoldChange >1
sig_res_LRT <- res_LRT %>%
                  data.frame() %>%
                  rownames_to_column(var="gene") %>% 
                  as_tibble() %>% 
                  filter(padj < 0.05 & abs(log2FoldChange) >1 )

# Extract gene list 
sigLRT_genes <- sig_res_LRT %>% pull(gene)

# Total number of genes identified based on  LRT test
length(sigLRT_genes)

# Subset results for faster cluster finding 
clustering_sig_genes <- sig_res_LRT %>% arrange(padj)

# Obtain vst values for significant genes
cluster_rlog <- Emp_NCoR1_vst[NCoR1_bound_genes_clusters_df$genes, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(cluster_rlog, metadata = coldata,time="condition",col = NULL)

head(clusters$df)

# Save RDS file of cluster
saveRDS(clusters,file = "results/NCoR1_LRT_cluster_genes.rds")

NCoR1_LRT_cluster = readRDS("results/NCoR1_LRT_cluster_genes.rds")
NCoR1_LRT_cluster_df = NCoR1_LRT_cluster$df
write.table(NCoR1_LRT_cluster_df$df,file = "results/NCoR1_LRT_cluster_genes.csv",sep = "\t",row.names = FALSE,quote = FALSE)

# Generate Expression boxplot for each cluster. 
# Below mentioned code has been implemented with slight modification from 
# https://github.com/lpantano/DEGreport/blob/master/R/clustering.R

table = NCoR1_LRT_cluster[["normalized"]]

time = coldata$condition
color = NULL
min_genes = 10
process = FALSE
points = TRUE
boxes = TRUE
smooth = TRUE
lines = TRUE
facet = TRUE
cluster_column = "cluster"
prefix_title = "Cluster "


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
                                 title = paste0(prefix_title,
                                               names(counts),
                                               "(n=" ,
                                               counts,")"),
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


splan <- length(unique(table$condition)) - 1L

table$condition = factor(table$condition,levels = c("Emp_Uns","NCoR1_Uns",
                                                    "Emp_IFNg","NCoR1_IFNg",
                                                    "Emp_6hr_CpG","NCoR1_6hr_CpG",
                                                    "Emp_6hr_pIC","NCoR1_6hr_pIC",
                                                    "Emp_6hr_all3", "NCoR1_6hr_all3"))

pdf("plots/NCoR1_LRT_clusters_1_to_11_boxplot.pdf",height = 30,width = 3.7)
table %>% filter(cluster %in% c(1:11)) %>% 
  mutate(cluster=as.character(cluster)) %>%
  mutate(cluster=factor(cluster,levels = c("1","7","11","2","3","6","10","9","4","5","8"))) %>% 
  ggplot(aes(x = condition, y = value,color = condition))+
  facet_wrap(~cluster+title,ncol = 1)+ 
  geom_boxplot(alpha = 0, outlier.size = 0,outlier.shape = NA)+ 
  #geom_point(alpha = 0.4, size = 0.1,position = position_jitterdodge(dodge.width = 0.9))+
  stat_smooth(aes_string(x = "condition", y = "value",group = color),
              se = FALSE,method = "lm", formula = y~poly(x, splan))+ 
  #geom_line(aes_string(group = "line_group"), alpha = 0.1)+ 
  #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ylab("Z-score of gene abundance") +
  xlab("")+
  gg_theme+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "bottom",
        legend.title = element_blank()) 
dev.off()
#ggplot(table, aes_string(x = "condition", y = "value",fill = color, color = color))+
#      facet_wrap(~title,ncol=1)+ 
#      geom_boxplot(alpha = 0, outlier.size = 0,outlier.shape = NA)+ 
#      #geom_point(alpha = 0.4, size = 1,position = position_jitterdodge(dodge.width = 0.9))+ 
#      #stat_smooth(aes_string(x = "condition", y = "value",
#       #                      group = color, color = color),
#       #           se = FALSE,
#       #           method = "lm", formula = y~poly(x, splan))+ 
#      geom_line(aes_string(group = "line_group"), alpha = 0.1)+ 
#      #theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#      ylab("Z-score of gene abundance") +
#      xlab("")+
#      gg_theme+
#      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#      scale_colour_manual(values = c("grey"))

################################################################################################
# Pathway Enrichment analysis for each Cluster 

# Create list of genes from each cluster
gene.list = split(table$genes, table$cluster)
# Convert all the gene name to upper case
gene.list_toupper =lapply(gene.list, toupper) 

# Read .gmt file of reactome pathway downlaoded from http://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#H
term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
compEnricher =compareCluster(gene.list_toupper, fun="enricher",TERM2GENE=term2gene)
write.table(file = "results/NCoR1_LRT_clusters_pathway.csv",
            compEnricher,sep = "\t",quote = FALSE,row.names = FALSE)

# Pathway barplot
compEnricher %>%
  as.data.frame() %>% 
  group_by(Cluster) %>%
  slice(1:6) %>%
  ungroup %>% 
  mutate(Description = reorder_within(Description, -pvalue, Cluster)) %>% 
  as.data.frame() %>% 
  ggplot(aes(x=Description, y=-log(pvalue), fill = p.adjust)) +
         geom_bar(stat="identity") +
         coord_flip() +
         facet_wrap(~Cluster, scales = "free",ncol = 1) +
         scale_x_reordered() +
         #scale_x_discrete(labels = function(x) str_wrap(x, width = 35))+
         scale_fill_gradient(low="red",high="blue")+
         scale_y_continuous(expand = c(0,0)) +
         gg_theme

################################################################################################
# Generate Heatmap of Cluster 1 to 11
NCoR1_LRT_cluster1_11 = NCoR1_LRT_cluster_df %>% filter(cluster %in% c(1:11))

NCoR1_LRT_cluster1_11_vst = Emp_NCoR1_vst_rep_merged.scaled[NCoR1_LRT_cluster1_11$genes,] %>% filter(!is.na(.))
NCoR1_LRT_cluster1_11 = NCoR1_LRT_cluster1_11[rownames(NCoR1_LRT_cluster1_11_vst),]

# split cluster
split <- factor(paste0("Cluster", NCoR1_LRT_cluster1_11$cluster), 
                  levels=c("Cluster1","Cluster7","Cluster11","Cluster2",
                           "Cluster3","Cluster6","Cluster10","Cluster9",
                           "Cluster4","Cluster5","Cluster8"))

NCoR1_LRT_cluster1_11_hmap <- Heatmap(NCoR1_LRT_cluster1_11_rlog[,c(1,6,2,7,3,8,4,9,5,10)],km = 1,
              name="Z-score",
              col= colorRamp2(c(-1.5,0,1.5),c("blue","white","red")),
              heatmap_legend_param=list(at=c(-1.5,0,1.5),color_bar="continuous", 
                                        legend_direction="horizontal", 
                                        legend_width=unit(5,"cm"),
                                        title_position="topcenter", 
                                        title_gp=gpar(fontsize=10, 
                                        fontface="bold")),
              
              #Split heatmap rows by gene family
              split = split,
              border=TRUE,
              # box colour
              #rect_gp = gpar(col = "black"),              
              #Row annotation configurations
              cluster_rows=TRUE,
              cluster_row_slices = FALSE,
              show_row_dend=FALSE,
              row_title_side="left",
              row_title_gp=gpar(fontsize=8),
              show_row_names=FALSE,
              row_names_side="right",
              #row_title_rot=30,

              #Column annotation configuratiions
              cluster_columns=FALSE,
              show_column_dend=TRUE,
              column_title="Samples",
              column_title_side="top",
              column_title_gp=gpar(fontsize=10, fontface="bold"),
              #column_title_rot=45,
              column_names_gp = gpar(fontsize = 10),
              
              show_column_names=TRUE,
              
              #Dendrogram configurations: columns
              clustering_distance_columns="euclidean",
              clustering_method_columns="single",
              column_dend_height=unit(10,"mm"),
              
              #Dendrogram configurations: rows
              clustering_distance_rows="euclidean",
              clustering_method_rows="single",
              row_dend_width=unit(10,"mm"))

pdf("plots/NCoR1_LRT_clusters_1_to_11.pdf",height = 9,width = 5)
draw(NCoR1_LRT_cluster1_11_hmap ,heatmap_legend_side="right")
dev.off()

################################################################################################
save(DE.list,
     Emp_NCoR1_vst_rep_merged,
     Emp_NCoR1_vst_rep_merged.scaled,
     NCoR1_LRT_cluster1_11,
     file = "results/NCoR1.KD_RNASeq.RData")
