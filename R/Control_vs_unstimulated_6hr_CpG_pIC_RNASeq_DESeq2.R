#count data
# Read Raw count of gene expression from RNAseq extracted from BAM file using featureCount
count <- read.csv("data/Emp_NCoR1_all_condition_gene_count.tsv",sep = "\t",header = T,row.names = 1)
saveRDS(count,file = "results/RNASeq_RAWCount_count.rds")

# Differential expression analysis between Control stimulated and Unstimulated sample
# select genes whose sum of expression across all the condition is greater than zero
countdata <- count[rowSums(count)>1,][,c(1:10)]
head(countdata)

###########################################################################################################################################################################################
# Prepare metadata of the Samples based on column name
replicate <- factor(c(rep(1,1),rep(2,1)))
condition <- gsub("_R[0-9]","",colnames(countdata))
coldata <- data.frame(row.names=colnames(countdata),condition)

# Run DESeq2 
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds, betaPrior = FALSE)
resultsNames(dds)

# Normalized counts
normalized_count = as.data.frame(counts(dds, normalized=TRUE))
head(normalized_count)

# Variance stabilized transformation values
vst <- varianceStabilizingTransformation(dds, blind=FALSE)
# PCA plot
vst$condition = factor(vst$condition,
                       levels = c("Emp_Uns","Emp_IFNg","Emp_6hr_CpG","Emp_6hr_pIC","Emp_6hr_all3"))
                       
pdf("plots/Emp_Sample_RNA-Seq_PCA.pdf",height = 4,width = 6)
plotPCA(vst , intgroup=c('condition'),ntop=500) +
        ylim(-25,25)+
        geom_point(size=6)+
        scale_color_discrete(name = "", labels = c("Uns","IFNg","CpG","pIC","CpG+pIC+IFNg"))+
        gg_theme
dev.off()

Emp_vst = assay(vst)
saveRDS(Emp_vst,file = "results/Emp_vst.rds")

# Average of variance stabilized transformation values 
Emp_vst_rep_merged <- merge_rep(Emp_vst)
Emp_vst_rep_merged.scaled <- t(apply(Emp_vst_rep_merged, 1, scale))
colnames(Emp_vst_rep_merged.scaled) <- colnames(Emp_vst_rep_merged)

# Average of variance stabilized transformation values
normalized_count_rep_merged <- merge_rep(normalized_count)

###########################################################################################################################################################################################
# Create list of object to store the results and plots

# Create list to store DEGs of comparison between Control stimulated and unstimulated condition
Emp_DE.list = list()
# Create list to store Volcano plot object
Emp.Volcano.plot.list = list()

label <- c("Il10","Il27","Ido1","Ido2","Cd274","Ctla4","Il6","Il12b",
           "Lag3", "Pparg", "Ifnb1", "Myd88", "Akt3", "Lag3","Cd83","Il12a","Ncor1")

# Extract DEGs for 6hr IFNg stimulation vs unstimulated 
res = results(dds, contrast=c("condition","Emp_IFNg","Emp_Uns")) %>% as.data.frame()
EI_vs_EU <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(1,2,3,4)]), by="row.names", sort=FALSE)
names(EI_vs_EU)[1] <- "Gene"

# Upregulated genes after 6hr IFNg stimulation
Emp_DE.list[[1]] = EI_vs_EU %>% 
                      filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[1]= "EI_vs_EU_up"
# Downregulated genes after 6hr IFNg stimulation
Emp_DE.list[[2]] = EI_vs_EU %>% 
                      filter(log2FoldChange <= -1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[2]= "EI_vs_EU_down"

# Extract DEGs for 6hr CpG stimulation vs unstimulated
res = results(dds, contrast=c("condition","Emp_6hr_CpG","Emp_Uns")) %>% as.data.frame()
EC_vs_EU <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(1,2,5,6)]), by="row.names", sort=FALSE)
names(EC_vs_EU)[1] <- "Gene"

# Upregulated genes after 6hr CpG stimulation
Emp_DE.list[[3]] = EC_vs_EU %>% 
                      filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[3]= "EC_vs_EU_up"

# Downregulated genes after 6hr  CpG stimulation
Emp_DE.list[[4]] = EC_vs_EU %>% 
                      filter(log2FoldChange <= -1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[4]= "EC_vs_EU_down"
Emp.Volcano.plot.list[[1]] = Volcano.plot(EC_vs_EU,label,"6hr CpG vs Uns")


# Extract DEGs for 6hr pIC stimulation vs unstimulated
res = results(dds, contrast=c("condition","Emp_6hr_pIC","Emp_Uns")) %>% as.data.frame()
EP_vs_EU <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(1,2,7,8)]), by="row.names", sort=FALSE)
names(EP_vs_EU)[1] <- "Gene"

# Upregulated genes after 6hr pIC stimulation
Emp_DE.list[[5]] = EP_vs_EU %>% 
                      filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[5]= "EP_vs_EU_up"

# Downregulated genes after 6hr pIC stimulation
Emp_DE.list[[6]] = EP_vs_EU %>% 
                      filter(log2FoldChange <= -1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[6]= "EP_vs_EU_down"
Emp.Volcano.plot.list[[2]] = Volcano.plot(EP_vs_EU,label,"6hr pIC vs Uns")


# Extract DEGs for 6hr pIC+CpG+IFNg stimulation vs unstimulated
res = results(dds, contrast=c("condition","Emp_6hr_all3","Emp_Uns")) %>% as.data.frame()
EA_vs_EU <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(1,2,9,10)]), by="row.names", sort=FALSE)
names(EA_vs_EU)[1] <- "Gene"

# Upregulated genes after 6hr pIC+CpG+IFNg stimulation
Emp_DE.list[[7]] = EA_vs_EU %>% 
                      filter(log2FoldChange >=1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[7]= "EA_vs_EU_up"

# Downregulated genes after 6hr pIC+CpG+IFNg stimulation
Emp_DE.list[[8]] = EA_vs_EU %>% 
                      filter(log2FoldChange <= -1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[8]= "EA_vs_EU_down"
Emp.Volcano.plot.list[[3]] = Volcano.plot(EA_vs_EU,label,"6hr CpG+pIC+IFNg vs Uns")

saveRDS(Emp_DE.list, "results/DEG_Control.Stimulated_vs_Unstimulated.rds")
###########################################################################################################################################################################################
pdf("plots/Emp_DE_genes_volcano_plot.pdf",height = 5,width = 18)
do.call(grid.arrange, c(Emp.Volcano.plot.list, ncol=3))
dev.off()
###########################################################################################################################################################################################

# Extract DEGs for between  6hr CpG and 6hr pIC to identify CpG and pIC specific genes
res = results(dds, contrast=c("condition","Emp_6hr_pIC","Emp_6hr_CpG")) %>% as.data.frame()
EP_vs_EC <- merge(res, as.data.frame(counts(dds, normalized=TRUE)[,c(7,8,5,6)]), by="row.names", sort=FALSE)
names(EP_vs_EC)[1] <- "Gene"
Emp_DE.list[[3]] = EP_vs_EC %>% 
                      filter(abs(log2FoldChange) >=1 & padj <= 0.05) %>% 
                      dplyr::select("Gene") %>% pull(.,"Gene")
names(Emp_DE.list)[3]= "EP_vs_EC"

# Write results
write.table(EP_vs_EC, file="DE_genes/Emp_6hr_CpG_vs_Emp_6hr_pIC_DE.txt", sep = "\t",quote = FALSE,row.names = FALSE)
###########################################################################################################################################################################################
# Extract CpG Specific genes
CpG_specific = EP_vs_EC %>% 
                  filter(log2FoldChange <=-1) %>% 
                  filter(Gene %in% EC_vs_EU[which(EC_vs_EU$log2FoldChange>=1 & EC_vs_EU$padj <=0.05),]$Gene | 
                         Gene %in% EP_vs_EU[which(EP_vs_EU$log2FoldChange>=1 & EP_vs_EU$padj <=0.05),]$Gene )
CpG_specific_genes = CpG_specific$Gene

# pIC specific genes
pIC_specific = EP_vs_EC %>% 
                  filter(log2FoldChange >=1) %>% 
                  filter(Gene %in% EC_vs_EU[which(EC_vs_EU$log2FoldChange>=1 & EC_vs_EU$padj <=0.05),]$Gene | 
                         Gene %in% EP_vs_EU[which(EP_vs_EU$log2FoldChange>=1 & EP_vs_EU$padj <=0.05),]$Gene )
pIC_specific_genes = pIC_specific$Gene

# CpG and pIC common genes
CpG_pIC_common = EP_vs_EC %>% 
                  filter(abs(log2FoldChange) <=0.58) %>% 
                  filter(Gene %in% EC_vs_EU[which(EC_vs_EU$log2FoldChange>=1 & EC_vs_EU$padj <=0.05),]$Gene | 
                         Gene %in% EP_vs_EU[which(EP_vs_EU$log2FoldChange>=1 & EP_vs_EU$padj <=0.05),]$Gene )
CpG_pIC_common_genes = CpG_pIC_common$Gene

CpG_pIC_specific_genes = list(CpG_specific_genes,CpG_pIC_common_genes,pIC_specific_genes)
names(CpG_pIC_specific_genes) = c("CpG","CpG_pIC","pIC")
saveRDS(CpG_pIC_specific_genes, file ="results/CpG_pIC_specific_genes.rds")

for (i in c(1:length(CpG_pIC_specific_genes))){
   write.xlsx(CpG_pIC_specific_genes[i], file="results/CpG_pIC_Specific_genes.xlsx", 
              sheetName=names(CpG_pIC_specific_genes)[i], append=T,,row.names = F)
}
###########################################################################################################################################################################################
###########################################################################################################################################################################################

# Generate Heatmap of CpG and pIC specific genes Cluster
CpG_pIC_specific_gene_exp.scaled = Emp_vst_rep_merged.scaled[c(CpG_specific_genes,
                                                               pIC_specific_genes,
                                                               CpG_pIC_common_genes),] %>% as.data.frame()
CpG_pIC_specific_gene_exp.scaled = CpG_pIC_specific_gene_exp.scaled %>%
                                    rownames_to_column("Gene") %>%
                                    mutate(Group = case_when(Gene %in% CpG_specific$Gene ~ "CpG",
                                                         Gene %in% pIC_specific$Gene ~ "pIC",
                                                         Gene %in% CpG_pIC_common$Gene ~ "CpG_pIC"))
row_split = CpG_pIC_specific_gene_exp.scaled$Group

labRow = c("C1qa","C1qb","C1qc","C2","Ccl3","Ccl4","Ccl6","Ccl9","Ccr4","H2Q-10","Il12a" ,"Il1b","Socs3","Il10", # CpG response genes
           "Ccr7","Il12b","Il6","Ebi3","Cd40","Cd44","Il1a", # Common genes
           "Irf7","Il27","Ifit1","Ifit3","Ifit2","Ifnb1","Cxcl9", # pIC specific
           "Cxcl10","Ccl5","Ddx58","Il15","Isg15","Isg20","Stat1","Oasl1") # pIC specific
subset = match(labRow,CpG_pIC_specific_gene_exp.scaled$Gene)
label =  CpG_pIC_specific_gene_exp.scaled$Gene[c(subset)]
ra <- anno_mark(at = subset, label,which = "row",labels_gp = gpar(fontsize=8))
CpG_pIC_specific_hmap = Heatmap(CpG_pIC_specific_gene_exp.scaled[,c(2:6)],km = 1,
              name="Z-score",
              col= colorRamp2(c(-1.5,0,1.5),c("blue","white","red")),
              heatmap_legend_param=list(at=c(-1.5,0,1.5),color_bar="continuous", 
                                        legend_direction="horizontal", legend_width=unit(5,"cm"),
                                        title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
              
              #Split heatmap rows by Cluster
              split = row_split,
              border=TRUE,
              
              #Row annotation configurations
              cluster_rows=TRUE,
              cluster_row_slices = FALSE,
              show_row_dend=TRUE,
              #row_title="",
              row_title_side="left",
              row_title_gp=gpar(fontsize=8),
              #row_names_gp = gpar(fontsize = 2),
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
              column_names_gp = gpar(fontsize = 10, fontface="bold"),
              column_names_rot = 45, 
              column_names_side = "top",
              show_column_names=TRUE,
              
              #Dendrogram configurations: columns
              clustering_distance_columns="euclidean",
              clustering_method_columns="single",
              column_dend_height=unit(10,"mm"),
              
              #Dendrogram configurations: rows
              clustering_distance_rows="euclidean",
              clustering_method_rows="single",
              row_dend_width=unit(10,"mm"))

pdf("plots/Emp_CpG_pIC_genes_exp_heatmap.pdf",height = 6,width = 3)
draw(CpG_pIC_specific_hmap+rowAnnotation(mark=ra) ,heatmap_legend_side="bottom")
dev.off()
###########################################################################################################################################################################################
# Generate Boxplot of CpG and pIC specific genes Cluster

# create list of comparsion for which statistical significance to be calculated
my_comparisons = list(c("Emp_6hr_CpG","Emp_6hr_all3"),
                      c("Emp_6hr_CpG","Emp_6hr_pIC"),
                      c("Emp_6hr_pIC","Emp_6hr_all3"))

facet_label = c(paste0("CpG(",length(CpG_specific_genes),")"),
                paste0("CpG & pIC(",length(CpG_pIC_common_genes),")"),
                paste0("pIC( ",length(pIC_specific_genes),")"))
names(facet_label) = c("CpG","CpG_pIC","pIC")
CpG_pIC_specific_gene_exp.scaled %>% 
  column_to_rownames('Gene') %>% 
  melt(id.vars=c("Group")) %>%    
  ggplot(aes(x=variable,y=value,color=variable)) + 
         geom_boxplot()+
         gg_theme +
         facet_wrap(~Group,ncol = 1,labeller = labeller(Group = facet_label))+
         theme(axis.text.x = element_blank())+
         ylab("scaled vst (DESeq2)")+
         xlab("")+
         #ylim(0,5)+
         scale_fill_discrete(name = "", labels = c("Uns","IFNg","CpG","pIC","CpG+pIC+IFNg"))+
         #stat_summary(fun.y=median, geom="line", aes(group=1), lwd=1,linetype=2)+
         stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",
                            y_position = 12,method = "wilcox",paired = FALSE,size=3)
###########################################################################################################################################################################################
# Pathway enrichment analysis for CpG and pIC specific genes

# Convert all the gene name to upper case
CpG_pIC_specific_genes_toupper =lapply(CpG_pIC_specific_genes, toupper) 


term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
compEnricher =compareCluster(CpG_pIC_specific_genes_toupper, fun="enricher",TERM2GENE=term2gene)
write.table(file = "results/Emp_LRT_clusters_pathway.csv",
            compEnricher,sep = "\t",quote = FALSE,row.names = FALSE)

# Pathway barplot
compEnricher %>%
  as.data.frame() %>% View()
  mutate(Description = gsub("REACTOME|_"," ",Description)) %>% 
  group_by(Cluster) %>%
  slice(1:6) %>%
  ungroup %>%
  as.data.frame(.) %>% 
  mutate(Cluster= factor(Cluster,levels = c("CpG","CpG_pIC","pIC"))) %>%

  ggplot(aes(x=reorder(Description,-pvalue), y=-log(pvalue), fill = p.adjust)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_wrap(~Cluster, scales = "free",ncol = 1) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 38))+
  scale_fill_gradient(low="red",high="blue")+
  scale_y_continuous(expand = c(0,0)) +
  gg_theme

###########################################################################################################################################################################################

# Identify Synergy genes
synergy = normalized_count_rep_merged %>% 
          rownames_to_column('gene') %>% filter(gene %in% unlist(CpG_pIC_specific_genes)) %>%
          dplyr::mutate(Synergy = Emp_6hr_all3/rowSums(.[,c(4,5)])) %>% 
          dplyr::filter(Synergy != "Inf" & Synergy != "NaN" & Synergy >= 1.2) %>% 
          column_to_rownames('gene') %>%
          dplyr::select(c(1:5)) %>%
          t() %>%
          scale(.) %>% 
          t() 
          
synergy_genes = rownames(synergy)
write.xlsx(synergy_genes, file="results/CpG_pIC_Specific_genes.xlsx", sheetName="Synergy Genes", append=T,row.names = F,col.names = F)

# Write list of synergy genes as text file 
write.table(synergy_genes,file = "results/Emp_Synergy_genes.txt",quote = FALSE,row.names = FALSE)

# Plot Heatmap of Synergy genes
labRow = c("Il1b","Il12b","Il6","Il10","Irf4","Cd274","Il12a","Ebi3","Socs3","Il23a","Il27","Ahr",
           "Ccr7","Cxcr5","Tnf","Ifnb1")
subset = match(labRow,rownames(synergy))
label =rownames(synergy)[c(subset)]
ra <- anno_mark(at = subset, label,which = "row",labels_gp = gpar(fontsize=10))

# Generate heatmap of Synergy genes
pdf("plots/Emp_Synergy.pdf",height = 4,width = 2.7)
draw(Heatmap(Emp_vst_rep_merged.scaled[rownames(synergy),],cluster_columns = FALSE,name = "Z-score",
        show_row_names = FALSE,
        column_names_rot = 45, 
        column_names_side = "top")+
        rowAnnotation(mark=ra),heatmap_legend_side="bottom")
dev.off()

# Plot Boxplot of normalized expression of Synergy genes
my_comparisons = list(c("Emp_6hr_CpG","Emp_6hr_all3"),c("Emp_6hr_pIC","Emp_6hr_all3"))
pdf("plots/Emp_Synergy_boxplot.pdf",height = 4,width = 6)
normalized_count_rep_merged %>% 
        filter(rownames(.) %in% synergy_genes) %>%
        #as.data.frame() %>%
        melt() %>% 
        ggplot(aes(x=variable,y=log(value+1),color=variable)) + 
        geom_boxplot()+
        gg_theme +
        theme(axis.text.x = element_blank())+
        ylab("Log(Noramlized count)")+
        xlab("")+
        ylim(0,16)+
        scale_fill_discrete(name = "", labels = c("Uns","IFNg","CpG","pIC","CpG+pIC+IFNg"))+
        #stat_summary(fun.y=median, geom="line", aes(group=1), lwd=1,linetype=2)+
        stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",
                           y_position = 12,method = "wilcox",paired = FALSE,size=5)
dev.off()

# Pathway enrichment for Synergy genes against Reactome database
term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
compEnricher = enricher(toupper(synergy_genes),TERM2GENE=term2gene)
# Pathway barplot

pdf("plots/Emp_Synergy_pathway.pdf",height = 4,width = 10)
compEnricher %>%
        as.data.frame() %>% 
        mutate(Description = gsub("REACTOME|_"," ",Description)) %>% 
        #arrange(p.adjust) %>% View()
        slice(1:6) %>% #View()
        #ungroup %>%
        #.[c(-4,-5,-7,-8,-10),] %>%
        #Description = reorder(Description, pvalue)) 
        as.data.frame(.) %>% 
        ggplot(aes(x=reorder(Description,-pvalue), y=-log(pvalue), fill = p.adjust)) +
        geom_bar(stat="identity") +
        coord_flip() +
        #facet_wrap(~Cluster, scales = "free",ncol = 1) +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 50))+
        scale_fill_gradient(low="red",high="blue")+
        scale_y_continuous(expand = c(0,0)) +
        gg_theme
 dev.off()   
###########################################################################################################################################################################################
# Identify  Antagonist genes
antagonist = normalized_count_rep_merged %>% 
                  rownames_to_column('gene') %>% filter(gene %in% unlist(CpG_pIC_specific_genes)) %>% 
                  dplyr::mutate(Antagonist = Emp_6hr_all3/rowSums(.[,c(4,5)])) %>% 
                  dplyr::filter(Antagonist != "Inf" & Antagonist != "NaN" & Antagonist <= 0.5) %>%
                  column_to_rownames('gene') %>%
                  t() %>%
                  scale(.) %>% 
                  t() %>%
                  as.data.frame()
# determine how many cluster you want, I specify 2 here
km <- kmeans(antagonist,2) 
# combine the cluster with the matrix 
m.kmeans <- cbind(antagonist, km$cluster) 
 
o <- order(m.kmeans[,5],decreasing = TRUE) # order the last column
m.kmeans <- m.kmeans[o,] # order the matrix according to the order of the last column
colnames(m.kmeans)[6] = "Cluster"
split <- factor(m.kmeans[,6], levels=c(1,2),labels = c("pIC","CpG"))
labRow = c("Cxcl9","Mx1","Mx2","Stat1","Il15","Oas1a","Oas1b","Oas1c","Oas1g",
           "Oas2","Cxcl10","Ifit3","Irf7","Oas3","Ddx60","Il10","C1qc","C2","Cfb",
           "C3ar1","C1qc","C1qb","C1qa","Cd36")
subset = match(labRow,rownames(m.kmeans))
label =rownames(m.kmeans)[c(subset)]
ra <- anno_mark(at = subset, label,which = "row",labels_gp = gpar(fontsize=10))
###########################################################################################################################################################################################
pdf("plots/Emp_Antagonist.pdf",height = 6,width = 3.4)
draw(Heatmap(m.kmeans[,c(1:5)],cluster_columns = FALSE,name = "Z-score",
          split = split,
          col= colorRamp2(c(-1,0,1),c("blue","white","red")),
          heatmap_legend_param=list(at=c(-1.5,0,1.5),color_bar="continuous", 
                                    legend_direction="horizontal", legend_width=unit(5,"cm"),
                                    title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
          show_row_names = FALSE,
          column_names_rot = 45, 
          column_names_side = "top")+
  rowAnnotation(mark=ra),heatmap_legend_side="bottom")
dev.off()
###########################################################################################################################################################################################
write.table(m.kmeans,file = "results/Emp_antagonist_genes.txt",quote = FALSE,sep = "\t")

CpG_antagonoist_genes = m.kmeans[which(m.kmeans$`km$cluster` ==1),] %>% rownames()
pIC_antagonoist_genes = m.kmeans[which(m.kmeans$`km$cluster` ==2),] %>% rownames()
write.xlsx(CpG_antagonoist_genes, file="results/CpG_pIC_Specific_genes.xlsx", 
                                  sheetName="CpG Antagonist Genes", append=T,row.names = F,col.names = F)
write.xlsx(pIC_antagonoist_genes, file="results/CpG_pIC_Specific_genes.xlsx", 
                                  sheetName="pIC Antagonist Genes", append=T,row.names = F,col.names = F)
###########################################################################################################################################################################################
pdf("plots/Emp_Antagonist_boxplot.pdf",height = 6,width = 5)
m.kmeans_count %>% 
        column_to_rownames('Row.names') %>% 
        melt(id.vars=c("Cluster")) %>% 

        ggplot(aes(x=variable,y=log(value+1),color=variable)) + 
        geom_boxplot()+
        gg_theme +
        facet_wrap(~Cluster,ncol = 1)+
        theme(axis.text.x = element_blank())+
        ylab("Log(Noramlized count)")+
        xlab("")+
        ylim(0,16)+
        scale_fill_discrete(name = "", labels = c("Uns","IFNg","CpG","pIC","CpG+pIC+IFNg"))+
        #stat_summary(fun.y=median, geom="line", aes(group=1), lwd=1,linetype=2)+
        stat_compare_means(comparisons = my_comparisons,tip.length=0.01,label = "p.format",
                           y_position = 12,method = "wilcox",paired = FALSE,size=5)
dev.off()
###########################################################################################################################################################################################
# Pathway enrichment for Antagonized genes against Reactome database
gene.list = split(rownames(m.kmeans), m.kmeans[,6])
gene.list_toupper =lapply(gene.list, toupper)

term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
compEnricher =compareCluster(gene.list_toupper, fun="enricher",TERM2GENE=term2gene)

# Pathway barplot
pdf("plots/Emp_Antagonist_pathway.pdf",height = 8,width = 10)
compEnricher %>% 
      as.data.frame() %>% #View()
      group_by(Cluster) %>% 
      slice(1:5) %>%
      #ungroup %>% 
      mutate(Cluster = factor(Cluster,levels=c(1,2),
                              labels = c("Antagonize pIC genes","Antagonized CpG genes")),
             Description = gsub("REACTOME|_"," ",Description)) %>% 
             #Description = reorder(Description, pvalue)) 
            as.data.frame(.) %>% 
      ggplot(aes(x=reorder(Description,-pvalue), y=-log(pvalue), fill = p.adjust)) +
      geom_bar(stat="identity") +
      #geom_text(aes(label =str_wrap(geneID, width = 130), y = 0), color = "black", hjust = 0,size=4)+
      facet_wrap(~Cluster, scales = "free",ncol = 1) +
      coord_flip() +

      #scale_x_reordered() +
      scale_x_discrete(labels = function(x) str_wrap(x, width = 60))+
      scale_fill_gradient(low="red",high="blue")+
      scale_y_continuous(expand = c(0,0)) +
      gg_theme
dev.off()
###########################################################################################################################################################################################
save(Emp_DE.list,
     synergy_genes,
     CpG_antagonoist_genes,
     pIC_antagonoist_genes,
     CpG_pIC_specific_genes,
     file="results/Control_RNASeq.RData")
