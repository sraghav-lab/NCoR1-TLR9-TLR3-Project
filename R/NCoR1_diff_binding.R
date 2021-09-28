################################################################################################################
load("results/Control_RNASeq.RData")

# List all the NCoR1 peaks
NCoR1_peak.list= list.files(path = "/results/NCoR1_ChIPSeq_peaks",
                            pattern = "*_NCoR.bed",full.names = TRUE)[-3][c(4,2,3,1)]

names(NCoR1_peak.list) = c("Uns","CpG","pIC","all3")
NCoR1_peakAnnoList <- lapply(NCoR1_peak.list, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000),
                       annoDb= "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))

# Plot distribution of NCoR1 peaks based on distance relative to TSS
plotDistToTSS(NCoR1_peakAnnoList)+
        gg_theme
ggsave("plots/NCoR1_peak_feature_dist.pdf",width = 16,height = 7,units = "cm",dpi=500)

NCoR1.peak.bed = lapply(NCoR1_peak.list, read.csv, sep="\t",header=F)

# Extract number of peaks in each stimulation condition
Z <- lapply(NCoR1.peak.bed, nrow) %>% unlist()
df2 <- data.frame(Condition=c("Uns","CpG","pIC","CpG+pIC+IFNg"),
                  value=as.numeric(Z))
df2$Condition <- factor(df2$Condition, levels=c("Uns","CpG","pIC","CpG+pIC+IFNg"))

# Bar plot to show number of peaks in each stimulation condition
ggplot(data=df2, aes(x=Condition, y=value, fill=Condition)) +
        geom_bar(stat = "identity") + #facet_grid(~Stimulation)+
        geom_text(aes(label=value,y=value), vjust=1.5, hjust=0.5,color = "black",position=position_stack(),size=5)+
        gg_theme+
        labs(x = "",y="No. of NCoR1 binding regions",label =FALSE,title = "")+
        scale_y_continuous(expand = c(0,0),limits = c(0,17000))
ggsave("plots/NCoR1_peak_count_Uns_CpG_pIC_all3.pdf",width = 16,height = 12,units = "cm",dpi=500)

# Read Differential NCoR1 peak 
dat = list.files(pattern = "2fold_.*density$",full.names = TRUE,
                 path = "/results/NCoR1_ChIPSeq_peaks/")

# NCoR1 peaks overlapping with H3K27ac overlap wit
head(dat)
plotlist =list()

for(i in 1:5){
        dat1 <- read.csv(dat[i],sep = "\t",header =T,na.strings = "NA")
        #df1 <- dat1[,c(1,2,5,8,11,14,17,20,23,26,29,32,35,38,41)]
        df1 <- dat1[,c(1,2,5,8,11)]
        colnames(df1) <- c("Dist","Uns","CpG","pIC","CpG+pIC+IFNg")
        head(df1)
        ptitle = gsub(".density","",basename(dat[i]))
        df2 <- melt(df1,id="Dist",variable.name = "Stimulation")
        head(df2)
        x = min(df2$value)
        y = max(df2$value)
        #pname <- paste0("Plot-12041_no_change",i)
        den1 <- df2 %>% 
                ggplot(aes(x=Dist  ,y=value,color=Stimulation))+
                geom_line(size=1) +
                gg_theme+
                #facet_wrap(~condition,ncol = 7)+
                theme_bw()+
                ylim(x,y)+
                #scale_color_manual(values = c("blue","red"))+
                theme(axis.line.x = element_line(color="black", size = 0.5),
                      plot.title = element_text(face="italic"),
                      axis.line.y = element_line(color="black", size = 0.5),
                      axis.text.y=element_text(size=10,colour = "black"),
                      axis.text.x=element_text(size=10,colour = "black"),
                      axis.title.y = element_text(size=10,colour = "black"),
                      axis.title.x = element_text(size=10,colour = "black"),
                      legend.text=element_text(size=15))+
                      #legend.position = "none")+
                ggtitle(ptitle)+
                labs(x="Distance from NCoR1 peaks summit (bp)",y="Normalized Tag count")
        plotlist[[i]] = den1
        # if( c < 15){
        #c = c +2 
        #}
}

pdf("plots/NCoR1_binding_diff_peaks.pdf",height = 14,width = 5)
do.call(grid.arrange,  c(plotlist,ncol=1))
dev.off()

NCoR1_diffPeak.list = list.files(pattern = "^CpG_pIC_all3.bed$|^CpG_all3.bed$|^Uns_pIC.bed$|^Uns_vs_CpG_pIC_all3.bed$|^no-change_uns_CpG_pIC_all3.bed$",full.names = TRUE,
                                 path = "/results/NCoR1_ChIPSeq_peaks/")
names(NCoR1_diffPeak.list) = c("CpG_all3","CpG_pIC_all3","No_Change","Uns_pIC","Uns")
diff_NCoR1.bed = lapply(NCoR1_diffPeak.list, import.bed)
resizeRanges <- lapply(diff_NCoR1.bed, resize ,width = 1000,fix = 'center')

write_bed= function(x,filename){
        bed_file = as.data.frame(x)
        write.table(bed_file[,c(1:3,6)],paste0(filename,"_2kb.bed"),
                    sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}
lapply(seq_along(resizeRanges), function(y, n, i) { write_bed(y[[i]],n[[i]]) }, y=resizeRanges, n=names(resizeRanges))

# Annotate Differetial NCoR1 peaks to nearest genes
peakAnnoList <- lapply(NCoR1_diffPeak.list, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000),
                       annoDb= "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))

peak_annotation = lapply(peakAnnoList, annotation)
peak_annotation = do.call(cbind,peak_annotation)[c(1,3,4,5),c(1,2,4,6,8,10)]
colnames(peak_annotation) = gsub("\\.Freq","",colnames(peak_annotation))
colnames(peak_annotation)[1] = "Features"
peak_annotation_melt = melt(peak_annotation)
ggplot(peak_annotation_melt,aes(x=variable,y=value,fill=Features))+
       geom_bar(stat = "identity")+
       coord_flip() +
       gg_theme

# Create data frame of annotated differential NCoR1 binding Cluster
diff_NCoR1.annotation.df = lapply(peakAnnoList, function(i) as.data.frame(i)) %>%
                                do.call(rbind,.) %>%
                                mutate(Clusters= gsub("\\..*","",rownames(.),""))

diff_NCoR1.annotation.df_10kb = diff_NCoR1.annotation.df %>% filter(abs(distanceToTSS) <= 10000)
diff_NCoR1.annotation.df.list = list()
for (i in unique(diff_NCoR1.annotation.df$Clusters)){
        
        cluster = diff_NCoR1.annotation.df %>% filter(Clusters ==i ) %>% dplyr::select("SYMBOL") %>% unique(.) %>% pull(.,"SYMBOL")
        diff_NCoR1.annotation.df.list[[i]]= cluster
}

Emp_DE.list = readRDS("results/DEG_Control.Stimulated_vs_Unstimulated.rds")
write_NCoR1.bed = function(x){
                        diff_NCoR1.annotation.df %>% 
                        filter(Clusters == x) %>% 
                        filter(SYMBOL %in% unlist(Emp_DE.list)) %>% 
                        dplyr::select(c(1:3)) %>%
                        write.table(.,file = paste0("results/NCoR1_ChIPSeq_peaks/",x,"_associated_with_DE_target.bed"),
                                    ,sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
}                       
write_NCoR1.bed("CpG_pIC_all3")
write_NCoR1.bed("CpG_all3") 
write_NCoR1.bed("Uns_pIC") 
write_NCoR1.bed("Uns") 
write_NCoR1.bed("No_Change") 

de_genes = read.csv("results/NCoR1_LRT_cluster_genes.csv",sep="\t",header = T)[,1]
gene.df = lapply(diff_NCoR1.annotation.df.list, function(i) i[i %in% cluster_df[which(cluster_df$cluster ==1),]$gene_symbol])
#subset gene list
        ncor_cluster_exp.list = list()
        ncor_cluster_gene_exp.list = list()
        for (j in 1:3){
                
                print(j)
                if( j  ==1){
                        common = Reduce(intersect, list(gene.df[[j]],gene.df[[j+1]],gene.df[[j+2]]))
                        diff = setdiff(gene.df[[j]],gene.df[[j+1]])
                        diff = setdiff(diff,gene.df[[j+2]])
                        com = Reduce(intersect, list(gene.df[[j]],gene.df[[j+1]]))
                        CpG_pIC_all3_vs_CpG_all3 = setdiff(com,common)
                        ncor_cluster_gene_exp.list[["CpG_pIC_all3_vs_CpG_all3"]] = CpG_pIC_all3_vs_CpG_all3
                        
                }
                print(length(diff))
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


# Pathway enrichment analysis of genes associated with differential NCoR1 binding Cluster
term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
diff_NCoR1.annotation.df.list_toupper =lapply(diff_NCoR1.annotation.df.list, toupper)
compEnricher =compareCluster(diff_NCoR1.annotation.df.list_toupper, fun="enricher",TERM2GENE=term2gene) %>% 
                        as.data.frame()

# Pathway barplot
compEnricher %>%
        as.data.frame() %>% 
        group_by(Cluster) %>% 
        slice(1:5) %>%
        ungroup %>% 
        mutate(Cluster = factor(Cluster,levels = sort(unique(Cluster))),
               Description = gsub("REACTOME|_"," ",Description)) %>% 
        #Description = reorder(Description, pvalue)) 
        as.data.frame(.) %>% 
        ggplot(aes(x=reorder(Description,-pvalue), y=-log(pvalue), fill = p.adjust)) +
        geom_bar(stat="identity") +
        #geom_text(aes(label =str_wrap(geneID, width = 130), y = 0), color = "black", hjust = 0,size=4)+
        facet_wrap(~Cluster, scales = "free",ncol = 1) +
        coord_flip() +
        
        #scale_x_reordered() +
        scale_x_discrete(labels = function(x) str_wrap(x, width = 30))+
        scale_fill_gradient(low="red",high="blue")+
        scale_y_continuous(expand = c(0,0)) +
        gg_theme

###########################################################################################
# NCoR1 peak annotation and binding enrichment at differential peak cluster
NCoR1.peak= "/home/imgsb/Gyan/NCOR1/Analysis_result_of_NCoR1_ChIPSeq/peak_file/NCoR1_all_cond_merged_peaks.bed"
NCoR1_peakAnnoList <- lapply(NCoR1.peak, annotatePeak, TxDb=txdb,tssRegion=c(-1000, 1000),
                       annoDb= "org.Mm.eg.db",
                       genomicAnnotationPriority = c("Promoter", "Exon", "Intron", "Intergenic"))
NCoR1.peak.annotation = lapply(NCoR1_peakAnnoList, function(i) as.data.frame(i)) %>%
        do.call(rbind,.)
NCoR1= read.csv("../../NCOR1/Analysis_result_of_NCoR1_ChIPSeq/peak_file/NCoR1_all_cond_merged_peaks.annotation",sep = "\t",header = T)[,c(1,16,20:23)]
colnames(NCoR1) = c("PeakID","Gene","Uns","CpG","pIC","all3")
################################################################################################
# density of Transcription factor ChIPseq data (Lajos group) at NCoR1 binding site and average at gene level. 
NCoR1_TF.density = read.csv("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/NCoR1_all_cond_merged_peaks_TFs_density_mod.txt",sep = "\t",header = T)
head(NCoR1_TF.density)
NCoR1_TF.density = NCoR1_TF.density %>% 
                        filter(PeakID %in% diff_NCoR1.annotation.df[which(diff_NCoR1.annotation.df$Clusters %in% c("CpG_pIC_all3","CpG_all3","Uns_pIC")),]$V4)

NCoR1.peak.annotation.density = merge(NCoR1_TF.density, NCoR1.peak.annotation,by.x="PeakID",by.y="V4",all.x=TRUE) %>% #head()
                                        dplyr::select(c(39,2:22)) %>% #head()
                                        dplyr::group_by(SYMBOL) %>% summarise_all(list(mean))

##########################################################################################
# Density of Transcription factor  ChIPseq data (Ido Amits's group) at NCoR1 binding site and average at gene level
NCoR1_TF_ido_amit.density = read.csv("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/NCoR1_all_cond_merged_peaks_Ido_amit_TFs_density_mod.txt",sep = "\t",header = T)
head(NCoR1_TF_ido_amit.density)
NCoR1_TF_ido_amit.density = NCoR1_TF_ido_amit.density %>% 
                                filter(PeakID %in% diff_NCoR1.annotation.df[which(diff_NCoR1.annotation.df$Clusters %in% c("CpG_pIC_all3","CpG_all3","Uns_pIC")),]$V4)

###########################################################################################
# Extract NCoR1 diff site for NCoR1_LRT_cluster
NCoR1.peak.annotation %>% 
        filter(V4 %in% NCoR1_TF.density$PeakID) %>% 
        filter(SYMBOL %in% NCoR1_bound_genes_clusters_df[NCoR1_bound_genes_clusters_df$cluster %in% c(1,7,11),]$genes) %>% 
        dplyr::select(c(1,2,3,6)) %>% 
        write.table(.,file = "NCoR1_LRT_cluster1.7.11_associated_diff_NCoR1_peaks.bed",sep="\t", row.names = FALSE,col.names = FALSE,quote = FALSE)
NCoR1.peak.annotation %>% 
        filter(V4 %in% NCoR1_TF.density$PeakID) %>% 
        filter(SYMBOL %in% NCoR1_bound_genes_clusters_df[NCoR1_bound_genes_clusters_df$cluster %in% c(2,3,6,10),]$genes) %>% 
        dplyr::select(c(1,2,3,6)) %>% 
        write.table(.,file = "NCoR1_LRT_cluster2.3.6.10_associated_diff_NCoR1_peaks.bed",sep="\t", row.names = FALSE,col.names = FALSE,quote = FALSE)
###########################################################################################
#############################################################################################################################################
#Combine TF and H3K27ac at NCoR1 sites

colnames(NCoR1_TF.density)[c(2,3,4,5)] = c("NCoR1_0hr","NCoR1_CpG","NCoR1_pIC","NCoR1_all3")
NCoR1_TF.H3K27ac.density = merge(NCoR1_TF.density,H3K27ac.density,by="PeakID")
pdf("plots/NCoR1_H3K27ac_correlation.pdf",width = 15,height = 5)
NCoR1_TF.H3K27ac.density.plot = list()
NCoR1_TF.H3K27ac.density.plot[[1]] = ggplot(NCoR1_TF.H3K27ac.density,aes(x=log2(NCoR1_CpG/NCoR1_0hr.x),y=log2(Emp_6hr_CpG/Emp_0hr)))+
        geom_point(color="grey",size=0.3)+
        geom_point(data =subset(NCoR1_TF.H3K27ac.density,log2(NCoR1_CpG/NCoR1_0hr.x) >=0.58 & log2(Emp_6hr_CpG/Emp_0hr) >=0.58),
                   aes(x=log2(NCoR1_CpG/NCoR1_0hr.x),y=log2(Emp_6hr_CpG/Emp_0hr)),color="red",size=0.3)+
        geom_point(data =subset(NCoR1_TF.H3K27ac.density ,
                                log2(NCoR1_CpG/NCoR1_0hr.x) <=-0.58 & log2(Emp_6hr_CpG/Emp_0hr) <=-0.58),
                   aes(x=log2(NCoR1_CpG/NCoR1_0hr.x),y=log2(Emp_6hr_CpG/Emp_0hr)),color="blue",size=0.3)+
        stat_cor()+
        geom_vline(xintercept = 0,linetype = "dashed",size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed",size = 0.5) +
        gg_theme+
        xlab("log2(6hr CpG/0hr)")+
        ylab("log2(6hr CpG/0hr)") 


NCoR1_TF.H3K27ac.density.plot[[2]] = ggplot(NCoR1_TF.H3K27ac.density,aes(x=log2(NCoR1_pIC/NCoR1_0hr.x),y=log2(Emp_6hr_pIC/Emp_0hr)))+
        geom_point(color="grey",size=0.3)+
        geom_point(data =subset(NCoR1_TF.H3K27ac.density,log2(NCoR1_pIC/NCoR1_0hr.x) >=0.58 & log2(Emp_6hr_pIC/Emp_0hr) >=0.58),
                   aes(x=log2(NCoR1_pIC/NCoR1_0hr.x),y=log2(Emp_6hr_pIC/Emp_0hr)),color="red",size=0.3)+
        geom_point(data =subset(NCoR1_TF.H3K27ac.density,log2(NCoR1_pIC/NCoR1_0hr.x) <=-0.58 & log2(Emp_6hr_pIC/Emp_0hr) <=-0.58),
                   aes(x=log2(NCoR1_pIC/NCoR1_0hr.x),y=log2(Emp_6hr_pIC/Emp_0hr)),color="blue",size=0.3)+
        stat_cor()+
        geom_vline(xintercept = 0,linetype = "dashed",size = 0.5) +
        geom_hline(yintercept = 0, linetype = "dashed",size = 0.5) +
        gg_theme+
        xlab("log2(6hr pIC/0hr)")+
        ylab("log2(6hr pIC/0hr)")
        

NCoR1_TF.H3K27ac.density.plot[[3]] = ggplot(NCoR1_TF.H3K27ac.density,aes(x=log2(NCoR1_all3/NCoR1_0hr.x),y=log2(Emp_6h_CpG_pIC/Emp_0hr)))+
        geom_point(color="grey",size=0.3)+
        geom_point(data =subset(NCoR1_TF.H3K27ac.density,log2(NCoR1_all3/NCoR1_0hr.x) >=0.58 & log2(Emp_6h_CpG_pIC/Emp_0hr) >=0.58),
                   aes(x=log2(NCoR1_all3/NCoR1_0hr.x),y=log2(Emp_6h_CpG_pIC/Emp_0hr)),color="red",size=0.3)+
        geom_point(data =subset(NCoR1_TF.H3K27ac.density,log2(NCoR1_all3/NCoR1_0hr.x) <=-0.58 & log2(Emp_6h_CpG_pIC/Emp_0hr) <=-0.58),
                   aes(x=log2(NCoR1_all3/NCoR1_0hr.x),y=log2(Emp_6h_CpG_pIC/Emp_0hr)),color="blue",size=0.3)+
        geom_vline(xintercept = 0,linetype = "dashed",size = 0.5) +
        stat_cor()+
        geom_hline(yintercept = 0, linetype = "dashed",size = 0.5) +
        gg_theme+
        xlab("log2(6hr CpG+pIC/0hr)")+
        ylab("log2(6hr CpG+pIC/0hr)")

png(filename = "plots/NCoR1_H3K27ac_correlation.png",width = 25,height = 8,units = "cm",res=400)
do.call(grid.arrange,c(NCoR1_TF.H3K27ac.density.plot,ncol=3))
ggsave("plots/NCoR1_H3K27ac_correlation.png",height = 5,width = 12)
dev.off()

############################################################################################################################
#Distribution of TFs (boxplot) on Diff NCoR1 bindnig sites ; (NCoR1, PU.1, RelA, RelB, cREl, IRF3)
tmp = merge(diff_NCoR1.annotation.df,NCoR1_TF.density,by.x="V4",by.y="PeakID") %>%
        dplyr::select(c(19:ncol(.)))
head(tmp)  

colnames(tmp)[c(2,3,4,5)] = c("NCoR1_0hr","NCoR1_CpG","NCoR1_pIC","NCoR1_all3")
tmp = tmp %>% melt(.) %>% separate(variable, c("Factor","Stimulation"), "_") 
head(tmp)
tmp$Clusters = factor(tmp$Clusters,levels = c("No_Change","Uns","CpG_all3","CpG_pIC_all3","Uns_pIC"))
tmp$Factor = factor(tmp$Factor, levels = c("NCoR1","PU.1","Rela","RelB","cREL","IRF3","p300","IRF5"))
tmp %>% filter(!(Stimulation == "all3")) %>%   filter(!(Factor %in% c("p300","IRF5"))) %>%  
        ggplot(.,aes(x=Clusters,y=log(value+1),fill=Stimulation))+
        #geom_violin(trim=FALSE) +
        geom_boxplot(aes(group=interaction(Factor,Stimulation,Clusters, fill="Stimulation")),width=0.5, position=position_dodge(0.5),outlier.shape=NA) +
        facet_grid(~Factor) +
        gg_theme_xrot +
        ylab("log(Normalized Tag Count)")
pdf("plots/IRF3_ChIP-density_at_NCoR1_sites.pdf",height = 5,width = 5.5)
tmp %>% filter(!(Stimulation == "all3")) %>%   filter((Factor %in% c("IRF3"))) %>%  
        
        ggplot(.,aes(x=Clusters,y=log(value+1),fill=Stimulation))+
        #geom_violin(trim=FALSE) +
        geom_boxplot(aes(group=interaction(Factor,Stimulation,Clusters, fill="Stimulation")),width=0.5, position=position_dodge(0.5),outlier.shape=NA) +
        facet_grid(~Factor) +
        gg_theme_xrot +
        ylab("log(Normalized Tag Count)")

dev.off()

################################################################################################################################
# Distribution of number of differential NCoR1 binding sites on CpG, pIC and common genes
df = diff_NCoR1.annotation.df %>% 
        dplyr::select(c(17,19)) %>% 
        filter(SYMBOL %in%  unique(unlist(CpG_pIC_specific_genes))) %>% 
        mutate(Gene_exp = case_when(SYMBOL %in% CpG_specific_genes ~ "CpG",
                                    SYMBOL %in% pIC_specific_genes ~ "pIC",
                                    SYMBOL %in% CpG_pIC_common_genes ~ "CpG_pIC")) %>% 
        filter(!is.na(Gene_exp ))
ncor1.dist = dcast(df, SYMBOL ~ Clusters, value.var = "Clusters") 

df = ncor1.dist %>% mutate(Gene_exp = case_when(SYMBOL %in% CpG_specific_genes ~ "CpG",
                                                SYMBOL %in% pIC_specific_genes ~ "pIC",
                                                SYMBOL %in% CpG_pIC_common_genes ~ "CpG_pIC"))
head(df)
ncor.dist.melt = melt(df,id.vars = c("Gene_exp","SYMBOL")) %>%  filter(value >0)

pdf("NCoR1_binding_distribution_CpG_pIC_specific_genes")
ggplot(ncor.dist.melt,aes(x=Gene_exp,y=value,color=variable))+geom_boxplot() +
        gg_theme+
        ylab("No. of NCoR1 binding per gene")+
        xlab("")+
        guides(color=guide_legend(title="NCoR1 binding \n category"))
theme_bw(30)

###########################################################################################################################

# Emperical cummulative distribution of NCoR1 binding on CpG and pIC specific and common genes
CpG_pIC_specific_genes = readRDS("results/CpG_pIC_specific_genes.rds")
pIC = merge(NCoR1_TF.density, diff_NCoR1.annotation.df,
            by.x="PeakID",by.y="V4",all.x=TRUE) %>% 
            filter(SYMBOL %in% CpG_pIC_specific_genes[[3]])%>% 
            dplyr::select(c(38,2,3,4,5)) %>% melt()
ks.test(pIC[which(pIC$variable=="pIC"),]$value,pIC[which(pIC$variable=="CpG"),]$value)
ks.test(pIC[which(pIC$variable=="pIC"),]$value,pIC[which(pIC$variable=="all3"),]$value)
ks.test(pIC[which(pIC$variable=="CpG"),]$value,pIC[which(pIC$variable=="all3"),]$value)

df_temp = ddply(pIC, .(variable), summarize,
          value = (value),
          ecdf = ecdf(value)((value)))
p1 = df_temp %>% #filter(value <=50)%>%
     ggplot(.,aes(x=log(value), y=ecdf, color = variable)) +
        geom_line(size=1)+
        #scale_colour_manual(values = c("black","blue","red","green"))+
        gg_theme+
        ylab("Cummulative Distribution")+
        xlab(" Log(Normalized Tag Count)")+
        ggtitle("pIC")

CpG_pIC = merge(NCoR1_TF.density, diff_NCoR1.annotation.df,,by.x="PeakID",by.y="V4",all.x=TRUE) %>% 
          filter(SYMBOL %in% CpG_pIC_specific_genes[[2]]) %>% 
          dplyr::select(c(38,2,3,4,5)) %>% melt()
ks.test(CpG_pIC[which(CpG_pIC$variable=="pIC"),]$value,CpG_pIC[which(CpG_pIC$variable=="CpG"),]$value)
ks.test(CpG_pIC[which(CpG_pIC$variable=="pIC"),]$value,CpG_pIC[which(CpG_pIC$variable=="all3"),]$value)
ks.test(CpG_pIC[which(CpG_pIC$variable=="CpG"),]$value,CpG_pIC[which(CpG_pIC$variable=="all3"),]$value)
df_temp = ddply(CpG_pIC, .(variable), summarize,
                value = (value),
                ecdf = ecdf(value)((value)))
p2 = df_temp %>% #filter(value <=50)%>%
      ggplot(.,aes(x=log(value), y=ecdf, color = variable)) +
      geom_line(size=1)+
      #scale_colour_manual(values = c("black","blue","red","green"))+
      gg_theme+
      ylab("Cummulative Distribution")+
      xlab(" Log(Normalized Tag Count)")+
      ggtitle("CpG and pIC")
p2

CpG = merge(NCoR1_TF.density, diff_NCoR1.annotation.df,,by.x="PeakID",by.y="V4",all.x=TRUE) %>% 
      filter(SYMBOL %in% CpG_pIC_specific_genes[[1]]) %>% 
      dplyr::select(c(38,2,3,4,5)) %>% melt()
ks.test(CpG[which(CpG$variable=="pIC"),]$value,CpG[which(CpG$variable=="CpG"),]$value)
ks.test(CpG[which(CpG$variable=="pIC"),]$value,CpG[which(CpG$variable=="all3"),]$value)
ks.test(CpG[which(CpG$variable=="CpG"),]$value,CpG[which(CpG$variable=="all3"),]$value)
df_temp = ddply(CpG, .(variable), summarize,
                value = (value),
                ecdf = ecdf(value)((value)))
p3 = df_temp %>% #filter(value <=50)%>%
     ggplot(.,aes(x=log(value), y=ecdf, color = variable)) +
        geom_line(size=1)+
        #scale_colour_manual(values = c("black","blue","red","green"))+
        gg_theme+
        ylab("Cummulative Distribution")+
        xlab(" Log(Normalized Tag Count)")+
        ggtitle("CpG")
p3
pdf("plots/NCoR1_binding_ecdf_dist.pdf",width = 14,height = 3)
grid.arrange(p3,p2,p1,ncol=3)
dev.off()

######################################################################################################################
save(diff_NCoR1.annotation.df,
     diff_NCoR1.annotation.df.list,
     NCoR1_TF.density,
     file="results/NCoR1_ChIP-Seq.RData")

