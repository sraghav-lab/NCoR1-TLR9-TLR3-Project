# Integrate NCoR1 binding and NCoR1 KD RNASeq data

diff_NCoR1.annotation.df.binding = diff_NCoR1.annotation.df %>% 
                                        filter(Clusters %in% c("CpG_pIC_all3","CpG_all3","Uns_pIC"))


NCoR1_bound_genes_clusters_df = NCoR1_LRT_cluster$df %>% 
                                        filter(cluster %in% c(1,2,3,6,7,10,11)) %>% 
                                        filter(genes %in% diff_NCoR1.annotation.df.binding$SYMBOL &
                                               genes %in% unlist(CpG_pIC_specific_genes))

NCoR1_bound_genes_clusters_df = NCoR1_bound_genes_clusters_df %>% 
                                  mutate(gene_set = case_when(genes %in% CpG_pIC_specific_genes[[1]] ~ "CpG",
                                                              genes %in% CpG_pIC_specific_genes[[2]] ~ "CpG_pIC",
                                                              genes %in% CpG_pIC_specific_genes[[3]] ~ "pIC",))

NCoR1_bound_genes_clusters_df$cluster = paste0("C",NCoR1_bound_genes_clusters_df$cluster)
NCoR1_bound_genes_clusters_df = NCoR1_bound_genes_clusters_df %>% arrange(cluster,gene_set)
NCoR1_bound_genes_clusters_df$gene_cluster = paste0(NCoR1_bound_genes_clusters_df$cluster,"_",NCoR1_bound_genes_clusters_df$gene_set)
NCoR1_bound_genes_clusters_df = NCoR1_bound_genes_clusters_df %>% arrange(gene_cluster)
NCoR1_bound_genes_clusters_df[which(NCoR1_bound_genes_clusters_df$gene_cluster == "C7_pIC"),][,c(2,4)] = c("C2","C2_pIC")
NCoR1_bound_genes_clusters_df[which(NCoR1_bound_genes_clusters_df$gene_cluster == "C3_CpG"),][,c(2,4)] = c(rep("C11",5),rep("C11_CpG",5))
NCoR1_bound_genes_clusters_df[which(NCoR1_bound_genes_clusters_df$gene_cluster == "C7_CpG"),][,c(2,4)] = c(rep("C1",34),rep("C1_CpG",34))

NCoR1_bound_genes_clusters_df = NCoR1_bound_genes_clusters_df %>% 
					filter(gene_cluster %in% c("C1_CpG","C1_CpG_pIC","C7_CpG_pIC","C2_pIC","C1_pIC","C6_pIC","C3_pIC"))
NCoR1_bound_genes_clusters_df.list = split(NCoR1_bound_genes_clusters_df$genes, NCoR1_bound_genes_clusters_df$gene_cluster)

for (i  in c(1,2,3,6,7,10,11)) {
  tmp = NCoR1_bound_genes_clusters_df %>%
          filter(cluster==i) %>%
          dplyr::select(1)
    write.table(tmp,paste0("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/MARGE_analysis/","NCoR1_bound_LRT_cluster",i,".txt"),quote = FALSE,row.names = FALSE,col.names = FALSE)
  
}

for (i  in c(1,2,3,6,7,10,11)) {
  merge(NCoR1_TF.density, NCoR1.peak.annotation,,by.x="PeakID",by.y="V4",all.x=TRUE) %>% 
        filter(SYMBOL %in% NCoR1_bound_genes_clusters_df[which(NCoR1_bound_genes_clusters_df$cluster == i),]$genes) %>% 
        dplyr::select(c(23,24,25,3)) %>%
  write.table(.,paste0("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/MARGE_analysis/","NCoR1_bound_LRT_cluster",i,".bed"),quote = FALSE,row.names = FALSE,col.names = FALSE,sep="\t")
  
}
                

NCoR1_bound_genes_clusters_df.exp = Emp_NCoR1_vst_rep_merged.scaled[NCoR1_bound_genes_clusters_df$genes,]

split = factor(NCoR1_bound_genes_clusters_df$gene_cluster,levels = c("C1_CpG","C1_CpG_pIC","C7_CpG_pIC","C2_pIC","C1_pIC","C6_pIC","C3_pIC"))

RowAnn <- HeatmapAnnotation(df=NCoR1_bound_genes_clusters_df$gene_cluster, which="row",col=list(df=c("C1_CpG" = "blue",
                                                                                                     "C1_CpG_pIC"  = "darkgreen",
                                                                                                     "C1_pIC" ="red",
                                                                                                     "C2_pIC"="red",
                                                                                                     "C3_pIC"="red",
                                                                                                     "C7_CpG_pIC"= "darkgreen",
                                                                                                     "C6_pIC"="red",
                                                                                                     "C3_CpG_pIC"="darkgreen")))

p2 <- Heatmap(NCoR1_bound_genes_clusters_df.exp[,c(1,6,2,7,3,8,4,9,5,10)],km = 1,
              name="vst gene expression (scaled)",
              col= colorRamp2(c(-1.5,0,1.5),c("blue","white","red")),
              heatmap_legend_param=list(at=c(-1.5,0,1.5),color_bar="continuous", 
                                        legend_direction="horizontal", legend_width=unit(5,"cm"),
                                        title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
              
              #Split heatmap rows by gene family
              split = split,
              border=TRUE,
              # box colour
              #rect_gp = gpar(col = "black"),
              
              #Row annotation configurations
              cluster_rows=FALSE,
              cluster_row_slices = FALSE,
              show_row_dend=TRUE,
              #row_title="Transcript", #overridden by 'split' it seems
              row_title_side="left",
              row_title_gp=gpar(fontsize=8),
              show_row_names=FALSE,
              row_names_side="right",
              #row_title_rot=30,
              
              #Column annotation configuratiions
              cluster_columns=FALSE,
              column_names_side = "top",
              column_names_rot = 45,
              show_column_dend=TRUE,
              column_title="RNASeq",
              column_title_side="top",
              #column_title_gp=gpar(fontsize=10, fontface="bold"),
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

labRow <- c("Il12b","Il12a","Il6","Ido2","Socs3","Stat1","Cd274","Ccr5","Myd88","Irf7","Ifit1","Ifit3","Cd80","Tnf",
            "Ifitm1","Ifit3","Il1a","Nfkbie","Oas1b","Zeb2","Znfx1","Oasl1","Oas2","Oas1g","Mx2","Ifit2","Ctla4",
            "Tlr9","Ebi3","Runx2","Nr4a3","Bst2","Socs1","Isg15",
            "Gbp3","Gbp7","Cxcl9","Cxcr5","Ddx58","Oas3","Traf2","Zbp1","Atf3","Ccl4","Cd40","Cd69","Ido2","Ifnlr1",
	    "Nr4a1","Tnfaip3","Jak2","Cxcl10","Il27","Traf1","Tnfip3",
            "Pparg","Lag3","Stat5a","Akt3","Il10","Tgfbr1","Socs2","Ccr7","Cd83","Ifnb1","Il15","Il1b","Irf4","Stat3")

subset = match(labRow,rownames(NCoR1_bound_genes_clusters_df.exp))
label =rownames(NCoR1_bound_genes_clusters_df.exp)[c(subset)]
ra <- anno_mark(at = subset, label,which = "row",labels_gp = gpar(fontsize=5))
pdf("plottest.pdf",height = 8,width = 7)
draw(p2+rowAnnotation(mark=ra)+RowAnn)
dev.off()
#################################################################
        
