###########################################################################################################################################################################################################################################
gs.RNASeq = 25000
# Overlap of CpG and pIC specific genes with NCoR1 KD genes list
NCoR1_bound_Emp_rnaseq = newGOM(diff_NCoR1.annotation.df.list,CpG_pIC_specific_genes, gs.RNASeq)

NCoR1_bound_Emp_rnaseq.inter.nl = getNestedList(NCoR1_bound_Emp_rnaseq,"intersection")
NCoR1_bound_Emp_rnaseq.overlap_matrix = getMatrix(NCoR1_bound_Emp_rnaseq,name = "intersection") %>% 
                                                  melt(.) %>% filter(Var1  !="No_Change") 
NCoR1_bound_Emp_rnaseq.overlap_matrix$Reg = gsub(".*_vs_.*_","",NCoR1_bound_Emp_rnaseq.overlap_matrix$Var2)
NCoR1_bound_Emp_rnaseq.overlap_matrix$Var2 = gsub("_up|_down","",NCoR1_bound_Emp_rnaseq.overlap_matrix$Var2)

ggplot(NCoR1_bound_Emp_rnaseq.overlap_matrix,aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity") +
  gg_theme
###########################################################################################################################################################################################################################################
# Overlap of CpG and pIC specific genes with NCoR1 KD genes list

NCoR1_Emp_rnaseq_cluster = newGOM(CpG_pIC_specific_genes,c(DE.list,NCoR1_KD_No.change.list), gs.RNASeq)

NCoR1_Emp_rnaseq_cluster.inter.nl = getNestedList(NCoR1_Emp_rnaseq_cluster,"intersection")
NCoR1_Emp_rnaseq_cluster.overlap_matrix = getMatrix(NCoR1_Emp_rnaseq_cluster,name = "intersection") %>% 
  melt(.) %>% filter(Var1  !="No_Change") 
NCoR1_Emp_rnaseq_cluster.overlap_matrix$Reg = gsub(".*_vs_.*_","",NCoR1_Emp_rnaseq_cluster.overlap_matrix$Var2)
NCoR1_Emp_rnaseq_cluster.overlap_matrix$Var2 = gsub("_up|_down|_noChange","",NCoR1_Emp_rnaseq_cluster.overlap_matrix$Var2)
NCoR1_Emp_rnaseq_cluster.overlap_matrix$Reg = factor(NCoR1_Emp_rnaseq_cluster.overlap_matrix$Reg,
                                                     levels = c("up","noChange","down"),
                                                     labels = c("Up-regulated","No Change","Down-regulated"))
NCoR1_Emp_rnaseq_cluster.overlap_matrix$Var2 = factor(NCoR1_Emp_rnaseq_cluster.overlap_matrix$Var2,
                                                      levels = c("NC_vs_EC","NP_vs_EP","NA_vs_EA"),
                                                      labels = c("NCoR1 CpG","NCoR1 pIC","NCoR1 CpG+pIC"))
ggplot(NCoR1_Emp_rnaseq_cluster.overlap_matrix,aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity") +
  facet_wrap(~Reg,scales = "free_x")+
  scale_y_continuous(limits = c(0, 950), breaks = seq(0, 950, by = 300))+
  scale_fill_manual(values = c("blue","darkgreen","red"))+
  gg_theme_xrot
ggsave("plots/NCoR1_KD_CpG_pIC_specific_gene_distribution.pdf",height = 5,width = 8,dpi = 300)
###########################################################################################################################################################################################################################################

CpG_pIC_specific_NCoR1_KD_genes_no. = data.frame(Var1=c("Others","Others","Others","Others","Others","Others"),
                                                 Var2=c("NCoR1 CpG","NCoR1 CpG","NCoR1 pIC","NCoR1 pIC","NCoR1 CpG+pIC","NCoR1 CpG+pIC"),
                                                 value= c(length(setdiff(DE.list[[1]],
                                                            c(NCoR1_Emp_rnaseq_cluster.inter.nl$NC_vs_EC_up$CpG,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NC_vs_EC_up$pIC,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NC_vs_EC_up$CpG_pIC))),
                                                          length(setdiff(DE.list[[2]],
                                                            c(NCoR1_Emp_rnaseq_cluster.inter.nl$NC_vs_EC_down$CpG,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NC_vs_EC_down$pIC,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NC_vs_EC_down$CpG_pIC))),
                                                          length(setdiff(DE.list[[3]],
                                                            c(NCoR1_Emp_rnaseq_cluster.inter.nl$NP_vs_EP_up$CpG,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NP_vs_EP_up$pIC,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NP_vs_EP_up$CpG_pIC))),
                                                          length(setdiff(DE.list[[4]],
                                                            c(NCoR1_Emp_rnaseq_cluster.inter.nl$NP_vs_EP_down$CpG,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NP_vs_EP_down$pIC,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NP_vs_EP_down$CpG_pIC))),
                                                          length(setdiff(DE.list[[5]],
                                                            c(NCoR1_Emp_rnaseq_cluster.inter.nl$NA_vs_EA_up$CpG,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NA_vs_EA_up$pIC,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NA_vs_EA_up$CpG_pIC))),
                                                          length(setdiff(DE.list[[6]],
                                                            c(NCoR1_Emp_rnaseq_cluster.inter.nl$NA_vs_EA_down$CpG,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NA_vs_EA_down$pIC,
                                                              NCoR1_Emp_rnaseq_cluster.inter.nl$NA_vs_EA_down$CpG_pIC)))),
                                                 Reg = c("Up-regulated","Down-regulated","Up-regulated","Down-regulated","Up-regulated","Down-regulated"))

###########################################################################################################################################################################################################################################
rbind(NCoR1_Emp_rnaseq_cluster.overlap_matrix,CpG_pIC_specific_NCoR1_KD_genes_no.)
ggplot(rbind(NCoR1_Emp_rnaseq_cluster.overlap_matrix,CpG_pIC_specific_NCoR1_KD_genes_no.),aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity") +
  facet_wrap(~Reg,scales = "free_x")+
  scale_y_continuous(limits = c(0, 1500), breaks = seq(0, 1400, by = 200))+
  scale_fill_manual(values = c("blue","darkgreen","red","grey"))+
  gg_theme
###########################################################################################################################################################################################################################################
# Overlap of NCoR1 bound genes with NCoR1 KD genes list
NCoR1_bound_NCoR1_KD_genes = newGOM(diff_NCoR1.annotation.df.list,DE.list,gs.RNASeq)
NCoR1_bound_NCoR1_KD_genes.inter.nl = getNestedList(NCoR1_bound_NCoR1_KD_genes,"intersection")
NCoR1_bound_NCoR1_KD_genes.overlap_matrix = getMatrix(NCoR1_bound_NCoR1_KD_genes,name = "intersection") %>% 
                                               melt(.) %>% filter(Var1  !="No_Change") 
NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Reg = gsub(".*_vs_.*_","",NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Var2)
NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Var2 = gsub("_up|_down","",NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Var2)
NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Reg = factor(NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Reg,
                                                      levels = c("up","down"),
                                                      labels = c("Up-regulated","Down-regulated"))
NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Var2 = factor(NCoR1_bound_NCoR1_KD_genes.overlap_matrix$Var2,
                                                      levels = c("NC_vs_EC","NP_vs_EP","NA_vs_EA"),
                                                      labels = c("NCoR1 CpG","NCoR1 pIC","NCoR1 CpG+pIC"))
ggplot(NCoR1_bound_NCoR1_KD_genes.overlap_matrix,aes(x=Var2,y=value,fill=Var1))+
        geom_bar(stat="identity") +
        scale_y_continuous(limits = c(0, 1400), breaks = seq(0, 1400, by = 200))+
        facet_wrap(~Reg,scales = "fixed")+
        #scale_fill_manual(values = c("blue","darkgreen","red"))+
        gg_theme
###########################################################################################################################################################################################################################################

# Overlap of NCoR1 bound genes with NCoR1 KD genes list
tmp.list= list(CpG_all3=NCoR1_bound_Emp_rnaseq.inter.nl$CpG_all3$CpG_pIC,
               CpG_pIC_all3=NCoR1_bound_Emp_rnaseq.inter.nl$CpG_pIC_all3$CpG_pIC,
               pIC=NCoR1_bound_Emp_rnaseq.inter.nl$Uns_pIC$CpG_pIC)

tmp.list= list(CpG_all3=NCoR1_bound_Emp_rnaseq.inter.nl$CpG_all3$CpG,
               CpG_pIC_all3=NCoR1_bound_Emp_rnaseq.inter.nl$CpG_pIC_all3$CpG,
               pIC=NCoR1_bound_Emp_rnaseq.inter.nl$Uns_pIC$CpG)

tmp.list= list(CpG_all3=NCoR1_bound_Emp_rnaseq.inter.nl$CpG_all3$pIC,
               CpG_pIC_all3=NCoR1_bound_Emp_rnaseq.inter.nl$CpG_pIC_all3$pIC,
               pIC=NCoR1_bound_Emp_rnaseq.inter.nl$Uns_pIC$pIC)

CpG_pIC_specific_NCoR1_KD_genes = newGOM(tmp.list,DE.list,gs.RNASeq)
CpG_pIC_specific_NCoR1_KD_genes.inter.nl = getNestedList(CpG_pIC_specific_NCoR1_KD_genes,"intersection")
CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix = getMatrix(CpG_pIC_specific_NCoR1_KD_genes,name = "intersection") %>% 
                                              melt(.) %>% filter(Var1  !="No_Change") 
CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Reg = gsub(".*_vs_.*_","",CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Var2)
CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Var2 = gsub("_up|_down","",CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Var2)
CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Reg = factor(CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Reg,
                                                       levels = c("up","down"),
                                                       labels = c("Up-regulated","Down-regulated"))
CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Var2 = factor(CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix$Var2,
                                                        levels = c("NC_vs_EC","NP_vs_EP","NA_vs_EA"),
                                                        labels = c("NCoR1 CpG","NCoR1 pIC","NCoR1 CpG+pIC"))
ggplot(CpG_pIC_specific_NCoR1_KD_genes.overlap_matrix,aes(x=Var2,y=value,fill=Var1))+
  geom_bar(stat="identity") +
  #scale_y_continuous(limits = c(0, 400), breaks = seq(0, 400, by = 200))+
  facet_wrap(~Reg,scales = "fixed")+
  #scale_fill_manual(values = c("blue","darkgreen","red"))+
  gg_theme_xrot


NCoR1_H3K27ac_boxplot(unique(c(NCoR1_bound_NCoR1_KD_genes.inter.nl$NC_vs_EC_up$CpG_all3,
                               NCoR1_bound_NCoR1_KD_genes.inter.nl$NC_vs_EC_up$CpG_pIC_all3,
                               NCoR1_bound_NCoR1_KD_genes.inter.nl$NC_vs_EC_up$Uns_pIC)))

NCoR1_H3K27ac_boxplot(NCoR1_bound_NCoR1_KD_genes.inter.nl$NC_vs_EC_down$CpG_all3)

gene.df = inter.nl$NP_vs_EP_down[c(2,3,4)]
gene.df = inter.nl$CpG[c(2,3,4)]
gene.df = inter.nl$CpG_pIC[c(2,3,4)]



NCoR1_H3K27ac_boxplot = function(gene_list,nam){

  #NCoR1 and H3K27ac
  tmp = merge(diff_NCoR1.annotation.df[which(diff_NCoR1.annotation.df$Clusters %in% c("CpG_all3","Uns_pIC","CpG_pIC_all3")),],
            NCoR1_TF.H3K27ac.density,by.x="V4",by.y="PeakID") %>% #head()
            filter(SYMBOL %in% gene_list) %>% 
            dplyr::select(c(17,19,20:23,42:ncol(.))) %>% melt(.)
  
  tmp = tmp %>% mutate(Data.type = case_when(tmp$variable %in% colnames(NCoR1_TF.density) ~ "NCoR1",
                                           tmp$variable %in% colnames(H3K27ac.density) ~ "H3K27ac",
                                           tmp$variable == "NCoR1_0hr.x" ~ "NCoR1",
                                           tmp$variable == "NCoR1_0hr.y" ~ "H3K27ac"))
  tmp$variable = gsub("\\.x|\\.y","",tmp$variable)
  #tmp$Factor= gsub("_.*","",tmp$variable)
  #tmp$Factor= factor(tmp$Factor,levels = c("NCoR1","Rela","cRel","IRF3","IRF5"))
  #tmp$Clusters = factor(tmp$Clusters,levels = c("Uns","CpG_all3","Uns_pIC","CpG_pIC_all3"))
  tmp$Condition = gsub("_.*","",tmp$variable)
  tmp$Stimulation = gsub("NCoR1_|Emp_","",tmp$variable)
  #tmp =tmp %>% mutate(Stimulation= case_when(tmp$Stimulation =="CpG_pIC" ~ "all3"))
  #tmp$Stimulation =factor(tmp$Stimulation,levels = c("0hr","CpG","pIC","all3","2hr_CpG","6hr_CpG","2hr_pIC","6hr_pIC","2h_CpG_pIC","6h_CpG_pIC"),
  #                        labels = c("0hr","6hr_CpG","6hr_pIC","6hr_CpG_pIC","2hr_CpG","6hr_CpG","2hr_pIC","6hr_pIC","2h_CpG_pIC","6h_CpG_pIC"))

  #tmp$Data.type = factor(tmp$Data.type,levels = c("NCoR1","H3K27ac"))
  tmp= tmp[,c(3:ncol(tmp))]
  
  
  
  # Expression
  Exp = Emp_NCoR1_rep_merged[gene_list,] %>% melt(.)
  Exp$Data.type= "Expression"
  Exp$Condition= gsub("_.*","",Exp$variable)
  Exp$Stimulation = gsub("Emp_|NCoR1_","",Exp$variable)
  Exp$Stimulation = factor(Exp$Stimulation,levels = c("Uns","IFNg","6hr_CpG","6hr_pIC","6hr_all3"),labels = c("0hr","IFNg","6hr CpG","6hr pIC","6hr CpG+pIC"))
    
  Exp_tmp = rbind(Exp,tmp)
  Exp_tmp$Data.type = factor(Exp_tmp$Data.type,levels = c("NCoR1","H3K27ac","Expression"))

  p1 = Exp_tmp %>% filter(Data.type=="NCoR1") %>%
                   mutate(Stimulation =factor(Stimulation,levels = c("Uns","CpG","pIC","all3"),
                                                          labels = c("0hr","6hr_CpG","6hr_pIC","6hr_CpG_pIC"))) %>%
                   ggplot(.,aes(x=Stimulation,y=log(value+0.1),color=Stimulation))+
                      geom_boxplot(lwd=0.5)+
                      facet_wrap(~Data.type,scales = "free")+
                      gg_theme_xrot+
                      ggtitle(nam)+
                      theme(legend.position = "none")+
                      stat_compare_means(comparisons = list(c("6hr_CpG","6hr_pIC"),c("6hr_pIC","6hr_CpG_pIC"),c("6hr_CpG","6hr_CpG_pIC")), 
                      label = "p.signif",  label.x = 1.5,method = "wilcox",paired = FALSE)
  
  p2 = Exp_tmp %>% filter(Data.type=="H3K27ac") %>%
          mutate(Stimulation = factor(Stimulation,levels = c("0hr","2hr_CpG","6hr_CpG","2hr_pIC","6hr_pIC","2h_CpG_pIC","6h_CpG_pIC"),
                                                  labels = c("0hr","2hr_CpG","6hr_CpG","2hr_pIC","6hr_pIC","2h_CpG_pIC","6h_CpG_pIC"))) %>%
          ggplot(.,aes(x=Stimulation,y=log(value+0.1),color=Condition))+
          geom_boxplot(lwd=0.5)+
          facet_wrap(~Data.type,scales = "free")+
          gg_theme_xrot+
          ggtitle(nam)+
          theme(legend.position = "none")+
          stat_compare_means(aes(group = Condition), label = "p.signif",  label.x = 1.5,method = "wilcox",paired = FALSE)
  
  p3 = Exp_tmp %>% filter(Data.type=="Expression") %>%
          ggplot(.,aes(x=Stimulation,y=value,color=Condition))+
          geom_boxplot(lwd=0.5)+
          facet_wrap(~Data.type,scales = "free")+
          gg_theme_xrot+
          ggtitle(nam)+
          theme(legend.position = "none")+
          stat_compare_means(aes(group = Condition), label = "p.signif",  label.x = 1.5,method = "wilcox",paired = FALSE) 
  
  return(do.call(grid.arrange,c(list(p1,p2,p3),ncol=3)))
  
}
pdf("/home/imgsb/Gyan/tmp/test2.pdf",width = 10,height = 5)
NCoR1_H3K27ac_boxplot(unlist(NCoR1_bound_genes_clusters_df.list[1]),"NcoR1 KD CpG")
NCoR1_H3K27ac_boxplot(unlist(NCoR1_bound_genes_clusters_df.list[c(3,4,5,6)]),"NCoR1 KD pIC")
#lapply(seq_along(NCoR1_bound_genes_clusters_df.list), function(y, n, i) { NCoR1_H3K27ac_boxplot(y[[i]],n[[i]]) }, y=NCoR1_bound_genes_clusters_df.list, n=names(NCoR1_bound_genes_clusters_df.list))
dev.off()
##################################################################################################################################################
# Transcription factor binding at NCoR1 binding sites associated with CpG, pIC and CpG & pIC response genes

TF_binding = function(gene_list,title){
  tmp = merge(NCoR1_TF.density, NCoR1.peak.annotation,by.x="PeakID",by.y="V4",all.x=TRUE) %>% 
            dplyr::select(c(39,9:11,15:17,18:23)) %>% filter(SYMBOL %in% gene_list) %>% melt(.)
  tmp$Factor = gsub("_.*","",tmp$variable)
  tmp$Condition = gsub(".*_","",tmp$variable)
  tmp$Factor = factor(tmp$Factor,levels = c("Rela","RelB","cREL","IRF3","IRF5"),
                      labels = c("RelA","RelB","cREL","IRF3","IRF5"))
  tmp %>% filter(!(Factor %in%  c("IRF5","RelB","RelA")))%>%
          filter(value !=0) %>% 
            ggplot(.,aes(x=Condition,y=log(value+1))) + 
            geom_violin(aes(color = Condition), position = position_dodge(0.8),size=1) +
            geom_boxplot(aes(group=interaction(Factor,Condition)),width=0.1, fill="white", position=position_dodge(0.8),outlier.shape=NA) +
            gg_theme_xrot +
            facet_grid(~Factor,scales = "free_x",space= "free_x") +
            #theme(axis.text.x = element_blank(),
             #     axis.ticks = element_blank()) +
            stat_compare_means(aes(group = Condition), comparisons = list(c("0hr","CpG"),c("0hr","pIC"),c("CpG","pIC")),
                               label = "p.signif",method = "wilcox",paired = FALSE) +
            ylab("Log(normalized count)")+
            ylim(0,11)+
            ggtitle(label = title)
}

pdf("plots/TF_binding_CpG_pIC_specific_genes.pdf",width = 6,height = 4)
TF_binding(unlist(NCoR1_bound_genes_clusters_df.list[c(1,2,7)]),"CpG , CpG and pIC")
TF_binding(unlist(NCoR1_bound_genes_clusters_df.list[c(3,4,5,6)]),"pIC")
TF_binding(CpG_specific_genes,"CpG specific genes")
TF_binding(pIC_specific_genes,"pIC specific genes")
TF_binding(CpG_pIC_common_genes,"CpG & pIC common genes")
dev.off()


merge(NCoR1_TF.density, NCoR1.peak.annotation,by.x="PeakID",by.y="V4",all.x=TRUE) %>% 
        dplyr::select(c(39,9:11,15:17,18:23)) %>% 
        filter(SYMBOL %in% unlist(CpG_pIC_specific_genes)) %>% 
  melt(.) %>% filter(grepl("IRF3_pIC",variable)) %>% 
  mutate(genes = case_when(SYMBOL %in% CpG_specific_genes ~ "CpG", 
                           SYMBOL %in% pIC_specific_genes ~ "pIC" , 
                           SYMBOL %in% CpG_pIC_common_genes ~ "CpG_pIC")) %>% 
  ggplot(aes(x=genes,y=log(value+1),color=genes)) + 
  geom_violin() +
  gg_theme +
  geom_boxplot(width=0.1, fill="white", position=position_dodge(0.8),outlier.shape=NA) +
  stat_compare_means(aes(group = genes), comparisons = list(c("CpG","pIC"),c("CpG","CpG_pIC"),c("pIC","CpG_pIC")),
                     label = "p.signif",method = "wilcox",paired = FALSE) +
  stat_summary(fun=median, geom="line", aes(group=1))
  
ggsave("plots/IRF3_binding_at_CpG_pIC_specific_genes.pdf",height = 5,width = 6)
###################################################################################################
#NCoR1 overlap with TF binding
NCor1_h3k27ac_overlap = read.csv("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/SRTFs_overlap_with_NCoR1/Number_of_H3K27ac_at_NCoR1_associated_with_DE_target.txt",
                                 sep="\t",header = FALSE)
SRTFs_overlap = read.csv("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/SRTFs_overlap_with_NCoR1/Number_of_SRTF_overlap.txt",
                         sep="\t",header = FALSE)
NcoR1_H3k27ac_SRTF_overlap = merge(SRTFs_overlap,NCor1_h3k27ac_overlap,by="V2")
#NcoR1_H3k27ac_SRTF_overlap =rbind(NCor1_h3k27ac_overlap,SRTFs_overlap)
NcoR1_H3k27ac_SRTF_overlap = NcoR1_H3k27ac_SRTF_overlap %>% dplyr::group_by(V2,V3.x) %>% mutate(percent = V1.x/V1.y*100)
NcoR1_H3k27ac_SRTF_overlap$stimulation = gsub("_Rep[0-9]|[0-9]*min_","",NcoR1_H3k27ac_SRTF_overlap$V4.x)
NcoR1_H3k27ac_SRTF_overlap$V2 = factor(NcoR1_H3k27ac_SRTF_overlap$V2,levels = c("No_change","Uns","CpG_pIC_all3","CpG_all3","Uns_pIC"))
NcoR1_H3k27ac_SRTF_overlap$V3.x = factor(NcoR1_H3k27ac_SRTF_overlap$V3.x,levels = c("PU1","JUNB","cREL","IRF3","IRF5","H4ac","p300","STAT1"))
TF_percent_plt =  NcoR1_H3k27ac_SRTF_overlap %>% 
  filter(!(V3.x %in% c("H4ac","p300","STAT1"))) %>% 
  filter(!(V4.x %in% c("90min_CpG","102min_pIC_Rep2"))) %>%  
  filter(V3.x %in% c("PU1","IRF3","JUNB","cREL")) %>% 
  #filter(V2 %in% c("CpG_pIC_all3","CpG_all3","Uns_pIC")) %>%
  ggplot(.,aes(x=stimulation,y=percent,fill=V3.x))+#geom_tile(colour = "black")+gg_theme_xrot+
  #scale_fill_gradient2(low = "white", high = "red") 
  geom_bar(stat="identity",size=1)+ 
  facet_grid(V2~V3.x,scales = "free_x")+gg_theme_xrot+
  ylab("Percent overlap with NCoR1 binding sites")+
  theme(strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18))
pdf("plots/TF_percent_overlap_with_NCoR1.pdf",height = 10,width = 6)
TF_percent_plt
dev.off()
#######################################################################################################

# Correlate gene expression and H3K27ac in NCoR1 KD samples

gs.RNASeq = 25000 # Total no. of genes associated with diff H3K27ac peaks
control_KD_overlap = newGOM(Emp_DE.list[c(5,6)], Emp_DE.list[c(3,4)],
                            
                            genome.size=gs.RNASeq)


control_KD_overlap_intersection = getMatrix(control_KD_overlap,name="intersection") 
-log10(.) %>%
  sapply(.,round,2) %>% 
  as.data.frame() %>%
  sapply(.,as.character)
control_KD_overlap.inter.nl = getNestedList(control_KD_overlap,"intersection")


Exp_H3K27ac.plot = function(exp,h3k27ac,tittle){
  tmp = exp %>% 
    filter(padj <=0.05 & abs(log2FoldChange) >=1 ) %>% #
    dplyr::select(1,3)
  tmp1 = h3k27ac %>% 
    #filter(padj <=0.05 & abs(log2FoldChange) >=0.58) %>% #
    merge(.,diff_h3K27ac_2k_summit_peak.bed.annotation.df,by.x = 0,by.y="PeakID") %>% 
    dplyr::select(c(1,3,23))
  tmp_merged = merge(tmp1,tmp,by.x="SYMBOL",by.y="Gene") %>% # View()
    dplyr::select(1,3,4) %>% 
    group_by(SYMBOL) %>% 
    dplyr::summarise_at(c("log2FoldChange.x", "log2FoldChange.y"), mean, na.rm = TRUE) %>%
    mutate(genes= case_when
           (SYMBOL %in% unlist(NCoR1_bound_genes_clusters_df.list[c(3,4,5,6)],use.names = F) ~ "pIC", 
             SYMBOL %in% unlist(NCoR1_bound_genes_clusters_df.list[1],use.names = F) ~ "CpG",
             SYMBOL %in% unlist(NCoR1_bound_genes_clusters_df.list[c(2,7)],use.names = F) ~ "CpG & pIC"
           )) %>% #View()
    filter(SYMBOL %in% unique(unlist(diff_NCoR1.annotation.df.list,use.names=FALSE))) #%>%
  data = subset(tmp_merged,SYMBOL %in% labRow & genes != "Na" & log2FoldChange.x < 0)
  ggplot(tmp_merged,aes(x=log2FoldChange.x,y=log2FoldChange.y)) +  
    geom_point(data=tmp_merged[is.na(tmp_merged$genes),], aes(x=log2FoldChange.x,y=log2FoldChange.y),alpha=1,color="#424040",size=0.6) +
    geom_point(data=tmp_merged %>% filter(!(is.na(tmp_merged$genes)) & log2FoldChange.x < 0), aes(x=log2FoldChange.x,y=log2FoldChange.y,color=genes),size=0.6)+
    stat_cor()+
    gg_theme+
    scale_color_manual(labels = c("CpG", "CpG & pIC", "pIC","Others"),values = c("blue","darkgreen","red","darkgrey"))+
    geom_hline(yintercept=c(1,-1), linetype="dashed", color = "black")+
    geom_vline(xintercept=c(0.58,-0.58), linetype="dashed", color = "black")+
    xlab("Log2 FC\nH3K27ac")+
    ylab("Gene Expression\nLog2 FC")+
    xlim(-3.2,3.2)+ylim(-5.5,10)+
    geom_text_repel(data =data,
                     aes(x=log2FoldChange.x,y=log2FoldChange.y,
                        label = SYMBOL,
                        #fill=factor(genes),
                        segment.colour = factor(genes)),
                        direction="y",
                        nudge_x = -0.5+data$log2FoldChange.x ,
                        nudge_y = 1+data$log2FoldChange.y,
                        color = "black" ,fontface="italic",box.padding = 0.5,max.overlaps = 1000,
                        size = 4, segment.size = 0.2)+
    scale_fill_manual(values = c("blue","darkgreen","red"),aesthetics = c("fill", "segment.color"))+
    
    ggtitle(tittle)
}

Exp_H3K27ac.plot.list = list()
Exp_H3K27ac.plot.list[[1]] = Exp_H3K27ac.plot(NC_vs_EC,NCoR1.6hrCpG_vs_Emp.6hrCpG, 
                                tittle = "NCoR1 vs Emp (6hr CpG)")

Exp_H3K27ac.plot.list[[2]] = Exp_H3K27ac.plot(NP2_vs_EP2,NCoR1.2hrpIC_vs_Emp.2hrpIC,
                                tittle = "NCoR1 vs Emp (2hr pIC)")

Exp_H3K27ac.plot.list[[3]] = Exp_H3K27ac.plot(NP_vs_EP,NCoR1.6hrpIC_vs_Emp.6hrpIC,
                                tittle = "NCoR1 vs Emp (6hr pIC)")


Exp_H3K27ac.plot.list[[4]] = Exp_H3K27ac.plot(NA_vs_EA,NCoR1.6hrCpG.pIC_vs_Emp.6hrCpG.pIC,
                                tittle = "NCoR1 vs Emp (6hr CpG+pIC)")

pdf("plots/NCoR1_vs_Emp_RNASeq_H3K27ac_Log2FC.pdf",height = 6.4,width = 14)
do.call(grid.arrange,c(Exp_H3K27ac.plot.list[c(1,2)],ncol=2))
dev.off()
