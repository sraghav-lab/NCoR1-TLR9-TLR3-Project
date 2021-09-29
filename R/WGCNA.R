library(flashClust)
enableresults/WGCNAThreads()
options(stringsAsFactors = FALSE)

dir.create("results/WGCNA")
########################################################################################################
#Data loading and cleaning
########################################################################################################
load("results/Control_RNASeq.RData")
load("results/NCoR1.KD_RNASeq.RData")

Emp_NCoR1_vst = readRDS("results/Emp_NCoR1_vst.rds")

NCoR1_LRT_cluster_df = readRDS("results/NCoR1_LRT_cluster_genes.rds")

exp_matrix <- Emp_NCoR1_vst %>% 
                  as.data.frame() %>% 
                  filter(rownames(.) %in% 
                           NCoR1_LRT_cluster_df$df$genes)
dim(exp_matrix)
names(exp_matrix)
femData <- exp_matrix
names(femData)
datExpr = as.data.frame(t(femData[,]))

dim(datExpr)
names(datExpr) = rownames(femData)


# Trait data
condition <- gsub("_R[0-9]","",rownames(datExpr))
coldata <- data.frame(row.names=rownames(datExpr),group=condition)
coldata$group <- as.numeric(factor(coldata$group))
datTraits <- coldata
########################################################################################################
# Soft power threshold calculation
########################################################################################################
powers = c(c(1:10), seq(from = 12, to=20, by=2))
#powers = c(1:20)  # in practice this should include powers up to 20.

sft=pickSoftThreshold(datExpr,dataIsExpr = TRUE,
                      powerVector = powers,corFnc = cor,
                      corOptions = list(use = 'p'),networkType = "signed hybrid")
sizeGrWindow(9, 7)
par(mfrow = c(1,2))
cex1 = 0.9


plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit, signed R^2",type="n", main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=cex1,col="red")

# Red line corresponds to using an R^2 cut-off
abline(h=0.7,col="red")

# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
##########################################################################################################




################################################################################################################
#Step-by-step network construction and module detection, including scaling of Topological Overlap Matrices
################################################################################################################
softPower =12
TOM=TOMsimilarityFromExpr(datExpr,
                          networkType = "signed hybrid", 
                          TOMType = "signed", 
                          power = softPower)


colnames(TOM) =rownames(TOM) = datExpr
dissTOM=1-TOM

#Hierarchical clustering of the genes based on the TOM dissimilarity measure
geneTree = flashClust(as.dist(dissTOM),method="average")
#Plot the resulting clustering tree (dendrogram)

#pdf("genetree.pdf",width = 10, height = 10)

plot(geneTree,  sub="",cex=0.0013)
geneTree$labels <- colnames(datExpr)
plot(geneTree,sub="")
dev.off()

minModuleSize <- 30
ds <- 2
cutHeight <- 0.98
dthresh <- 0.15 

# Module identification using dynamic tree cut

#dynamicMods = cutreeDynamic(dendro = geneTree,  method="tree", minClusterSize = minModuleSize);
dynamicMods = cutreeHybrid(dendro = geneTree,cutHeight = cutHeight,distM = dissTOM,
                           deepSplit = 1, pamRespectsDendro = TRUE,pamStage = TRUE, 
                           minClusterSize = minModuleSize);


table(dynamicMods[[1]])

merged <- mergeCloseModules(exprData = datExpr,colors =dynamicMods$labels,cutHeight = dthresh)
summary(merged)

mColorh <- labels2colors(merged$colors)

mLabelh <- "DS=2,MMS=200,DCOR=0.15"

table(mColorh)
plotDendroAndColors(dendro = geneTree, mColorh, groupLabels = mLabelh,
                    addGuide=TRUE,dendroLabels=FALSE,main="Dendrogram With Defined Modules")
#the following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
#table(tree1)
dynamicColors = labels2colors(merged$colors)
table(dynamicColors)

plotDendroAndColors(geneTree, dynamicColors,
                    "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.2, main = "Gene dendrogram and module colors")


restGenes= (dynamicColors != "grey")
diss1=1-TOMsimilarityFromExpr(datExpr[,restGenes], power = softPower)
collectGarbage()

par(mfrow = c(1,2))
colnames(diss1) =rownames(diss1) = datExpr[restGenes]
hier1=flashClust(as.dist(diss1), method="average" )
plotDendroAndColors(hier1, dynamicColors[restGenes], 
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, 
                    guideHang = 0.2, main = "HCLUST_Gene dendrogram and module colors")
abline(h=5,col="red")

diag(diss1) = NA
#sizeGrWindow(7,7)
collectGarbage()
png("TOMplot1.png",width =20,height=20,unit="cm",res=500)
TOMplot(diss1, hier1, as.character(dynamicColors[restGenes]))
dev.off()
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module=datExpr[which(dynamicColors==color)]
  write.table(t(module), paste("results/WGCNA/tmod_",color, ".txt",sep=""), sep="\t", 
              row.names=TRUE, col.names=TRUE,quote=FALSE)
}

module.order <- tapply(1:ncol(datExpr),as.factor(dynamicColors),I) %>% .[names(.)!="grey"] %>% unlist()

m<-t(t(datExpr[,module.order])/apply(datExpr[,module.order],2,max))
#m<-t(apply(t(datExpr[,module.order]),2,scale))
names1 <- rownames(datExpr)
rownames(m) <- names1
heatmap(t(m),zlim=c(0,1),col= redWhiteGreen(256),Rowv=NA,Colv=NA,labRow=NA,scale="row",RowSideColors=dynamicColors[module.order])
#######################
library(ComplexHeatmap)
library(circlize)
p2  <-   Heatmap(t(m), name = "Row Z score"
                 #col=colorRamp2(c(0.2,0.5,1), c("yellow","#fa5f5f","white","darkblue","blue","blue","#ffd14c"))
                 ,col=colorRamp2(c(0.8,0.85,0.95,1),c("darkblue","blue","yellow","red"))
                 ,cluster_rows = FALSE
                 ,clustering_distance_rows = "euclidean"
                 ,clustering_method_rows = "complete"
                 ,cluster_columns = FALSE
                 ,clustering_distance_columns = "euclidean"
                 ,clustering_method_columns = "complete"
                 #,rect_gp = gpar(col = "black"
                 ,na_col = "white"
                 ,km =1,
                 column_names_side = "top",
                 column_names_rot = 60,
                 #,combined_name_fun = NULL
                 #,heatmap_legend_param = list(at=c(0.3,1),color_bar = "continuous")
                 ,show_row_names = FALSE
                 ,show_row_dend = FALSE
                 #,column_names_gp = gpar(cex=0.1, font=0.1, col= "green"),
                 ,row_names_gp =  gpar(fontsize=12,fontcolor= "red"),
                 column_names_gp = gpar(fontsize = 10, fontface="bold")
                 ,show_column_names = TRUE)

module_color <- as.data.frame(dynamicColors[module.order])
colnames(module_color) <- "color"
c <- list(color=unique(module_color$color))
names(c$color) <- unique(module_color$color)
ha_row = rowAnnotation(df = module_color, col= c,width = unit(1, "cm"))
p2+ha_row
################################################################################################################
#Quantifying module–trait associations
################################################################################################################
# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, dynamicColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits$group, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = gsub("ME","",names(MEs)),
               colorLabels = TRUE,
               colors = blueWhiteRed(300),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
Heatmap(moduleTraitCor,
        name = "Cor",
        column_title="Module-trait relationships",
        rect_gp = gpar(col = "black"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf(textMatrix[i, j]), x, y, gp = gpar(fontsize = 10))})

#Gene relationship to trait and important modules: Gene Significance and Module Membership
#We quantify associations of individual genes with our trait of interest by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. For each module, we also define a quantitative
#measure of module membership MM as the correlation of the module eigengene and the gene expression profile. This
#allows us to quantify the similarity of all genes on the array to every module.
# Define variable weight containing the weight column of datTrait
group = as.data.frame(datTraits$group);
names(group) = "group"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

#Figure 1: Module-trait associations. Each row corresponds to a module eigengene, column to a trait. Each cell
#contains the corresponding correlation and p-value. The table is color-coded by correlation according to the color
#legend.
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, group, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(group), sep="");
names(GSPvalue) = paste("p.GS.", names(group), sep="");


#Intramodular analysis: identifying genes with high GS and MM
#Using the GS and MM measures, we can identify genes that have a high significance for weight as well as high module
#membership in interesting modules. As an example, we look at the brown module that has the highest association
#with weight. We plot a scatterplot of Gene Significance vs. Module Membership in the brown module:
module = unique(dynamicColors)
pdf("results/results/WGCNA/Gene_significance_vs_Module_membership.pdf")
gs.mm = list()
for (i in 1:length(module)){
  print(module[i])
  column = match(module[i], modNames);
  moduleGenes = dynamicColors==module[i];
  #sizeGrWindow(7, 7);
  #par(mfrow = c(1,1));
  #png(paste(module[i],"GS","MM",".png",sep="_"),height = 12,width = 12,res=400,units = "cm")
  print(verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]), 
                           abs(geneTraitSignificance[moduleGenes, 1]),
                           xlab = paste("Module Membership"),
                           ylab = "Gene significance",
                           main = paste(module[i],"\n"),
                           cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module[i]))
}
dev.off()
################################################################################
#Module genes heatmap
pdf("results/WGCNA/Module_heatmap.pdf")
whichmodule = c("darkred","green","tan")
for (i in c(1:length(whichmodule))){
  
  ME= MEs[,paste("ME",whichmodule[i],sep = "")]
  barplot(ME,col = whichmodule[i],main=whichmodule[i])
  Emp_NCoR1_Normalized_count[rownames(datExpr[, dynamicColors == whichmodule[i]] %>% t()),] %>% write.table(., file=paste0("/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/results/WGCNA/tmod_",whichmodule[i],"_normalized_count.txt"),sep = "\t",quote = FALSE)
  datExprModule = datExpr[, dynamicColors == whichmodule[i]]
  # set the margins of the graphics window
  par(mfrow = c(1, 1), mar = c(0.3, 5.5, 3, 2))
  # create a heatmap whose columns correspond to the arrays and whose rows
  # correspond to genes
  #plotMat(t(scale(datExprModule)), cex.axis = 2, nrgcols = 30, rlabels = T,clabels = col_label, rcols = whichmodule,main = paste("heatmap", whichmodule, "module"))
  p2  <-   Heatmap(t(scale(datExprModule)), name = "Row Z score"
                   #,col=colorRamp2(c(-1,0,1), c("yellow","#fa5f5f","white","darkblue","blue","blue","#ffd14c"))
                   ,col=colorRamp2(c(-1,0,1),c("blue","white","red"))
                   ,cluster_rows = FALSE
                   ,clustering_distance_rows = "euclidean"
                   ,clustering_method_rows = "complete"
                   ,cluster_columns = FALSE
                   ,column_names_side = "top"
                   ,column_names_rot = 45
                   ,clustering_distance_columns = "euclidean"
                   ,clustering_method_columns = "complete"
                   #,rect_gp = gpar(col = "black"
                   ,na_col = "white"
                   ,km =1
                   #,combined_name_fun = NULL,
                   ,heatmap_legend_param = list(at=c(-1,0,1),color_bar = "continuous")
                   ,show_row_names = FALSE
                   ,show_row_dend = FALSE
                   #,column_names_gp = gpar(cex=0.1, font=0.1, col= "green"),
                   ,row_names_gp =  gpar(fontsize=12,fontcolor= "red")
                   ,show_column_names = TRUE)
  print(p2)
}
dev.off()
#cor.test(MEs$MEblue,MEs$MEturquoise,method = "spearman")
################################################################################################################



#To get a sense of how related the modules are one can summarize each module by its eigengene (first principal component).
datME=MEs
signif(cor(datME, use="p"), 2)
#We define a dissimilarity measure between the module eigengenes that keeps track of the sign of the correlation
#between the module eigengenes, and use it to cluster the eigengene:
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")

sizeGrWindow(8,9)
y <- datTraits$group
plotMEpairs(datME,y=y)

## Gene significance
GS1=as.numeric(cor(y,datExpr, use="p"))
GeneSignificance=abs(GS1)
ModuleSignificance=tapply(GeneSignificance, dynamicColors, mean, na.rm=T)

#To plot module significance, one can use
sizeGrWindow(8,7)
par(mfrow = c(1,1))
plotModuleSignificance(GeneSignificance,dynamicColors)


##################################################################################################################
##################################################################################################################
################################################################################################################
#Intramodular connectivity, module membership, and screening for intramodular hub genes
################################################################################################################
# We begin by calculating the intramodular connectivity for each gene. 
# (In network literature, connectivity is often #referred to as ”degree”.) 
# The function intramodularConnectivity computes the whole network connectivity kTotal,
# the within module connectivity kWithin, kOut=kTotal-kWithin, and kDiff=kIn-kOut=2*kIN-kTotal

#Intramodular connectivity
ADJ1=abs(cor(datExpr,use="p"))^12
Alldegrees1=intramodularConnectivity(ADJ1, dynamicColors)
head(Alldegrees1)

datKME=signedKME(datExpr, MEs, outputColumnName="MM.")  

#Relationship between gene significance and intramodular connectivity
colorlevels=unique(dynamicColors)
sizeGrWindow(9,6)
par(mfrow=c(2,as.integer(0.5+length(colorlevels)/2)))
par(mar = c(4,5,3,1))
pdf("GeneSignificance_IntramodularConnectivity.pdf")
for (i in c(1:length(colorlevels)))
{
  whichmodule=colorlevels[[i]];
  restrict1 = (dynamicColors==whichmodule);
  #png(paste(whichmodule,"GS","MC",".png",sep="_"),height = 10,width = 10,res=400,units = "cm")
  verboseScatterplot(Alldegrees1$kWithin[restrict1],cex.main = 1,cex = 0.5,
                           GeneSignificance[restrict1], col=dynamicColors[restrict1],
                           main=whichmodule,
                           xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
  
}
dev.off()

#Finding genes with high gene significance and high intramodular connectivity in specific modules
#abs(GS1)>.8 #adjust parameter based on actual situations
#abs(datKME$MM.black)>.8 #at least larger than 0.8
FilterGenes= GeneSignificance >0.8 & abs(datKME$MM.green)>0.7
table(FilterGenes)

hubgenes <- rownames(datKME)[FilterGenes]
hubgenes
write.table(hubgenes,file="black_hub_genes",quote = FALSE,row.names = FALSE)

hubs <- chooseTopHubInEachModule(
  datExpr, 
  mColorh, 
  omitColors = "grey", 
  power = 2, 
  type = "signed hybrid")
################################################################################################################
################################################################################################################




################################NETWORK_EXPORT####################################################
TOM=TOMsimilarityFromExpr(datExpr,
                          networkType = "signed hybrid", 
                          TOMType = "unsigned", 
                          power = softPower)
modules = c("green", "red")
modules = names(table(dynamicColors)[c(-3,-6,-13)])
probes = names(datExpr)
inModule = is.finite(match(dynamicColors, modules))
modProbes = probes[inModule] # This should probably be modProves = names(probes)[inModule]

#modProves = names(probes)[inModule]

modTOM = TOM[inModule, inModule]
dimnames(modTOM) = list(modProbes,modProbes)
dim(modTOM)
cyt = exportNetworkToCytoscape(modTOM, 
                               #edgeFile = paste("results/WGCNA/tmod_CytoscapeInput-edges_signed-",paste(modules, collapse="-"), ".txt", sep=""), 
                               #nodeFile = paste("results/WGCNA/tmod_CytoscapeInput-nodes_signed-",paste0(modules, collapse="-"), ".txt",sep=""),
                               weighted = TRUE, 
                               threshold = 0.3, nodeNames = modProbes, 
                               nodeAttr = dynamicColors[inModule]);






#https://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/Tutorials/anRichment-Tutorial1.pdf
options(stringsAsFactors = FALSE)
library("anRichment")
GOcollection = buildGOcollection(organism = "mouse")
symbol.0 = probes;
moduleColor = dynamicColors;
table(moduleColor)
entrez = convert2entrez(organism = "mouse", symbol = symbol.0);
table(is.finite(entrez))
#GOcollection = buildGOcollection(organism = "human")

GOenrichment = enrichmentAnalysis(
  classLabels = moduleColor, identifiers = entrez,
  refCollection = GOcollection,
  useBackground = "given",
  threshold = 1e-4,
  thresholdType = "Bonferroni",
  getOverlapEntrez = TRUE,
  getOverlapSymbols = TRUE,
  ignoreLabels = "grey",
  maxReportedOverlapGenes=1000);
table.display = GOenrichment$enrichmentTable;
table.display$overlapGenes = shortenStrings(table.display$overlapGenes, maxLength = ,
                                            split = "|");
head(table.display)
names(GOenrichment$dataSetDetails)
write.csv(GOenrichment$enrichmentTable, file = "GOenrichment-enrichmentTable_1.csv",
          row.names = FALSE,quote = FALSE);



active = entrez[dynamicColors=="yellow"];
all = entrez;
GOenrichment.brown = enrichmentAnalysis(active = active, inactive = all,
                                        refCollection = GOcollection,
                                        useBackground = "given",
                                        threshold = 1e-4,
                                        thresholdType = "Bonferroni",
                                        maxReportedOverlapGenes=1000,
                                        getOverlapSymbols = TRUE);

table.display <- GOenrichment.brown$enrichmentTable
View(table.display)
y <- table.display[c(1,2,3,8,9),]  #brown module
y <- table.display[c(1,2,3,4,5),] #green yellow
y <- table.display[c(1,2,4,5,10),]
ggplot(y, # you can replace the numbers to the row number of pathway of your interest
       aes(y = -log(pValue), x = reorder(dataSetName, -log(pValue)))) + 
  geom_bar(stat="identity") +
  coord_flip()+
  theme_bw(base_size = 14) +
  #scale_colour_gradient(limits=c(0, 0.10), low="red") +
  xlab("-log(p-value")+
  ggtitle("GO BP enrichment") +
  theme(axis.text.x=element_text(size=12,color="black")
        , axis.title.x=element_text(size=12)
        , axis.text.y=element_text(size=12,color="black")
        , axis.title.y=element_blank(),
        ,legend.justification=c(0,1)
        ,legend.box.background = element_rect(colour = "black"))
##########################################################################################################
names(probes) = dynamicColors
probes.list = split(probes,names(probes)) 
probes.list =lapply(probes.list, unname)
probes.list_toupper =lapply(probes.list, toupper) 

term2gene = read.gmt("/home/imgsb/tools/GSEA/c2.cp.reactome.v7.1.symbols.gmt")
compEnricher =compareCluster(probes.list_toupper, fun="enricher",TERM2GENE=term2gene)
write.table(file = "results/WGCNA/Module_Reactome_pathway.csv",compEnricher,sep = "\t",quote = FALSE,row.names = FALSE)

# Pathway barplot
library(tidytext)
compEnricher %>%
  as.data.frame() %>% 
  group_by(Cluster) %>%
  slice(1:6) %>%
  ungroup %>% 
  mutate(Description = gsub("REACTOME|_"," ",Description),
         Description = reorder_within(Description, -pvalue, Cluster)) %>% 
  
  as.data.frame() %>% 
  #filter(Cluster %in% c(1,2,6,7)) %>%
  #.[c(-4,-5,-7,-8,-10),] %>%
  #Description = reorder(Description, pvalue))
  ggplot(aes(x=reorder(Description,-log(pvalue)), y=-log(pvalue), fill = p.adjust)) +
  geom_bar(stat="identity") +
  coord_flip() +
  facet_wrap(~Cluster, scales = "free",ncol = 1) +
  #scale_x_reordered() +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 35))+
  scale_fill_gradient(low="red",high="blue")+
  scale_y_continuous(expand = c(0,0)) +
  gg_theme
########################################################################################################
save(datExpr,MEs, dynamicColors, geneTree,
     file = "NCoR1_SMRT-networkConstruction-auto.RData")



########################################################################################################
# Exporting results of the network analysis. We now put together a data frame that summarizes the results of network analysis, 
# namely the gene significances(GS) and module memberships (also known as kME) of all probes.  


consMEs.unord = multiSetMEs(datExpr, universalColors = moduleLabels, excludeGrey = TRUE)
GS = list();
kME = list();
for (set in 1:nSets)
{
  GS[[set]] = corAndPvalue(multiExpr[[set]]$data, Traits[[set]]$data);
  kME[[set]] = corAndPvalue(multiExpr[[set]]$data, consMEs.unord[[set]]$data);
}
# We  perform  a  very  simple  ”meta-analysis”  by  combining  the  Z  scores  of  correlations  
# from  each  set  to  form  a meta-Z score and the corresponding p-value.
GS.metaZ = (GS[[1]]$Z + GS[[2]]$Z)/sqrt(2);
kME.metaZ = (kME[[1]]$Z + kME[[2]]$Z)/sqrt(2);
GS.metaP = 2*pnorm(abs(GS.metaZ), lower.tail = FALSE);
kME.metaP = 2*pnorm(abs(kME.metaZ), lower.tail = FALSE);

# Next we form matrices holding the GS and kME. 
# We use a simple re-shaping trick to put the values and theassociated p-values and meta-analysis results next to one another.
GSmat = rbind(GS[[1]]$cor, GS[[2]]$cor, GS[[1]]$p, GS[[2]]$p, GS.metaZ, GS.metaP);
nTraits = checkSets(Traits)$nGenes
traitNames = colnames(Traits[[1]]$data)
dim(GSmat) = c(nGenes, 6*nTraits)
rownames(GSmat) = probes;
colnames(GSmat) = spaste(
  c("GS.set1.", "GS.set2.", "p.GS.set1.", "p.GS.set2.", "Z.GS.meta.", "p.GS.meta"),rep(traitNames, rep(6, nTraits)))
# Same code for kME:
kMEmat = rbind(kME[[1]]$cor, kME[[2]]$cor, kME[[1]]$p, kME[[2]]$p, kME.metaZ, kME.metaP);
MEnames = colnames(consMEs.unord[[1]]$data);
nMEs = checkSets(consMEs.unord)$nGenes
dim(kMEmat) = c(nGenes, 6*nMEs)
rownames(kMEmat) = probes;
colnames(kMEmat) = spaste(
  c("kME.set1.", "kME.set2.", "p.kME.set1.", "p.kME.set2.", "Z.kME.meta.", "p.kME.meta"),
  rep(MEnames, rep(6, nMEs)))

# Finally we put together the full information data frame and write it into a plain text CSV file that can be read 
# bystandard spreadsheet programs.  Note that the probes are not sorted in any particular way; many sort orders 
# are possible and we leave it to the reader to either modify the code or to perform the sort in a spreadsheet software.

info = data.frame(Probe = probes, 
                  GeneSymbol = annot$gene_symbol[probes2annot],
                  EntrezID = annot$LocusLinkID[probes2annot],
                  ModuleLabel = moduleLabels,
                  ModuleColor = labels2colors(moduleLabels),
                  GSmat,kMEmat);
write.csv(info, file = "consensusAnalysis-CombinedNetworkResults.csv",row.names = FALSE, quote = FALSE)

##########################################################################################################

TFs = read.csv("../../Gyan/mm10/Mus_musculus_TF_list.txt",header = F) 

Module_TF = merge(TFs,cyt[[2]],by.x="V1",by.y="nodeName" ) %>% 
                    filter(V1 %in% NCoR1_LRT_cluster_df[which(NCoR1_LRT_cluster_df$cluster %in% c(1,2,3,6,7,10,11)),]$genes)
Module_TF_clusters_df.exp = Emp_NCoR1_vst_rep_merged[Module_TF$V1,]
Module_TF_clusters_df.exp_scaled= t(apply(Module_TF_clusters_df.exp, 1, scale))
colnames(Module_TF_clusters_df.exp_scaled) = colnames(Module_TF_clusters_df.exp)

draw(Heatmap(Module_TF_clusters_df.exp_scaled[,c(1,6,2,7,3,8,4,9,5,10)],
             name="Z-score",
             k=1,
             split=Module_TF$`nodeAttr[nodesPresent, ]`,
             heatmap_legend_param=list(at=c(-1.5,0,1.5),color_bar="continuous", 
                                       legend_direction="horizontal", legend_width=unit(5,"cm"),legend_height=unit(3,"cm"),
                                       title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
             #rect_gp = gpar(col = "grey"),
             cluster_columns = FALSE,
             cluster_rows = FALSE,
             column_names_side = "top",
             column_names_rot = 70,
             row_names_gp = gpar(fontsize=7),
             column_names_gp = gpar(fontsize = 10)),
     heatmap_legend_side="bottom")

Module_TF_target = list()
target = ""
for (i  in Module_TF$V1) {
  print(i)
  c=1
  print(c)
  target=""
  cyt_target = cyt[[1]] %>% filter(i == cyt[[1]]$fromNode | i == cyt[[1]]$toNode)
  for (j in 1:dim(cyt_target)[1]) {
    
      target[c] = ifelse(cyt_target[j,]$fromNode ==i, cyt_target[j,]$toNode, cyt_target[j,]$fromNode)
      
      c=c+1            
      }
      target =  na.omit(unique(target))
      Module_TF_target[[i]] = target[target != ""]
}

Module_TF_target_vs_DE = newGOM(Emp_DE.list[c(3,5,7)],Module_TF_target,genome.size=gs.RNASeq)

Module_TF_target_vs_DE_pval = getMatrix(Module_TF_target_vs_DE, name="pval") %>% 
                                      t() %>% 
                                      as.data.frame()
Module_TF_target_vs_DE_pval = Module_TF_target_vs_DE_pval %>%
                              rownames_to_column("TF") %>%
                              mutate(EC_vs_EU_up_Rank= rank(-EC_vs_EU_up),
                                     EP_vs_EU_up_Rank=rank(-EP_vs_EU_up),
                                     EA_vs_EU_up_Rank=rank(-EA_vs_EU_up)) #%>%
                              mutate(Rank= rowMeans(cbind(EC_vs_EU_Rank,EP_vs_EU_Rank,EA_vs_EU_Rank)))
Module_TF_target_vs_DE_pval$Total_Target = lengths(Module_TF_target_vs_DE@gsetB)

#Total DE targets 
inter.nl <- getNestedList(Module_TF_target_vs_DE, name="intersection")
Total_TF_hits = lapply(inter.nl,lengths) %>% do.call(rbind,.) %>% as.data.frame()
Module_TF_target_vs_DE_pval = merge(Module_TF_target_vs_DE_pval,Total_TF_hits,by.x="TF",by.y=0)
write.table(Module_TF_target_vs_DE_pval,
          file="/home/imgsb/Gyan/NCOR1/Emp_NCoR1_CpG_pIC_C+P_DE_analysis/results/WGCNA/Module_TF_target_vs_DE_pval.tsv",
          sep = "\t",quote = FALSE,row.names = FALSE)

Heatmap(-log10(Module_TF_target_vs_DE_pval[,c(2,3,4)]+10^-381),
        col= colorRamp2(c(0,10,100),c("grey","yellow","red")))
        #heatmap_legend_param=list(at=c(0,5,10,15,20),color_bar="continuous", 
         #                         legend_direction="horizontal", legend_width=unit(5,"cm"),
          #                        title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),)


green_module_NP_vs_EP= Module_TF_clusters_df.exp %>% 
                              filter(rownames(.) %in% probes.list$green &
                                     rownames(.) %in% DE.list$NP_vs_EP_up |
                                     rownames(.) %in% DE.list$NP_vs_EP_down) %>% 
                              dplyr::select(1,6,2,7,3,8,4,9,5,10)
green_module_NC_vs_EC= Module_TF_clusters_df.exp[Module_TF_target_vs_DE_pval$TF[Module_TF_target_vs_DE_pval$TF %in% probes.list$green & 
                                                 Module_TF_target_vs_DE_pval$TF %in% Emp_DE.list$NC_vs_EC],][,c(1,6,2,7,3,8,4,9,5,10)] %>% 
                                      filter(!(rownames(.) %in% rownames(green_module_NP_vs_EP)))
draw(Heatmap(green_module_NP_vs_EP,
             name="vst expression",
             k=2,
             #rect_gp = gpar("grey"),
             #split=Module_TF$`nodeAttr[nodesPresent, ]`,
             col=colorRamp2(c(5,10,15),c("blue","white","red")),
             heatmap_legend_param=list(at=c(5,10,15),color_bar="continuous", 
                                       legend_direction="horizontal", legend_width=unit(5,"cm"),legend_height=unit(3,"cm"),
                                      title_position="topcenter", title_gp=gpar(fontsize=10, fontface="bold")),
             #rect_gp = gpar(col = "grey"),
             cluster_columns = FALSE,
             cluster_rows = TRUE,
             column_names_side = "top",
             column_names_rot = 70,
             row_names_gp = gpar(fontsize=10),
             column_names_gp = gpar(fontsize = 10)),
              heatmap_legend_side="bottom")

