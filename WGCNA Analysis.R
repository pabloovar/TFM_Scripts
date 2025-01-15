# 2024 - Pablo Vargas Rodríguez
# RNA seq analysis of oligodendrocytes WT and KO for Cort-/-

# WGCNA analysis, this analysis is based in this two: 
# https://bioinformaticsworkbook.org/tutorials/wgcna.html#gsc.tab=0
# https://fuzzyatelin.github.io/bioanth-stats/module-F21-Group1/module-F21-Group1.html#Data_Cleaning_and_Removing_Outliers



library(WGCNA)
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(pheatmap)
library(GeneOverlap)
library(enrichplot)
library(patchwork)
library(corrplot)
library(RColorBrewer)
library(viridis)

setwd('C:/Users/usuario/Desktop/RESULTADOS/RNAseq_OLG/')
WGCNA_dir <- 'WGCNA'

head(data)
data <- read.csv('datasets/Raw_counts_genes.txt', sep = '\t') #

rownames(data) <- data$Gene_ID
data$Gene_ID <- NULL
data$Length <- NULL

data <- data[rowSums(data) > 0, ]

data[1:5,]


#Identifying Outiers Gene

goodGenes <- goodSamplesGenes(t(data))

summary(goodGenes)
goodGenes$allOK #TRUE

data <- read.csv('datasets/norm_counts.txt', sep = '\t')

#SOFT Threshold

softPower <- pickSoftThreshold(norm_counts, networkType = 'signed')  #10
write_tsv(softPower$fitIndices, file.path(WGCNA_dir,'pickSoftThreshold.tsv'))

#Threshold graphs---------------------------------------------------------------------------------------------
pdf(file.path(WGCNA_dir,'soft_threhold.pdf'))

plot(softPower$fitIndices[,1], softPower$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(softPower$fitIndices[,1],softPower$fitIndices[,2],col="red")
abline(h=0.80,col="red")


plot(softPower$fitIndices[,1], softPower$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(softPower$fitIndices[,1], softPower$fitIndices[,5], labels= softPower$fitIndices[,1],col="red")

dev.off()



allowWGCNAThreads() #Allow use more threads of the system


# Convert norm_counts to numeric values using sapply
norm_counts[] <- sapply(norm_counts, as.numeric)

# Define the soft-thresholding power
power <- 10

# Calculate the dissimilarity matrix using TOM (Topological Overlap Matrix) similarity from adyacency matrix
TOMdissim <- 1 - TOMsimilarity(adjacency(norm_counts, power = power))

geneTree <- hclust(as.dist(TOMdissim), method = 'average') #hierarchical clustering dendrogram based on dissimilarity

Modules <- cutreeDynamic(dendro = geneTree, distM = TOMdissim, deepSplit = 2, pamRespectsDendro = FALSE,
                         minClusterSize = 30) # minimum size of a module



moduleColors <- labels2colors(Modules) # Assign a unique color to each module


writeLines(moduleColors, file.path(WGCNA_dir, 'modulecolors.txt')) #save



#Module fusion---------------------------------------------------------------------------------------------------------
ME.dissim = 1-cor(MElist$eigengenes, use="complete") #corr between ME
METree = hclust(as.dist(ME.dissim), method = "average")

pdf(file.path(WGCNA_dir,'MEplot.pdf'))
par(cex = 0.3)
plot(METree)
abline(h=.25, col = "red")
abline(h=.3, col = "green")
dev.off()

#merge modulos with 70% correlation
merge <- mergeCloseModules(norm_counts, moduleColors, cutHeight = 0.3) #altura de 0.25 se corresponde con una correlación de 0.75

Merged_colors =merge$colors #lista con los modulos a los que se asigna cada gen
write.table(Merged_colors, file.path(WGCNA_dir,"moduleColors_merged.txt"), sep = "\t")

module_eigengens = merge$newMEs #lista con ME de cada muestra
write.table(Merged_MEs, file.path(WGCNA_dir,"moduleEigenges_merged.txt"), sep = "\t")

pdf(file.path(WGCNA_dir,'ME_merged_plot.pdf'))
plotDendroAndColors(geneTree, cbind(moduleColors, Merged_colors),
                    c("Original Module", "Merged Module"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors for original and merged modules (70%cor)")
dev.off()

# TRAIT Modules correlation ----------------------------------------------------
pheno <- data.frame(
  Genotype = c(rep("ko", 4), rep("wt", 6)),
  Condition = c(rep("differentiation", 2), rep("proliferation", 2), rep("differentiation", 3), rep("proliferation", 3)))  # Segunda anotación: wt vs ko 
nSamples <- nrow(norm_counts)
nGenes <- ncol(norm_counts)


#Construcción de dataframe con las condiciones de cada muestr
sampletype <- c('kodif', 'kodif', 'koprol', 'koprol', 'wtdif', 'wtdif', 'wtdif', 'wtprol', 'wtprol',  'wtprol')

trait <- cbind(pheno,sampletype )
trait$sampletype <- factor(trait$sampletype, levels = c('wtprol', 'koprol','wtdif','kodif' ))


diff <- binarizeCategoricalColumns(trait$sampletype,
                                   includePairwise = FALSE,
                                   includeLevelVsAll = TRUE,
                                   minCount = 1)

traita <- cbind(trait, diff) 

traita <- traita %>% mutate(Diff = ifelse(grepl('differentiation', Condition), 1, 0)) %>%
  mutate(KO = ifelse(grepl('ko', Genotype), 1, 0)) %>%
  dplyr::select(-c(Genotype, Condition, sampletype))

traita

module_eigengens
traita

# correlation 
module_trait_correlation <- cor(module_eigengens, traita, use = 'p')

module_trait_correlation_pval <- corPvalueStudent(module_trait_correlation, nSamples)
module_trait_correlation_pval

# just significant
signif_modules <- apply(module_trait_correlation_pval, 1, function(row) any(row<0.05)) #just 18!!

module_trait_correlation <- module_trait_correlation[signif_modules, , drop = FALSE]
module_trait_correlation_pval <- module_trait_correlation_pval[signif_modules, , drop = FALSE]

signif_modules_v <- names(signif_modules[signif_modules, drop = FALSE]) %>%sub("^ME", "", .)  
signif_modules_v

#heatmap

rd_bu_palette <- colorRampPalette(brewer.pal(11, "RdBu"))

pdf(file = 'WGCNA/corr_plot.pdf', width = 15, height = 10)

corrplot(module_trait_correlation, p.mat = module_trait_correlation_pval, insig = 'label_sig', pch.cex = 0.9,
         sig.level = c(0.001, 0.01, 0.05),method = 'color',
         tl.col = 'black', col = rev(rd_bu_palette(200)),
)

dev.off()


signif_modules_v

writeLines(signif_modules_v, 'WGCNA/signif_modules.txt')

MElist

#eigengens Network

plotEigengeneNetworks(Merged_MEs, "", plotHeatmap = F,
                      marDendrograms = c(3, 4, 2, 1), colorLabels = T)


# Save txt with gene list of each module and variables----------------------------------------------------------
Merged_colors <- read.table(file.path(WGCNA_dir,"moduleColors_merged.txt"))
Merged_colors <- as.vector(t(Merged_colors))

Merged_MEs <- read.table(file.path(WGCNA_dir,"moduleEigenges_merged.txt"))

genes <- colnames(norm_counts)
gen_color <- data.frame(genes, Merged_colors) 
genesbyColor <- gen_color %>% group_by(Merged_colors)%>% summarise(genes = list(genes)) # create dataframe with genes by module


color_dir <- 'WGCNA/Modules'
for (i in 1:nrow(genesbyColor)){
  module <- genesbyColor$Merged_colors[i] #get the module
  gene_list <- unlist(genesbyColor$genes[i]) #get the genes
  
  assign(module, gene_list)
  
  filename <- paste0(module,'_genes.txt')
  writeLines(gene_list, file.path(color_dir, filename) )
}

# RESULTs PLOTS-----------------------------------------------------------------------------------------------------------------------------------------

# Gene expressions plots

data_expr <- read.table('ReadCount/tpm.txt') #tpm


data_expr <- as.data.frame(t(data_expr))
data_expr <- data_expr[, apply(data_expr, 2, function(col) any(col != 0))]%>% #delete genes with 0 expression in all samples
  t() %>% mutate(Module = Merged_colors )

data_expr[,1:5]


module <- as.data.frame(Merged_colors)
data_plot <- cbind(data_expr,module)

norm_expression <- t(norm_counts)
norm_plot <- cbind(norm_expression,module)

sampletype_order <- c("wt_prol", "ko_prol", "wt_diff", "ko_diff")
sample_order <- c("wt_prol1","wt_prol2", "wt_prol3","ko_prol2","ko_prol3",
                  "wt_diff1", "wt_diff2","wt_diff3","ko_diff2", "ko_diff3")

summary(signif_modules_v %in% Merged_colors)
str(Merged_colors)

norm_plot_group <- norm_plot %>%
  filter(Merged_colors %in% signif_modules_v)%>% #select signif modules with correlation with traits
  group_by(Merged_colors) %>% 
  summarise(across(all_of(sample_order), mean))%>% #medium  value of expression
  pivot_longer(cols = -Merged_colors, names_to = 'Sample', values_to = 'Value') #longer for plotting

norm_plot_group$Sample <- factor(norm_plot_group$Sample, levels = sample_order) #factor to choose the order of the sample


x_positions <- c(3.5, 5.5, 8.5) #values to draw lines  (between wtprol y ko prol, etc.)

norm <- ggplot(norm_plot_group, aes(x=Sample, y= Value, group = Merged_colors, color = Merged_colors))+
  geom_line()+
  geom_point()+
  theme_minimal()+
  labs(title="Average Expression by Module Across Samples",
       y = 'Normalized Expression'
  )+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text.x = element_text(size = 16, face = 'bold'))+
  facet_wrap(~Merged_colors, scales = "free_y")+ # split by modules
  geom_vline(xintercept = x_positions)+ # draw lines
  theme(legend.position = "none")

norm

ggsave(file.path(WGCNA_dir,'Module_Expression_NORM.pdf'), norm,  
       width = 16,            
       height = 10)



#Plot of interest modules

filtered_data <- norm_plot %>%
  filter(Merged_colors %in% signif_modules_v)

pdf('WGCNA/modulos.pdf')
for (module in signif_modules_v){
  norm_plot_group <- norm_plot_group[norm_plot_group$Merged_colors == module,]
  x <- ggplot(norm_plot_group, aes(x=Sample, y= Value, group = module, color = module))+
    geom_line()+
    geom_point()+
    theme_minimal()+
    labs(title= paste("Average Expression of", module),
         y = 'Normalized Expression'
    )+
    theme(axis.text.x = element_text(angle = 90, hjust = 1),
          strip.text.x = element_text(size = 16, face = 'bold'))+
    facet_wrap(~Merged_colors, scales = "free_y")+ #split the modules
    geom_vline(xintercept = x_positions)+ #draw the lines
    theme(legend.position = "none")       
  print(x)
}

dev.off()


# Overlaps of modules with genesets of DEGs----------------------------------------------------------------------------------------------------
# https://rdrr.io/bioc/GeneOverlap/man/GeneOverlap.html
over_dir <-'results/Overlaps'

up_diff <- as.vector(t(read.table('results/UP_diff.txt')))
down_diff <- as.vector(t(read.table('results/DOWN_diff.txt')))
up_prol <- as.vector(t(read.table('results/UP_prol.txt')))
down_prol <- as.vector(t(read.table('results/DOWN_prol.txt')))

difprocess_onlyKO_UP <- as.vector(t(read.table('results/Overlaps/Differentiation_Process/UP_Diff_onlyKO.txt')))
difprocess_onlyKO_DOWN <- as.vector(t(read.table('results/Overlaps/Differentiation_Process/DOWN_Diff_onlyKO.txt')))

difprocess_onlyWT_UP <- as.vector(t(read.table('results/Overlaps/Differentiation_Process/UP_Diff_onlyWT.txt')))
difprocess_onlyWT_DOWN <- as.vector(t(read.table('results/Overlaps/Differentiation_Process/DOWN_Diff_onlyWT.txt')))

difprocess_common_UP <- as.vector(t(read.table('results/Overlaps/Differentiation_Process/UP_Diff_Common.txt')))
difprocess_common_DOWN <- as.vector(t(read.table('results/Overlaps/Differentiation_Process/DOWN_Diff_Common.txt')))

KOvsWT_common_UP<- as.vector(t(read.table('results/Overlaps/KOvsWT/UP_KOvsWT_Common.txt')))
KOvsWT_common_DOWN<- as.vector(t(read.table('results/Overlaps/KOvsWT/DOWN_KOvsWT_Common.txt')))

degs <- list(up_diff=up_diff, down_diff=down_diff,up_prol= up_prol,down_prol= down_prol,difprocess_onlyKO_UP= difprocess_onlyKO_UP, 
             difprocess_onlyKO_DOWN=difprocess_onlyKO_DOWN, difprocess_onlyWT_UP=difprocess_onlyWT_UP, difprocess_onlyWT_DOWN=difprocess_onlyWT_DOWN, 
             difprocess_common_UP=difprocess_common_UP,difprocess_common_DOWN=difprocess_common_DOWN, 
             KOvsWT_common_UP=KOvsWT_common_UP, KOvsWT_common_DOWN=KOvsWT_common_DOWN )

ordenado_geneset <- c("down_diff", "down_prol", "difprocess_onlyKO_DOWN", "difprocess_onlyWT_DOWN", "difprocess_common_DOWN", "KOvsWT_common_DOWN",
                      "up_diff", "up_prol", "difprocess_onlyKO_UP", "difprocess_onlyWT_UP", "difprocess_common_UP", "KOvsWT_common_UP")


#FUNCION DE ODSS RATIO PLOT
getOverlap <- function(list, degs, name = 'Geneset', x.title = T, combined = F, max_log10pval = NULL, max_odds_ratio = NULL, print = F, DEGS=T){
  #list = gene list of the module
  #degs = gene list to compare
  #name = axis y and graph tittle
  #x title = if x axis tittle should appear or no (important for combined graph)
  #combined. Si T, the scale are previouly calculated values to plot all equal
  #print - boolean to print the GeneOVerlap values
  #DEGS - the odss ration ir with the differents degs genesets (overlaps, UP, down etc)
  
  ods <- c() 
  pval <- c()
  for (deg in degs){ #diferentes genelist of degs
    over <- newGeneOverlap(list,deg, genome.size = 31201) #31201, total number de genes 
    over <- testGeneOverlap(over)
    ods <- append(ods,over@odds.ratio)
    pval <- append(pval,over@pval)
  }
  
  ods_df <- data.frame(Ods = ods, geneset = names(degs), pval = pval )
  ods_df <- ods_df %>% mutate(log10pval = -log10(pval))%>%
    mutate(log10pval = pmin(-log10(pval), 25)) # Trunc log10(p-value) to 25 max
  
  if (DEGS) {ods_df$geneset <- factor(ods_df$geneset, levels = ordenado_geneset)}
  
  if (print){print(ods_df)}
  
  if (DEGS) {color_eje_x <- c(rep('blue',6), rep('red',6))} else {color_eje_x <- 'black'} #6 down blue, up en red
  
  ggplot(ods_df, aes(x=geneset, y=1))+
    geom_point(aes(size= log10pval, color = Ods))+
    scale_size_continuous(name = "-log10(p-value)", range = c(2, 10), limits = c(1, max_log10pval)) + # point size by log10pval
    scale_color_gradientn(
      colors = c("white", viridis(10)),
      values = scales::rescale(c(0, 0.01, max_odds_ratio)), # white is close to 0
      name = "Odds Ratio",
      limits = c(0, max_odds_ratio),
      oob = scales::squish)+ # color gradient to Odds Ratio
    
    scale_y_continuous(expand = c(0,0))+
    labs(x = "", y = name) +
    theme_bw() +
    theme(axis.text.x = if (x.title) element_text(angle = 90, hjust = 1, color = color_eje_x) else element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_text(angle = 0, vjust = 0.5),
          aspect.ratio = 0.1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = 'top',
          legend.text = element_text(size = 14), # Make legend text larger
          legend.title = element_text(size = 14)) 
}


max_log10pval <- 75
max_odds_ratio <- 77.63283

#all modules plot
plots <- list()
for (i in 1:nrow(genesbyColor)){
  module <- genesbyColor$Merged_colors[i] #get the modulo
  gene_list <- unlist(genesbyColor$genes[i]) #get the gene
  
  assign(module, gene_list)
  if (i<nrow(genesbyColor)){plot <- getOverlap(gene_list, degs, name = module, x.title = F, combined = T,
                                               max_log10pval = max_log10pval, max_odds_ratio = max_odds_ratio)}
  #just the last module have the x axis
  else {plot <- getOverlap(gene_list, degs, name = module, combined = T,
                           max_log10pval = max_log10pval, max_odds_ratio = max_odds_ratio)}
  
  plots[[i]] <- plot
}

combined_plot <- wrap_plots(plots, ncol = 1) + plot_layout(guides = "collect")

combined_plot

pdf('WGCNA/Overlaps_ModuleDEGs.pdf', height = 16, width = 8)
print(combined_plot)

dev.off()


#PLOT DEGS Overlaps WITH Cell Types-------------------------------------------------------

cellular_type_gsea_list <- tapply(cellular_type_gsea$gene_symbol, cellular_type_gsea$gs_name, c)
cellular_type_gsea_list

plots_CELLTYPE <- list()
i=0
for (i in 1:nrow(genesbyColor)){
  module <- genesbyColor$Merged_colors[i] #get the module
  gene_list <- unlist(genesbyColor$genes[i]) #get the genes
  print(module)
  assign(module, gene_list)
  if (i<nrow(genesbyColor)){plot <- getOverlap(gene_list, cellular_type_gsea_list, name = module, x.title = F, combined = T,
                                               max_log10pval = 10, max_odds_ratio = 15, print = T, DEGS = F)}
  #just the last module have the x axis
  else {plot <- getOverlap(gene_list, cellular_type_gsea_list, name = module, combined = T,
                           max_log10pval = 10, max_odds_ratio = 15, print = T, DEGS = F)}
  
  plots_CELLTYPE[[i]] <- plot
}

combined_plot <- wrap_plots(plots_CELLTYPE, ncol = 1) + plot_layout(guides = "collect")

combined_plot

pdf('WGCNA/Overlaps_Modules_CELLTYPES.pdf', height = 25, width = 10)
print(combined_plot)
dev.off()


# Module membership Calculation (cor: gen and ME)---------------------------------------------------------------
module_names = substring(names(Merged_MEs), 3)

geneModuleMembership <- cor(module_eigengens, norm_counts, use ='p')
geneModuleMembership_pval <- corPvalueStudent(geneModuleMembership, nSamples) #pval matrix

which.max(geneModuleMembership[,'Cort']) #MEcyan 0.70

geneModuleMembership[61,'Cort'] #skyblue3 0.68

gene_corr <- cor(norm_counts, traita$data.KO.vs.all)
gene_corr_pval <- corPvalueStudent(gene_corr, nSamples)

gene_corr_pval%>% as.data.frame() %>% arrange(V1) %>% head(25)


cort_module <- norm_plot %>%
  filter(moduleColors == 'MEcyan')
cort_module

print(rownames(cort_module))

interest_modules <- genesbyColor %>% filter(sapply(genes, function(gene_list) 'Cort' %in% gene_list)) #esta en lightsteelblue1 
interest_modules$Merged_colors
blueMM <- geneModuleMembership_pval['MElightsteelblue1',] 
blueMM <- data.frame(MM = blueMM) 
blueMM %>% arrange(MM) %>% tail(100)

blueMM['Cort',]

