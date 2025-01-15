# 2024 - Pablo Vargas Rodríguez
# RNA seq analysis of oligodendrocytes WT and KO for Cort-/-
# basic analysis and functional enrichment 


library(readr)
library(ggplot2)
library(VennDiagram) #calculos and overlaps
library(ggVennDiagram) #nice venns plots
library(ggrepel)
library(magrittr) #usar los %>%
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork) 
library(enrichplot)
library(stringr) 
library(pheatmap)
library(ggbreak)
library(viridis)
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)

#settings the paths
setwd('C:/Users/usuario/Desktop/RESULTADOS/RNAseq_OLG/')

plots_dir <- 'plots'
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}
results_dir <- 'results'
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

results_dir_go <- 'results/GO'
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

results_dir_overlaps <- 'results/Overlaps'
if (!dir.exists(results_dir_overlaps)) {
  dir.create(results_dir_overlaps)
}

set.seed(2024)

# data import - output files from reanalyzerGSE
koDvswtD <- read_tsv('datasets/kodiffvswtdiff.txt')
koPvswtP <- read_tsv('datasets/koprolvswtprol.txt')
koDvskoP <- read_tsv('datasets/kodiffvskoprol.txt')
wtDvswtP <- read_tsv('datasets/wtdiffvswtprol.txt')

# Sign change, Cort gene like a downregulated gene
koDvswtD <- koDvswtD %>% mutate(logFC__ko_diff_VS___wt_diff= -logFC__ko_diff_VS___wt_diff)
koPvswtP <- koPvswtP %>% mutate(logFC__ko_prol_VS___wt_prol= -logFC__ko_prol_VS___wt_prol)

# Sign change, proliferation as reference
koDvskoP <- koDvskoP %>% mutate(logFC__ko_diff_VS___ko_prol = -logFC__ko_diff_VS___ko_prol)
wtDvswtP <- wtDvswtP %>% mutate(logFC__wt_diff_VS___wt_prol = -logFC__wt_diff_VS___wt_prol)




#DEGS CALCULTOR FUNCTION-------------------------------------------------------------------------------------------------------

degs_calculator <- function(data, padj_threshold= 0.05, log2FC_threshold =0.65, print= TRUE){
  
  # this function calculate the up and down DEGS from DEG analysis data
  # data - DEG analysis data (reanalyzerGSE output) 
  # padj_threshold & log2FC_threshold
  # output - list with up and down DEGs
  
  up_degs <- data %>% filter(FDR < padj_threshold, data[[3]]>log2FC_threshold) %>%
    mutate(Direction = "UP")
  down_degs <- data %>% filter(FDR < padj_threshold, data[[3]] < -log2FC_threshold)%>%
    mutate(Direction = "DOWN")
  
  degs <- 
    if (print){
      print(paste(nrow(up_degs), ' UP and ',nrow(down_degs), ' DOWN', ' have been found.', ))
    }
  
  return(list(UP = up_degs, DOWN = down_degs)) 
}

# VOLCANO PLOT FUNCION ------------------------------------------------------------------------------------------------

volcano_plot <- function(data, title = 'Volcano', top_up= 5, top_down= 5, padj_threshold= 0.05, log2FC_threshold =0.65){
  
  # this function plot the volcano's from DEG analysis data
  # data - DEG analysis data (reanalyzerGSE output) 
  # top_up & top_down - number of genes to be write in the plot
  # padj_threshold & log2FC_threshold for degs_calculator function
  # output - volcano plot
  
  logFC = colnames(data)[3] #log name normalization, change in each file
  
  # DEGs calculation
  degs <- degs_calculator(data, padj_threshold, log2FC_threshold, print=FALSE)
  top10_up <- degs$UP %>% dplyr::arrange(desc(.[[3]]))%>% head(top_up) 
  top10_down <- degs$DOWN %>% dplyr::arrange(.[[3]])%>% head(top_down) 
  
  topfdr_up <- degs$UP %>% dplyr::arrange(.[[6]])%>% head(top_up)
  topfdr_down <- degs$DOWN %>% dplyr::arrange(.[[6]])%>% head(top_down)
  
  top_up <- bind_rows(top10_up, topfdr_up) %>% distinct()
  top_down <- bind_rows(top10_down, topfdr_down) %>% distinct()
  
  # Plot generation
  ggplot(data, aes_string(x= logFC, y= "-log10(FDR)"))+ #aes_string, column info is a string
    
    geom_point(aes(color = ifelse(FDR < padj_threshold & .data[[logFC]] > log2FC_threshold , "red",   #change point colour if is down, up or no deg
                                  ifelse(FDR < padj_threshold & .data[[logFC]] < -log2FC_threshold, "blue", 'grey'))))+
    theme_bw() +
    
    labs(
      title = title,
      subtitle = paste("Up: ", nrow(degs$UP)," Down: ", nrow(degs$DOWN), "Total: ", nrow(degs$UP)+ nrow(degs$DOWN)),
      x = "log2(Fold-Change)",
      y = "-log10(FDR)",
      color = "FDR<0.05 & LogFC>1") +
    
    # Top upregulated genes annotation 
    geom_text_repel(data = top_up,
                    aes(label = Gene_ID),
                    size = 3, max.overlaps = Inf,
                    box.padding = 0.5) +
    
    
    # Top downregulated genes annotation 
    geom_text_repel(data = top_down,
                    aes(label = Gene_ID),
                    size = 3, max.overlaps = Inf,
                    box.padding = 0.5)+
    
    #geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Lineas verticales en log2FC de -1 y 1
    #geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # Linea horizontal en p-adj 0.05
    
    scale_color_manual(values = c("dodgerblue", "grey", "firebrick2"), labels = c('DOWN', '', 'UP'))+
    theme(legend.position = 'none',
          plot.title = element_text(face = 'bold', size = 15))
}


# Volcano PLOTs

# KOvsWT volcanos, with Cort gene annotation
p1 <- volcano_plot(koDvswtD, title = 'Cort-/- vs Cort+/+ (Differentiation)') + geom_text_repel(data = filter(koDvswtD, Gene_ID=='Cort'),
                                                                                               aes(label = 'Cort'),
                                                                                               size = 3, max.overlaps = Inf,
                                                                                               box.padding = 0.5) 

p2 <- volcano_plot(koPvswtP, title = 'Cort-/- vs Cort+/+ (Proliferation)')+ geom_text_repel(data = filter(koPvswtP, Gene_ID=='Cort'),
                                                                                            aes(label = 'Cort'),
                                                                                            size = 3, max.overlaps = Inf,
                                                                                            box.padding = 0.5)
p1

# DiffvsProliferation
p3 <- volcano_plot(koDvskoP, title = 'Cort-/- Dif vs Cort-/- Prolif', log2FC_threshold =1)
p4 <- volcano_plot(wtDvswtP, title = 'Cort+/+ Dif vs Cort+/+ Prolif', log2FC_threshold =1)

#save
ggsave("plots/Volcano/KoDvsWtD.png", plot = p1, width = 5, height = 5, dpi = 300)
ggsave("plots/Volcano/KoPvsWtP.png", plot = p2, width = 5, height = 5, dpi = 300)
ggsave("plots/Volcano/KoDvsKoP.png", plot = p3, width = 5, height = 5, dpi = 300)
ggsave("plots/Volcano/WtDvsWtP.png", plot = p4, width = 5, height = 5, dpi = 300)


# CALCULO DE CPM----------------------------------------------------------------------------------
# reanalyzerGSE doesnt calculate CPM, 

gene_size <- read.table('ReadCount/str-Size.tab', header = T, row.names = 1)
gene_count <- read.table('ReadCount/str-ReadCount.tab', header = T, row.names = 1)
names(gene_count) <- gsub('_1.fastq.gz_nat', '', names(gene_count))


gene_size1 <- gene_size[rownames(gene_count),, drop=F] #ordenar los tamaños para que esten en el mismo orden que gene count

rpk <- gene_count/gene_size1$Length
rpk_sum <- colSums(rpk)

tpm <- (rpk/rpk_sum)*1e6

write_tsv(tpm, file = 'ReadCount/tpm.tsv')
write.table(tpm, file = 'ReadCount/tpm.txt', sep = '\t')

universe <- tpm %>% filter(!if_all(everything(), ~ . == 0)) #quitar todos los genes donde no se detecta expresion
universe <- rownames(universe)
head(universe)
write_lines(universe, 'ReadCount/universe.txt')


#FUNCION GO_Calculator------------------------------------------------------------

plots_dir_go <- 'plots/GO'
if (!dir.exists(plots_dir_go)) {
  dir.create(plots_dir_go)}


go_calculator <- function(degs, filename, list=FALSE, result=F){

  # degs - output data from degs_calculator, with UP and DOWN
  # filename - name of the file to generate
  #list -  boolean which allows to use simple list of degs (no degs_calculator output)

  if (list){
    print(paste('Analizando: ', filename))
    BP_degs <- enrichGO(gene = degs,
                        keyType = 'SYMBOL',
                        OrgDb         = org.Mm.eg.db,
                        ont           = "BP",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        universe = universe)
    print(paste(nrow(BP_degs), ' BP GO'))
    
    CC_degs<- enrichGO(gene = degs,
                       keyType = 'SYMBOL',
                       OrgDb         = org.Mm.eg.db,
                       ont           = "CC",
                       pAdjustMethod = "BH",
                       pvalueCutoff  = 0.01,
                       qvalueCutoff  = 0.05,
                       readable      = TRUE,
                       universe = universe)
    print(paste(nrow(CC_degs), ' CC GO'))
    MF_degs <- enrichGO(gene = degs,
                        keyType = 'SYMBOL',
                        OrgDb         = org.Mm.eg.db,
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        pvalueCutoff  = 0.01,
                        qvalueCutoff  = 0.05,
                        readable      = TRUE,
                        universe = universe)
    print(paste(nrow(MF_degs), ' MF GO'))
    
    
    if (result){
      return(list(BP =BP_degs, CC=CC_degs, MF = MF_degs))  
    }
    #sSAVE
    results_dir_go<- 'results/GO'
    if (!dir.exists(results_dir_go)) {
      dir.create(results_dir_go)
    }
    write_tsv(BP_degs@result, file.path(results_dir_go,paste('BP_', filename, "_GO.tsv", sep='')))
    write_tsv(CC_degs@result, file.path(results_dir_go,paste('CC_', filename, "_GO.tsv", sep='')))
    write_tsv(MF_degs@result, file.path(results_dir_go,paste('MF_', filename, "_GO.tsv", sep='')))
    
    go_plot(BP_degs, filename = filename, type = 'BP')
    go_plot(MF_degs, filename = filename, type = 'MF')
    go_plot(CC_degs, filename = filename, type = 'CC')
    
    
  } 
  
  else{#gene list
    print(paste('Analizando: ', filename))
    
    for (i in 1:length(degs)){
      if (i==1){ #UP---------------------------------------------------------------------
        BP_degs_UP <- enrichGO(gene = degs$UP$Gene_ID,
                               keyType = 'SYMBOL',
                               OrgDb         = org.Mm.eg.db,
                               ont           = "BP",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               universe = universe)
        print(paste(nrow(BP_degs_UP), 'BP GO UP'))
        
        CC_degs_UP <- enrichGO(gene = degs$UP$Gene_ID,
                               keyType = 'SYMBOL',
                               OrgDb         = org.Mm.eg.db,
                               ont           = "CC",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               universe = universe)
        print(paste(nrow(CC_degs_UP), 'CC GO UP'))
        MF_degs_UP <- enrichGO(gene = degs$UP$Gene_ID,
                               keyType = 'SYMBOL',
                               OrgDb         = org.Mm.eg.db,
                               ont           = "MF",
                               pAdjustMethod = "BH",
                               pvalueCutoff  = 0.01,
                               qvalueCutoff  = 0.05,
                               readable      = TRUE,
                               universe = universe)
        print(paste(nrow(MF_degs_UP), 'MF GO UP'))
        
        #save
        results_dir_go_up <- 'results/GO/UP'
        if (!dir.exists(results_dir_go_up)) {
          dir.create(results_dir_go_up)
        }
        write_tsv(BP_degs_UP@result, file.path(results_dir_go_up,paste('BP_', filename, "_GO_UP.tsv", sep='')))
        write_tsv(CC_degs_UP@result, file.path(results_dir_go_up,paste('CC_', filename, "_GO_UP.tsv", sep='')))
        write_tsv(MF_degs_UP@result, file.path(results_dir_go_up,paste('MF_', filename, "_GO_UP.tsv", sep='')))
        
      }else if(i==2){ #DOWN-----------------------------------------------------------------
        BP_degs_DOWN <- enrichGO(gene = degs$DOWN$Gene_ID,
                                 keyType = 'SYMBOL',
                                 OrgDb         = org.Mm.eg.db,
                                 ont           = "BP",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE,
                                 universe = universe)
        print(paste(nrow(BP_degs_DOWN), 'BP GO DOWN'))
        
        CC_degs_DOWN <- enrichGO(gene = degs$DOWN$Gene_ID,
                                 keyType = 'SYMBOL',
                                 OrgDb         = org.Mm.eg.db,
                                 ont           = "CC",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE,
                                 universe = universe)
        print(paste(nrow(CC_degs_DOWN), 'CC GO DOWN'))
        
        MF_degs_DOWN <- enrichGO(gene = degs$DOWN$Gene_ID,
                                 keyType = 'SYMBOL',
                                 OrgDb         = org.Mm.eg.db,
                                 ont           = "MF",
                                 pAdjustMethod = "BH",
                                 pvalueCutoff  = 0.01,
                                 qvalueCutoff  = 0.05,
                                 readable      = TRUE,
                                 universe = universe)
        print(paste(nrow(MF_degs_DOWN), 'MF GO DOWN'))
        
        #sAVE
        results_dir_go_DOWN <- 'results/GO/DOWN'
        if (!dir.exists(results_dir_go_DOWN)) {
          dir.create(results_dir_go_DOWN)
        }
        write_tsv(BP_degs_DOWN@result, file.path(results_dir_go_DOWN,paste('BP_', filename, "_GO_DOWN.tsv", sep='')))
        write_tsv(CC_degs_DOWN@result, file.path(results_dir_go_DOWN,paste('CC_', filename, "_GO_DOWN.tsv", sep='')))
        write_tsv(MF_degs_DOWN@result, file.path(results_dir_go_DOWN,paste('MF_', filename, "_GO_DOWN.tsv", sep='')))
      }
    }
    
  } 
  
}

go_plot <- function(go_result, filename, type, sign = ''){ 
  # go_result - output from go_calculation function, use the enrichplot function to generate fast graphs
  
  # filaname, type and sign are just for create the output file name
  # filename - name of the comparison
  # type - BP, MF or CC
  # sign - UP or DOWN. both as default 
  
  
  plots_dir_go_unic <- paste('plots/GO/', filename, sep='')
  if (!dir.exists(plots_dir_go_unic)) {
    dir.create(plots_dir_go_unic)}
  
  print(paste(filename, type, '.pdf', sep=''))
  
  pdf(file.path(plots_dir_go_unic,paste(type, '_',filename, '_',sign, '.pdf', sep ='')), width = 12, height = 6)
  
  # all the graphs are error controlled, beacause is easy to fail
  tryCatch({
    bar <- barplot(go_result, showCategory = 10)+
      ggtitle(paste(sign, 'DEG'))+
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(bar)
  }, error =function(e){message("Error al crear el treeplot: ", e$message)})
  
  tryCatch({ 
    enriches <- pairwise_termsim(go_result)
    print(
      treeplot(enriches) +
        ggtitle(paste(type, '- Treeplot', filename, type)) + labs(color='NES')+
        theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))+
        scale_color_gradient(low = "blue", high = "red", limits = c(-2, 2),  oob = scales::squish)
    )}, error = function(e) {message("Error al crear el treeplot: ", e$message)})
  
  tryCatch({
    emap <- emapplot(
      enriches, 
      showCategory = 40,           
      node_label = "group", 
      shadowtext = TRUE,
      layout = 'nicely',
      label_size = 1
    ) + 
      ggtitle(paste(type, "- Enrichment Map", type)) + 
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    
    print(emap)
    
  }, error = function(e) {message("Error al crear el emapplot: ", e$message)})
  
  tryCatch({
    options(enrichplot.colours = c('yellow', 'violetred4'))
    print(dotplot(go_result,showCategory=10, split=".sign", color = 'p.adjust') +  
            facet_grid(.~.sign))
    
  }, error = function(e) {message("Error al crear el dotplot: ", e$message)})
  
  dev.off()
}

#CALCULO DE DEGS---------------------------------------------------------------------------

degs_koDvswtD <- degs_calculator(koDvswtD)
degs_koPvswtP <- degs_calculator(koPvswtP)
degs_koDvskoP <- degs_calculator(koDvskoP, log2FC_threshold =1)
degs_wtDvswtP <- degs_calculator(wtDvswtP, log2FC_threshold =1)


max(degs_koDvswtD$UP$logFC__ko_diff_VS___wt_diff)
min(degs_koDvswtD$DOWN$logFC__ko_diff_VS___wt_diff)

# write de degs of KOvsWT, splitting UP and DOWN
writeLines(degs_koDvswtD$UP$Gene_ID, file.path(results_dir, 'UP_Diff.txt'))
writeLines(degs_koDvswtD$DOWN$Gene_ID, file.path(results_dir, 'DOWN_Diff.txt'))
writeLines(degs_koPvswtP$UP$Gene_ID, file.path(results_dir, 'UP_Prol.txt'))
writeLines(degs_koPvswtP$DOWN$Gene_ID, file.path(results_dir, 'DOWN_Prol.txt'))


# overlap genes calculation
diffprocess_UP <- calculate.overlap(list(degs_koDvskoP$UP$Gene_ID,degs_wtDvswtP$UP$Gene_ID ))
diffprocess_DOWN <- calculate.overlap(list(degs_koDvskoP$DOWN$Gene_ID,degs_wtDvswtP$DOWN$Gene_ID ))

WTvsKO_UP <- calculate.overlap(list(degs_koDvswtD$UP$Gene_ID,degs_koPvswtP$UP$Gene_ID ))
WTvsKO_DOWN <- calculate.overlap(list(degs_koDvswtD$DOWN$Gene_ID,degs_koPvswtP$DOWN$Gene_ID ))

# Export overlaps groups - a1 onlyKO, a2 onlyWT, a3 all
results_dir_overlapsdiff <- 'results/Overlaps/Differentiation_Process'
if (!dir.exists(results_dir_overlapsdiff)) {
  dir.create(results_dir_overlapsdiff)
}

writeLines(diffprocess_UP$a1, file.path(results_dir_overlapsdiff, 'UP_Diff_onlyKO.txt'))
writeLines(diffprocess_UP$a2, file.path(results_dir_overlapsdiff, 'UP_Diff_onlyWT.txt'))
writeLines(diffprocess_UP$a3, file.path(results_dir_overlapsdiff, 'UP_Diff_Common.txt'))

writeLines(diffprocess_DOWN$a1, file.path(results_dir_overlapsdiff, 'DOWN_Diff_onlyKO.txt'))
writeLines(diffprocess_DOWN$a2, file.path(results_dir_overlapsdiff, 'DOWN_Diff_onlyWT.txt'))
writeLines(diffprocess_DOWN$a3, file.path(results_dir_overlapsdiff, 'DOWN_Diff_Common.txt'))

results_dir_overlapsKOvsWT <- 'results/Overlaps/KOvsWT'
if (!dir.exists(results_dir_overlapsKOvsWT)) {
  dir.create(results_dir_overlapsKOvsWT)
}

writeLines(WTvsKO_UP$a1, file.path(results_dir_overlapsKOvsWT, 'UP_KOvsWT_onlyKO.txt'))
writeLines(WTvsKO_UP$a2, file.path(results_dir_overlapsKOvsWT, 'UP_KOvsWT_onlyWT.txt'))
writeLines(WTvsKO_UP$a3, file.path(results_dir_overlapsKOvsWT, 'UP_KOvsWT_Common.txt'))

writeLines(WTvsKO_DOWN$a1, file.path(results_dir_overlapsKOvsWT, 'DOWN_KOvsWT_onlyKO.txt'))
writeLines(WTvsKO_DOWN$a2, file.path(results_dir_overlapsKOvsWT, 'DOWN_KOvsWT_onlyWT.txt'))
writeLines(WTvsKO_DOWN$a3, file.path(results_dir_overlapsKOvsWT, 'DOWN_KOvsWT_Common.txt'))



#Venn Diagrmams----------------------------------------------------------------------------------------------------------------------------------------------------

# Common and differents genes in differentiation. How KO and WT changes
Venn_ProcessDiff_UP <- ggVennDiagram(list(degs_koDvskoP$UP$Gene_ID,degs_wtDvswtP$UP$Gene_ID ), 
                                     category.names = c("Cort-/-","WT"), set_size = 6,edge_size = 1.2, label_size = 8)+
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  scale_fill_gradient(low = '#ffcccc', high = 'firebrick2') +
  theme(legend.position = 'none')

Venn_ProcessDiff_DOWN <- ggVennDiagram(list(degs_koDvskoP$DOWN$Gene_ID,degs_wtDvswtP$DOWN$Gene_ID ), 
                                       category.names = c("Cort-/-","Cort+/+"), set_size = 6, label_size = 8)+
  scale_x_continuous(expand = expansion(mult = 0.1)) + 
  scale_fill_gradient(low = '#d2e8ff', high = 'dodgerblue') +
  theme(legend.position = 'none')

combined_diff <- (Venn_ProcessDiff_UP + Venn_ProcessDiff_DOWN) + 
  plot_annotation(title = 'DEG during Differentiation', 
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))

combined_diff


ggsave('plots/Venn_ProcessDiff.png', plot = combined_diff, device = 'png', width = 10) #export


# DEGS KOvsWT
Venn_KOVSWT_UP <- ggVennDiagram(list(degs_koDvswtD$UP$Gene_ID,degs_koPvswtP$UP$Gene_ID ), 
                                category.names = c("Diff","Prol"), set_size = 6,edge_size = 1.2, label_size = 8)+
  scale_x_continuous(expand = expansion(mult = .1)) + 
  scale_fill_gradient(low = '#ffcccc', high = 'firebrick2') +
  theme(legend.position = 'none')

Venn_KOVSWT_DOWN <- ggVennDiagram(list(degs_koDvswtD$DOWN$Gene_ID,degs_koPvswtP$DOWN$Gene_ID ), 
                                  category.names = c("Diff","Prol"), set_size = 6,edge_size = 1.2, label_size = 8)+
  scale_x_continuous(expand = expansion(mult = .1)) + 
  scale_fill_gradient(low = '#d2e8ff', high = 'dodgerblue') +
  theme(legend.position = 'none')

combined_kovswt <- (Venn_KOVSWT_UP + Venn_KOVSWT_DOWN) + 
  plot_annotation(title = 'DEG KOvsWT', 
                  theme = theme(plot.title = element_text(size = 16, face = "bold")))
combined_kovswt

ggsave('plots/Venn_KOVSWT.png', plot = combined_kovswt, device = 'png', width = 10) #export


# Functional enrichment GOs---------------------------------------------------------------------------------------------

#Fusion Up Y DOWN DEGs
degs_koDvswtDC<- bind_rows(degs_koDvswtD$UP %>% mutate(Regulation = "UP"),
                           degs_koDvswtD$DOWN %>% mutate(Regulation = "DOWN"))
degs_koPvswtPC<- bind_rows(degs_koPvswtP$UP %>% mutate(Regulation = "UP"),
                           degs_koPvswtP$DOWN %>% mutate(Regulation = "DOWN"))

degs_koDvskoPC<- bind_rows(degs_koDvskoP$UP %>% mutate(Regulation = "UP"),
                           degs_koDvskoP$DOWN %>% mutate(Regulation = "DOWN"))

degs_wtDvswtPC<- bind_rows(degs_wtDvswtP$UP %>% mutate(Regulation = "UP"),
                           degs_wtDvswtP$DOWN %>% mutate(Regulation = "DOWN"))


degs_listC <- list(koDvswtD = degs_koDvswtDC, koPvswtP=degs_koPvswtPC, 
                   koDvskoP= degs_koDvskoPC, wtDvswtP= degs_wtDvswtPC)

overlaps_list <- list(diffprocess_UP_onlyKO =diffprocess_UP$a1, diffprocess_UP_onlyWT =diffprocess_UP$a2,
                      diffprocess_UP_Common =diffprocess_UP$a3,
                      diffprocess_DOWN_onlyKO =diffprocess_DOWN$a1, diffprocess_DOWN_onlyWT =diffprocess_DOWN$a2,
                      diffprocess_DOWN_Common =diffprocess_DOWN$a3,
                      WTvsKO_UP_common = WTvsKO_UP$a3, WTvsKO_DOWN_common = WTvsKO_DOWN$a3)


#Functional enrichment of WGCNA interest modules (WGCNA)
cyan_genes <- as.vector(t(read.table("WGCNA/Modules/cyan_genes.txt")))
darkolivegreen_genes <- as.vector(t(read.table("WGCNA/Modules/darkolivegreen_genes.txt")))
darkturquoise_genes <- as.vector(t(read.table("WGCNA/Modules/darkturquoise_genes.txt")))
green_genes <- as.vector(t(read.table("WGCNA/Modules/green_genes.txt")))
lightsteelblue_genes <- as.vector(t(read.table("WGCNA/Modules/lightsteelblue1_genes.txt")))

modules <- list(cyan= cyan_genes, darkolivegreen= darkolivegreen_genes,darkturquoise= darkturquoise_genes,
                green= green_genes, lightsteelblue= lightsteelblue_genes)


result_list <- list()

# DEGs GO

for (degs_name in names(degs_listC)){
  #  degs data must be a file from degs_calculator, with up and down splitted
  degs <- degs_listC[[degs_name]]$Gene_ID
  result_list[[degs_name]] <- go_calculator(degs,degs_name, list=T, result=T)
}


#FUNCION PARA CALCULAR PROPORCION DE UP Y DOWNS------------------------------------------------------------------------------------------

proportionUP_cal <- function(data, degs){
  
  # this funciton calculate the proportion of up and down (up/total) degs of GO terms
  # data - GO output file
  # degs - degs_calculation output 
  
  # Create new columns in the GO output file with nuymber of UP, DOWN DEGs and the proportion
  data@result$numUP <- sapply(data@result$geneID, function(gene_list) {
    genes <- unlist(strsplit(gene_list, "/")) # split the gene list of the file
    sum(genes %in% degs$UP$Gene_ID)  # Count genes UP
  })
  
  data@result$numDOWN <- sapply(data@result$geneID, function(gene_list) {
    genes <- unlist(strsplit(gene_list, "/")) # split the gene list of the file
    sum(genes %in% degs$DOWN$Gene_ID)  # Count genes DOWN
  })
  
  data@result$proportionUP <- with(data@result, {
    numUP / (numUP + numDOWN)  # UP degs proportion calculation
  })
  
  return(data) #GO output data modified
}


# GO output data isolation
GO_koDvswtD <- result_list$koDvswtD$BP 
GO_koPvswtP <- result_list$koPvswtP$BP

# -log10(padjust) calculation
GO_koDvswtD@result$log10_p.adjust <- -log10(GO_koDvswtD@result$p.adjust)
GO_koPvswtP@result$log10_p.adjust <- -log10(GO_koPvswtP@result$p.adjust)

# proportionUP uses
GO_koDvswtD <- proportionUP_cal(GO_koDvswtD,degs_koDvswtD)
GO_koPvswtP <- proportionUP_cal(GO_koPvswtP,degs_koPvswtP)


#SAVE BP GO files 
results_dir_go <- 'results/GO'
if (!dir.exists(results_dir_go)) {
  dir.create(results_dir_go)
}

write_tsv(GO_koDvswtD@result, file.path(results_dir_go,'BP_koDvswtD_GO.tsv'))
write_tsv(GO_koPvswtP@result, file.path(results_dir_go,'BP_koPvswtP_GO.tsv'))


dotbarplot <- function(data, title= 'Plot Title', showCategory=10, down = F){
  # make a plot combining a dotplot of the log10(pval) and the barplot of nu,ber of genes
  # data - GO with proportionUp_cal and -log10p.adjust previously calculated
  # showCategory - number of GO to be plotted 
  # down - boolean to plot the GOs with the low proportionUP (mainly down DEGs)
  
  
  # order by significance and if down is T select the Go terms with more down DEGs
  data_sort <- data@result %>%
    {if (down) filter(., proportionUP <0.5) else .}  %>% arrange(desc(log10_p.adjust)) 
  
  #dotplot
  dot <- ggplot(head(data_sort, showCategory), aes(x=reorder(str_wrap(Description, width = 30), log10_p.adjust), y= log10_p.adjust))+
    geom_point(aes(size = Count, color = log10_p.adjust))+
    scale_color_gradient(low = "pink", high = "red")+
    scale_size(range = c(3, 10))+
    labs(x='', y='-log10(p.adjust)', size='Gene Count', color ='-log10(p.adjust)')+
    theme_minimal()+
    theme(axis.text.y = element_text(size = 14, color = 'black'),
          axis.text.x = element_text(color = 'black'),
          panel.grid.major.x = element_line(size = 1),
          panel.grid.minor.x = element_blank())+
    coord_flip()
  
  #barplot
  bar <- ggplot(head(data_sort, showCategory), aes(x=reorder(str_wrap(Description, width = 10) , log10_p.adjust), y =log10_p.adjust))+
    geom_bar(aes(y=numUP), fill = 'firebrick2', stat = 'identity', position = 'stack')+
    geom_bar(aes(y = numDOWN), fill = "dodgerblue", stat = "identity", position = "stack") +
    scale_y_continuous(labels = abs) +  # absolute values in Y axis
    labs(x = NULL, y = "Number of Genes", fill = "Expression") +
    theme_minimal() +
    theme(axis.text.y = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_line(size = 1, linetype = 2, color = 'darkgrey'),
          axis.text.x = element_text(color = 'black'),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.ontop = TRUE)+
    coord_flip()
  
  combined_plot <- dot + bar + plot_layout(ncol = 2, widths = c(1, 1.5), guides = 'collect') + 
    plot_annotation(title = title, theme = theme(plot.title = element_text(face='bold', size=16)))
  
  return(combined_plot)
}


DB_dif <- dotbarplot(GO_koDvswtD, title = 'Cort-/- vs Cort+/+ (Differentation)', showCategory = 15)
DB_dif

DB_dif_DOWN <- dotbarplot(GO_koDvswtD, title = 'Cort-/- vs Cort+/+ (Differentation) MAIN DOWN', showCategory = 10, down = T)
DB_dif_DOWN

DB_prol <-dotbarplot(GO_koPvswtP, title = 'Cort-/- vs Cort+/+ (Proliferation)', showCategory = 15)
DB_prol

ggsave('plots/GO/dotbar_dif.png', plot = DB_dif, height = 9, width = 12)
ggsave('plots/GO/dotbar_prol.png', plot = DB_prol, height = 9, width = 12)
ggsave('plots/GO/dotbar_diffdown.png', plot = DB_dif_DOWN, height = 9, width = 12)

combined_plot


# Overlaps GOs
for (over_name in names(overlaps_list)){
  # degs data must be a file from degs_calculator, with up and down splitted
  over <- overlaps_list[[over_name]]
  go_calculator(over,over_name, list=TRUE)
}


# WGCNA modules GO calculation
for (module_name in names(modules)){
  #  degs data must be a file from degs_calculator, with up and down splitted

  module <- modules[[module_name]]
  go_calculator(module,module_name, list=TRUE)
}


# CLUSTER  STRING-----------------------------------------------------------------------------------------------------------------------------------------

cluster_1 <- c("Serpinb9b", "Serpinb9", "Acta2", "Pdgfb", "Shh", "Tgfb3", "Igf1r", "Ccl4", "Ogn", "Loxl2", 
               "Pdgfrb", "Col1a2", "Bgn", "Col4a1", "Col4a2", "Mmp2", "Thbs1", "Cxcl10", "Ltbp2", "Col8a1", 
               "Fbln1", "Fbln2", "Klf4", "Hspg2", "Igf1", "Bcl2", "Spp1", "Tgfbr2", "Fbn1", "Itgb1", "Postn", 
               "Fn1", "Nid1", "Il7", "Tlr3", "Bdnf", "Fcgr4", "Aurkb", "Hells", "Spc25", "Kif20b", "Kifc1", 
               "Prc1", "Kif18b", "Spag5", "Smc4", "Ect2", "Cep55", "Kif20a", "Kif23", "Cenpf", "Cenpe", 
               "Trem2", "Stat6", "Casp8", "Cx3cr1", "Irf8", "Fcer1g", "Tlr7", "Col8a2", "Plod2")

cluster_2 <- c("Ccl2", "Cav1", "Eng", "Cybb", "Sparc", "Timp3", "Fbln5", "Edn1", "Col5a1", "Mki67", 
               "Ptgs2", "Serpine1", "Tlr4", "Fstl1", "Prom1", "Loxl1", "Ccnd1", "Cxcr4", "Itgam", 
               "Icam1", "Nes", "Kdr", "Ptprc")

cluster_3 <- c( "Col1a1", "Vtn", "Col2a1", "Itga5", "Lox", "Twist1", "Igf2", "Csf1r", "Col5a2", "Itgb5", 
                "Itga3", "Col5a3", "Adamts4", "Hbegf", "Cd33", "Ets1", "Erbb3", "Serpinf1", "Fkbp10", 
                "Copz2", "Nid2", "Lama2", "Col12a1", "Cd48", "Lama4", "Col4a4", "Thbs3", "Emilin1")

cluster_4 <- c("Cldn11", "Tgfb1i1", "Mbp", "Ugt8a", "Cnp", "Sox10", "Mag", "Nfatc4", "Vgll3", "Gjc2", 
               "Myrf", "Mobp")
cluster_5 <- c( "Cnn1", "Tpm4", "Actn1", "Pdgfd", "Tagln2", "Cald1", "Tpm1", "Lmod1", "Tpm2", "Calml4", 
                "Epha2", "Nrp1", "Ephb4")
cluster_6 <- c("Runx1", "Cdk6", "Plk3", "Brca2", "Fzd7", "Fzd6", "Fzd5", "Ccnd2", "Tcf7", "Wnt9a", 
               "Wnt3", "Sfrp4", "Sfrp5")

string_cluster <- list('cluster_1'=cluster_1, 'cluster_2'=cluster_2, 'cluster_3'=cluster_3, 'cluster_4'=cluster_4, 'cluster_5'=cluster_5, 'cluster_6'=cluster_6)

result_list_gostring <- list()

# Go calculation of differents clusters 
result_list_gostring[['cluster_1']]<- go_calculator(cluster_1, 'cluster_1', list=T, result = T)
result_list_gostring[['cluster_2']]<- go_calculator(cluster_2, 'cluster_2', list=T, result = T)
result_list_gostring[['cluster_3']]<- go_calculator(cluster_3, 'cluster_3', list=T, result = T)
result_list_gostring[['cluster_4']]<- go_calculator(cluster_4, 'cluster_4', list=T, result = T)
result_list_gostring[['cluster_5']]<- go_calculator(cluster_5, 'cluster_5', list=T, result = T)
result_list_gostring[['cluster_6']]<- go_calculator(cluster_6, 'cluster_6', list=T, result = T)


for (result in names(result_list_gostring)){
  export <- result_list_gostring[[result]]$BP@result %>% head(20) %>% dplyr::select(Description, p.adjust)
  write.table(export, file=paste('STRING/BP', result, '.txt', sep = ''), quote = F, sep='\t', row.names = F)
  
}


