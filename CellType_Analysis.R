# 2024 - Pablo Vargas Rodríguez
# RNA seq analysis of oligodendrocytes WT and KO for Cort-/-
# Cell Type Analysis - GSEA and plot of tpm plots


library(readr)
library(ggplot2)
library(VennDiagram) #calculosVenn
library(ggVennDiagram) #Venn Bonitos
library(ggrepel)
library(magrittr) #usar los %>%
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork) #juntar graficas

library(enrichplot)
library(clusterProfiler)
library(msigdbr)# Package that contains MSigDB gene sets in tidy format
library(org.Mm.eg.db)


#Overlap with cell type marker--------------------------------------------------------------------------------------------------------------
# https://www.jneurosci.org/content/34/36/11929/tab-figures-data
# astro - astrocyte, OPC - oligodendrocyte progenitor cell, NFO - Newly formed oligodendrocyte, MO- Myelinating oligodendrocyte
# micro - microglia, , endo - endothelial, peri - pericyte

astro <- c('Hgf', 'Aqp4', 'Itih3', 'Bmpr1b', 'Itga7', 'Plcd4', 'Grm3', 'Slc14a1', 'Phkg1',
           'Pla2g3', 'Cbs', 'Paqr6', 'Aldh1l1', 'Cth', 'Ccdc80', 'Fmo1', 'Slc30a10', 'Slc6a11',
           'Fgfr3', 'Slc4a4','Gdpd2', 'Ppp1r3c', 'Grhl1', 'Entpd2', 'Egfr', 'Al464131', 'Otx1',
           'Nwd1', 'Atp13a4', 'Kcnn3', 'Ptx3', 'Sorc2', 'Tnc', 'Sox9', 'Abcd2', 'Fzd10', 'Lrig1',
           'Mlc1', 'Chrdl1', 'Aifm3',
           'Gli1', 'Gli2', 'Gli3', 'Otx1','Grhl1', 'Sox9', 'Hes5', 'Rfx4', 'Pax6', 'Dbx2')
astro <- unique(astro)


OPC <- c('Pdgfr1', 'Lnx1', 'Dcn', 'Mmp15', 'Cdo1', 'Sapcd2', 'Kcnk1', 'Rasgrf1', ' Pcdh15', 'Chrna4',
         'Dll3', 'Col1a2', 'Fam70b', 'Sstr1', 'Pnlip', 'Cspg4', 'Lppr1', 'Ppapdc1a', 'Nxph1', 'Pid1',
         'Ugdh', 'Slitrk1', 'Shc4', 'Smoc1', 'Emid1', 'Rlbp1', 'Dcaf12l1', 'Lypd6', 'Lhfpl3', 'Myt1',
         'Gfia3', 'C1ql1', 'Tmem179', 'Megf11', 'Ncald', 'Sdc3', 'Rprm', 'Cacng4', 'Grin3a', 'Fam5c',
         'Sox10', 'Gsx1', 'Olig1', 'Myt1', 'Pou3f1', 'Sox8', 'Olig2', 'Sox3', 'Nkx2.2', 'Sox6')
OPC <- unique(OPC)

NFO <- c('Gp1bb', 'Tmem108', 'Fyn', 'Ust', 'Mical3', 'Kif19a', '1810041L15Rik', '9630013A20Rik',
         'Nfasc', 'Ssh3', 'Pik3r3', 'Enpp6', 'Tns3', 'Bmp4', 'Mcl1', 'Cdv3', 'Tmem163', 'Rap2a',
         'Tmem2', 'Cnksr3', 'Cyfip2', 'Fmd4a', 'Slc12a2', ' Itpr2', 'Rnf122', 'Lims2', 'Samd4b',
         'Chn2', 'Pp2r3a', 'Strn', 'Glrb', 'Rras2', 'Fmnl2', 'Sema5a', 'Fam3c', 'Cdc37l1', 'Fam73a',
         'Elovl6', 'Atrn', 'Lrrc42', 
         'Myrf', 'Nkx6.2', 'Sox10', 'Barx2', 'Olig1', 'Nkx2.2', 'Sox8', 'Olig2', 'Sox3', 'Mycl1')
NFO <- unique(NFO)

MO <- c('Gjb1', 'Ndrg1', 'Ppp1r14a', 'Adssl1', 'Aspa', 'Acy3', 'Trp53inp2', 'Pla2g16', 'Efhd1', 
        'Itgb4', 'Hapln2', 'Mbp', 'Hcn2', 'Nmra1', 'Cdc42ep2', 'Mal', 'Mog', 'Slco3a1', 'Apod',
        'Gsn', 'Pdlim2', 'Prr18', 'Inf2', 'Tppp3', 'Tbc1d9b', 'Nol3', 'Cenpb', 'Slc45a3', 'Carns1',
        'Opalin', 'Arsg', 'Rftn1', 'Adap1', 'Plekhb1', 'Trf', 'Insc', 'Cryab', 'Kif5a', 'Trak2', 'Cldn11',
        'Nkx6.2', 'Myrf', 'Sox10', 'Sp7', 'Barx2', 'Olig1', 'Pou3f1', 'Sox8', 'Carhsp1', 'Nfe2l3')
MO <- unique(MO)

micro <- c('Slfn2', 'Gpr84', 'Ccr7', 'Bcl2a1d', 'Tnf', 'Ncf1', 'Gdf15', 'Osm', 'Lrrc25', 
           'H2-Oa', 'Cd83', 'Ccl3', 'Slamf8', 'Ccl4', 'Gna15', 'Il1b', 'Plau', 'Ccl9',
           'Tmem119', 'C1qa', 'Irf8', '1810011H11Rik', 'Pla2g15', 'Cxcl16', 'Ch25h', 'Hck',
           'Ccl12', 'Ptafr', 'Cd300a', 'Irf5', 'Sfpi1', 'Selplg', 'Sash3', 'Pltp', 'Trem2',
           'Tlr2', 'P2ry6', 'Cdl4', 'Bcl2a1a', 'Bcl2a1c',
           'Sfpi1', 'Irf8', 'Irf5', 'Irf4', 'Batf', 'Runx1', 'Ikzf1', 'Cebpa', 'Mlxipl', 'Batf3')
micro <- unique(micro)

endo <- c('Cldn5', 'Ttr', 'Ly6a', 'Madcam1', '8430408G22Rik', 'Akr1c14', 'Ly6c2', 'Meox1',
          'Car4', 'Bsg', 'Aplnr', 'Sigirr', 'Slco1a4', 'Slc16a1', 'Icam2', 'Kank3', 'Slc19a3',
          'Farm101b', 'Slc16a4', 'Nostrin', 'Sdpr', 'Ptgis', 'Myct1', 'Vwa1', 'Ankrd37', 'Sox18',
          'Prnd', 'Dok4', 'Serpinb6b', 'Efna1', 'Cd34', 'Egfl7', 'Pglyrp1', 'Slc35f2', 'Cdkn2b',
          'Fam129a', 'Sgms1', 'Flt1', 'Tie1',
          'Erg', 'Sox17', 'Foxq1', 'Mecom', 'Foxf2', 'Sox18', 'Bcl6b', 'Sox7', 'Meox1', 'Zic3')
endo <- unique(endo)

peri <- c('Fmod', 'Rps2', 'Igf2', 'Gpc3', 'Ogn', 'Lrrc32', 'Flnc', 'Gjb2', 'Itih2', 'Rdh10', 'Bmp6',
          'Aldh1a2', 'Postn', 'Sidt1', 'Lamc3', 'Slc22a6', 'Clec3b', 'Slc6a13', 'Bicc1', 'S100a10', 'Rps18',
          'Serping1', 'Col1a1', 'Dcn', 'Col1a2', 'Pcolce', 'Cyp1b1', 'Cited1', 'Emp1', 'C4b', 'Ahnak',
          'S1pr3', 'Col3a1', 'Fstl1', 'Col4a5', 'Vtn', 'Lama2', 'Mfap4', 'Kcne4', 'Errfi1',
          'Tbx15', 'Foxc2', 'Twist1', 'Tbx18', 'Foxd1', 'Fosl1', 'Heyl', 'Hic1', 'Foxc1', 'Prrx2')
peri <- unique(peri)

cellular_type <- list(OPC= OPC, NFO=NFO, MO= MO,astro = astro, micro= micro,
                      endo=endo, peri=peri)

# stack the list and get a dataframe with the name of cell type and a list of the genes

cellular_type <- stack(cellular_type)
colnames(cellular_type)<- c('gene_symbol', 'gs_name')

cellular_type_gsea <- cellular_type %>% group_by(gs_name)%>% summarise(gene_symbol = list(gene_symbol))

cellular_type_gsea

# gs_name     gene_symbol
# <fct>       <list>     
#   1 astrocytes  <chr [47]> 
#   2 OPC         <chr [49]> 
#   3 NFO         <chr [50]> 
#   4 MO          <chr [50]> 
#   5 microglia   <chr [47]> 
#   6 endothelial <chr [47]> 
#   7 pericyte    <chr [50]> 


setwd('C:/Users/usuario/Desktop/RESULTADOS/RNAseq_OLG/')

plots_dir <- 'plots'
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir)
}
results_dir <- 'results/GSEA'
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
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

datasets <- list(koDvswtD=koDvswtD, koPvswtP=koPvswtP, koDvskoP=koDvskoP, wtDvswtP=wtDvswtP)



plots_dir_go <- 'plots/GSEA'
if (!dir.exists(plots_dir_go)) {
  dir.create(plots_dir_go)}

go_plot <- function(go_result, filename, type='', sign = ''){
  # go_result - output from GSEA@result function, use the enrichplot function to generate fast graphs
  
  # filaname, type and sign are just for create the output file name
  # filename - name of the comparison
  # type - BP, MF or CC
  # sign - UP or DOWN. both as default
  
  options(enrichplot.colours = c("blue","red")) #change enrichplot colours
  

  plots_dir_go_unic <- paste('plots/GSEA/', filename, sep='')
  if (!dir.exists(plots_dir_go_unic)) {
    dir.create(plots_dir_go_unic)}
  
  print(paste(filename,'_', type, '.pdf', sep=''))
  
  pdf(file.path(plots_dir_go_unic,paste(type, '_',filename, '_',sign, '.pdf', sep ='')), width = 13, height = 8)
  
  tryCatch({
    bar <- barplot(go_result, showCategory = 10)+
      ggtitle(paste(sign, 'DEG'))+
      theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
    print(bar)
  }, error =function(e){message("Error al crear el barplot: ", e$message)})
  
  tryCatch({
    enriches <- pairwise_termsim(go_result)
    print(
      treeplot(enriches, color= 'NES') + labs(color='NES')+
        ggtitle(paste(type, '- Treeplot', filename, type)) +
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
  options(enrichplot.colours = c("blue","red"))
}



# Cell Type GSEA -------------------------------------------------------------------------------------------------------------------------------

results_dir <- 'results/GSEA/CelullarType'
if (!dir.exists(results_dir)) {
  dir.create(results_dir)
}

results_list <- list()

x=0 #contador para acceder al nombre del dataset
for (dataset in datasets){ #hacer todos los datasets a la vez
  x=x+1
  name <- names(datasets)[x]
  print(name)
  
  #Crear un vector con los rankings en orden decreciente con los Gene ID como nombres
  
  # Create a vector and prerank the genes in descending order by logFC
  rank <- dataset[[3]]
  names(rank) <- koDvswtD$Gene_ID
  rank <- sort(rank, decreasing = TRUE) 
  
  #GSEA
  gsea_results <- GSEA(
    geneList = rank, # Ordered ranked gene list
    minGSSize = 5, # Minimum gene set size
    maxGSSize = 500, # Maximum gene set set
    pvalueCutoff = 1, # p-value cutoff , 1 to get all
    nPerm = 1000,
    eps = 0, 
    seed = TRUE, 
    pAdjustMethod = "BH", # Benjamini-Hochberg correction
    TERM2GENE = cellular_type_gsea )
  
  gsea_result_dataset <- data.frame(gsea_results@result)
  
  go_plot(gsea_result_dataset, filename = name, type = 'Cell Type')
  results_list[[name]] <- gsea_results
  
  readr::write_tsv(
    gsea_result_dataset,
    file.path(
      results_dir, 
      paste(name,'_gsea_cell_results.tsv',sep = '')
    )
  )
}

# GSEA plot
pdf('plots/GSEA/Diff_Celltype.pdf')
gseaplot(results_list[["koDvswtD"]], geneSetID = '3', color.line = 'red', title = 'Newly Formed OLG') 
gseaplot(results_list[["koDvswtD"]], geneSetID = '4', color.line = 'red', title = 'Myelinating OLG') 
gseaplot(results_list[["koDvswtD"]], geneSetID = '7', color.line = 'green', title = 'Pericyte') 
gseaplot(results_list[["koDvswtD"]], geneSetID = '5', color.line = 'green', title = 'Microglia')
gseaplot(results_list[["koDvswtD"]], geneSetID = '1', color.line = 'black', title = 'Astrocyte')
gseaplot(results_list[["koDvswtD"]], geneSetID = '2', color.line = 'black', title = 'OLG precursor cell')
gseaplot(results_list[["koDvswtD"]], geneSetID = '6', color.line = 'black', title = 'Endothelial')
dev.off()

pdf('plots/GSEA/Prol_Celltype.pdf')
gseaplot(results_list[["koPvswtP"]], geneSetID = '3', color.line = 'green', title = 'Newly Formed OLG') 
gseaplot(results_list[["koPvswtP"]], geneSetID = '4', color.line = 'green', title = 'Myelinating OLG') 
gseaplot(results_list[["koPvswtP"]], geneSetID = '7', color.line = 'green', title = 'Pericyte') 
gseaplot(results_list[["koPvswtP"]], geneSetID = '5', color.line = 'green', title = 'Microglia')
gseaplot(results_list[["koPvswtP"]], geneSetID = '1', color.line = 'green', title = 'Astrocyte')
gseaplot(results_list[["koPvswtP"]], geneSetID = '2', color.line = 'green', title = 'OLG precursor cell')
gseaplot(results_list[["koPvswtP"]], geneSetID = '6', color.line = 'green', title = 'Endothelial')
dev.off()



# Dot plot calculation with the NES as color and -log p val as size

graphs <- list()
x=0
for (result in results_list){
  x=x+1
  name <- names(results_list)[x]
  print(name)
  result_data <-  results_list[[name]]@result[,1:7] %>%dplyr::select(ID, NES, p.adjust)%>% 
    mutate(ID = case_when(ID == 1 ~ "astrocyte",ID == 2 ~ "OPC",ID == 3 ~ "NFO",ID == 4 ~ "MO",ID == 5 ~ "microglia", ID == 6 ~ "endothel",
                          ID == 7 ~ "pericyte")) %>% mutate(log10Padj = -log10(p.adjust))
  
  p <- ggplot(result_data, aes(x=ID, y=1))+
    geom_point(aes(size= log10Padj, color = NES))+
    scale_size_continuous(name = "-log10(p-value)", range = c(2, 10), limits = c(1, 5)) + # Point size as log10pvalue 
    scale_color_gradient(low = "blue", high = "red", name = "NES",limits = c(-3, 3)) + # Color gradient as NES
    scale_y_continuous(expand = c(0,0))+
    ggtitle(name) +
    theme_bw() +
    theme(axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          aspect.ratio = 0.1,
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "top")
  
  graphs <- append(graphs, list(p))
}


pdf('plots/GSEA/NES_Overlap.pdf')
graphs[[1]]
graphs[[2]]
graphs[[3]]
graphs[[4]]
dev.off()

#heatmap de genes de tipos celulares

#genes_celltype <- unique(c(astrocytes, OPC, NFO, MO, microglia, endothelial, pericyte))

genes_celltype <- unique(c(OPC, NFO, MO))

tpm <- read.table('ReadCount/tpm.txt') #tpm

tpm_cell <- tpm %>% filter(gene %in% genes_celltype) #get the interest genes
head(tpm_cell)

cellular_type_nd <-  cellular_type %>% filter(gene_symbol %in% tpm_cell$gene) %>% distinct(gene_symbol, .keep_all = TRUE) %>%
  column_to_rownames(var='gene_symbol')%>%rename(gs_name= 'cell_type') 

tpm_cell$gene <- NULL
tpm_cell_log <- log2(tpm_cell + 1) #normalización

tpm_cell_log <- tpm_cell_log[match(rownames(cellular_type_nd), rownames(tpm_cell_log)), ] # order the genes by the order of cell type
tpm_cell_log <- tpm_cell_log[, c("wt_prol1","wt_prol2", "wt_prol3", 
                                 "ko_prol2", "ko_prol3", 
                                 "wt_diff1", "wt_diff2", "wt_diff3", 
                                 "ko_diff2", "ko_diff3")]

tpm_cell_plot <- tpm_cell[match(rownames(cellular_type_nd), rownames(tpm_cell)), ] # order the genes by the order of cell type
tpm_cell_plot <- tpm_cell_plot[, c("wt_prol1","wt_prol2", "wt_prol3", 
                                   "ko_prol2", "ko_prol3", 
                                   "wt_diff1", "wt_diff2", "wt_diff3", 
                                   "ko_diff2", "ko_diff3")]

annotation_colors <- list()

gaps_row <- cellular_type_nd %>%
  mutate(row = row_number()) %>%          # Añadir el número de fila
  group_by(cell_type) %>%
  summarise(last_row = max(row))  %>%     # Agrupar por `cell_type`
  pull(last_row)                          # Obtener la última fila de cada grupo


head(gaps_row)

annotation_col <- data.frame(
  Genotype = c(rep("wt", 3), rep("ko", 2), rep("wt", 3), rep("ko", 2)),
  Condition = c(rep("proliferation", 5), rep("differentiation", 5)))  # Segunda anotación: wt vs ko 
rownames(annotation_col) <- colnames(tpm_cell_plot)

annotation_colors <- list(
  Genotype = c("wt" = "#66c2a5", "ko" = "red"),
  Condition = c("proliferation" = "#60dbe8", "differentiation" = "#8bd346"),
  cell_type = c(
    # "astrocytes" = "#9b5fe0",
    "OPC" = "#16a4d8",
    "NFO" = "#60dbe8",
    "MO" = "#8bd346",
    # "microglia" = "#efdf48",
    # "endothelial" = "#d64e12",
    # "pericyte" = "#f9a52c"
  ))


p <- pheatmap(tpm_cell_log,
              scale = 'row',
              cluster_rows = FALSE,
              annotation_row = cellular_type_nd,
              annotation_col = annotation_col,
              annotation_colors = annotation_colors,
              cluster_cols = FALSE,
              show_rownames = FALSE,
              show_colnames = FALSE,
              gaps_col= c(5),
              gaps_row = gaps_row, 
              border_color = NA,
              color =viridis(10, option = 'inferno'), name = "log2 (TPM)")+ ggtitle('log2TPM+1')

ggsave('plots/HEATMAP_CELLTYPE.png', plot = p, height = 7, width = 7)


pdf('plots/HEATMAP_CELLTYPE.pdf')
print(p)
dev.off()



