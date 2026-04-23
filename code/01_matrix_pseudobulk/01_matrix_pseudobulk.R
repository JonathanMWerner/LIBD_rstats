
#Code for using matrix operations for quick pseudobulking of single cell data

#Using SingleCellExperiment objects
#All you need is the expression data and the cluster labels/annotations


#Remember to order the Variable Pane by size to see the memory usage of the objects we're creating


library(SingleCellExperiment)
library(qs2)
library(here)

here()

plot_path = here('plots', '01_matrix_pseudobulk')


multiome_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/05_03_annotation_adjustments/06_refined_annotations'
multiome_sce = qs_read(paste0(multiome_path, '/refined_annotation_multiomeHab_SCE.qs2'))

#Drop the ATAc for memory
altExp(multiome_sce, "ATAC") <- NULL
gc()

#What is the dataset size here
multiome_sce

#dim: 36601 54193 


#What annotation do we want to use for pseudobulking?
colnames(colData(multiome_sce))



#Go ahead with the refined_mid_cluster annotation
#And the refined_cluster_ann annotation

table(multiome_sce$refined_mid_cluster)

# 17 cluster annotations
#
#    Astrocyte          Endo     Ependymal    Excit.Thal Inhib_LHb_4.1 Inhib_LHb_4.2 
#         1343           343          1341          9738          1123           932 
#   Inhib.Thal     LHb.1.3.4       LHb.2.7         LHb.4         MHb.1       MHb.1.2 
#         5934          2990          2771         10735          4370          1607 
#        MHb.2         MHb.3     Microglia         Oligo           OPC 
#         4822           145           663          4698           638 


table(multiome_sce$refined_cluster_ann)

# 37 cluster annotations
#
#C.01.Inhib.Thal      C.02.Oligo C.03.Excit.Thal      C.04.LHb.4    C.05.LHb.2.7 
#           3906            3371            3271            2317            2771 
#     C.06.LHb.4      C.07.MHb.2      C.08.LHb.4      C.09.LHb.4      C.10.MHb.1 
#           2012            2610            2513            2256            2344 
#   C.11.MHb.1.2 C.12.Excit.Thal      C.13.LHb.4      C.14.MHb.1 C.15.Excit.Thal 
#           2212            2187            1455            2026            1625 
#   C.16.MHb.1.2 C.17.Excit.Thal  C.18.LHb.1.3.4 C.19.Inhib.Thal  C.20.Astrocyte 
#           1607            1587            1535            1380            1343 
# C.21.Ependymal      C.22.Oligo      C.23.LHb.1 C.25.Excit.Thal        C.26.OPC 
#           1341            1327            1269             707             638 
# C.27.Microglia C.28.Inhib.Thal       C.29.Endo      C.31.LHb.4 C.32.Excit.Thal 
#            587             543             343             182             196 
#   C.33.LHb.1.3 C.35.Excit.Thal      C.36.MHb.3 C.38.Inhib.Thal  C.41.Microglia 
#            186             165             145             105              76 
#  Inhib_LHb_4.1   Inhib_LHb_4.2 
#           1123             932 

#Get the one-hot encoding of the cluster annotations for the multiome data, in this case the midcluster metadata column
cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_mid_cluster, data = colData(multiome_sce))

dim(cell_annot_matrix)
# 54193    17

head(cell_annot_matrix)
#Column names are hidden, but shows how cells are encoded as 1s for their matched cluster annotaion
#
#S04_AAACAGCCAGAATGAC-1 . . . 1 . . . . . . . . . . . . .
#S04_AAACAGCCAGCAAGGC-1 . . . . . . . . . 1 . . . . . . .
#S04_AAACATGCACCTGGTG-1 . . . . . . . . . 1 . . . . . . .
#S04_AAACATGCAGGATGGC-1 . . . . 1 . . . . . . . . . . . .
#S04_AAACATGCAGTAATAG-1 . . . . . . 1 . . . . . . . . . .
#S04_AAACATGCATAAGTCT-1 . . . . . . 1 . . . . . . . . . .

#Drops the 'mid_cluster' prefix on the column names
colnames(cell_annot_matrix) <- gsub("refined_mid_cluster", "", colnames(cell_annot_matrix))

#Get the pseudobulk counts
start = Sys.time()
pseudobulk_counts = assay(multiome_sce, 'counts') %*% cell_annot_matrix
Sys.time() - start

#Time difference of 0.4109297 secs

dim(pseudobulk_counts)
# 36601    17

head(pseudobulk_counts)

#MIR1302-2HG . .  .   8  3  .  1  1  .   7  2  .  . . . 1 .
#FAM138A     . .  .   .  .  .  .  .  .   .  .  .  . . . . .
#OR4F5       . .  .   .  .  .  .  .  .   2  .  .  . . . . .
#AL627309.1  5 . 14 253 23 15 50 49 29 187 29 12 27 . 2 7 .
#AL627309.3  . .  2   6  2  2  .  2  .   7  1  1  1 . . . .
#AL627309.2  . .  .   .  .  .  1  .  .   1  .  .  . . . . .

#
#Sanity check
#
gene_subset = c('MIR1302-2HG','FAM138A','OR4F5','AL627309.1','AL627309.3', 'AL627309.2')
gene_index = rownames(multiome_sce) %in% gene_subset

#LHb.4 cell-type
cell_subset = colnames(pseudobulk_counts)[10]
cell_index = multiome_sce$refined_mid_cluster == cell_subset

rowSums(assay(multiome_sce, 'counts')[gene_index, cell_index])
#MIR1302-2HG     FAM138A       OR4F5  AL627309.1  AL627309.3  AL627309.2 
#          7           0           2         187           7           1 




#
# And now for the fine resolution annotations
#

cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_cluster_ann, data = colData(multiome_sce))

dim(cell_annot_matrix)
# 54193    37

colnames(cell_annot_matrix) <- gsub("refined_cluster_ann", "", colnames(cell_annot_matrix))

#Get the pseudobulk counts
start = Sys.time()
pseudobulk_counts = assay(multiome_sce, 'counts') %*% cell_annot_matrix
Sys.time() - start

#Time difference of 0.4491587 secs


#
# Now make a large sce object and see how long it takes to pseudobulk with increasing the number of cells
#



#Make a large SCE Object by adding more cells to the multiome SCE object by rbind-ing it to itself multiple times
append_sce_to_itself <- function(sce, times = 3L) {
  stopifnot(times >= 1L)

  sce_list <- lapply(seq_len(times), function(i) {
    x <- sce
    colnames(x) <- paste0(colnames(sce), "_rep", i)  # keep cell IDs unique
    x
  })

  do.call(cbind, sce_list)
}

# Example: make object 4x larger (same genes, 4x cells)
sce_big <- append_sce_to_itself(multiome_sce, times = 4L)

dim(multiome_sce)
dim(sce_big)

# 36601 54193
#  36601 216772


num_cells = seq(10000, 201000, by = 10000)

time_vec = vector(mode = 'numeric', length = length(num_cells))

for(i in 1:length(num_cells)){
  sce_subset <- sce_big[, 1:num_cells[i]]
  
  cell_annot_matrix <- Matrix::sparse.model.matrix(~ 0 + refined_cluster_ann, data = colData(sce_subset))
  
  colnames(cell_annot_matrix) <- gsub("refined_cluster_ann", "", colnames(cell_annot_matrix))
  
  start = Sys.time()
  pseudobulk_counts = assay(sce_subset, 'counts') %*% cell_annot_matrix
  time_vec[i] = as.numeric(Sys.time() - start, units = "secs")

  #Free up memory
  rm(pseudobulk_counts, cell_annot_matrix, sce_subset)
  gc()
  print(i)
}


plot(num_cells, time_vec, type = 'b', xlab = 'Number of Cells', ylab = 'Time (seconds)', main = 'Time taken for pseudobulking with increasing number of cells')


pdf(file.path(plot_path, 'pseudobulk_time_plot.pdf'), width = 6, height = 4)
plot(num_cells, time_vec, type = 'b', xlab = 'Number of Cells', ylab = 'Time (seconds)', main = 'Time taken for pseudobulking with increasing number of cells')
dev.off()




