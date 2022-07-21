## ---- include = FALSE---------------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.dpi = 1200
)

## ---- eval=T, echo=T, warning=F, message=F------------------------------------
# Load the required packages
require(Biostrings)
require(CProtMEDIAS)

# Some sample data
data("Homeobox_small")
data("Homeobox_small_annot")

## ----eval=T, results = 'hide',echo=T, warning=F, message=F,fig.width = 5,fig.height = 5,fig.align = "center"----
##Combine alignment results in order of family homology 
Homeobox_small_all <- combine_alignment(alignment_list = Homeobox_small,tree = T)

## ----eval=T,echo=T, warning=F, message=F,fig.width = 5,fig.height = 5,fig.align = "center"----
##Show family homology order 
Homeobox_small_all$tree$labels <- names(Homeobox_small)

plot(Homeobox_small_all$tree)

##Combined alignment results
Homeobox_small_all$combined_alignment


## ----echo=T, fig.align="center",fig.height=3,fig.width=6,message=FALSE, warning=FALSE----

##Alignment score matrix
alignment_score <- get_alignment_score(alignment = Homeobox_small_all$combined_alignment,type = "AA")

require(pheatmap)

labels_row <- rep("",500)

labels_row[c(50,150,250,350,450)] <- names(Homeobox_small)

##Heatmap of Alignment score matrix
pheatmap(alignment_score,labels_row = labels_row,show_colnames = F,cluster_rows = F,cluster_cols = F)


## ----eval=T, echo=T, warning=F, message=F-------------------------------------
##Seurat workflow
my_Seurat <- Seurat_workflow(alignment_score,dims_UMAP = 1:5,resolution = 0.05)

##Import metadata from data frame (annotation)
my_Seurat <- import_metadata(my_Seurat,Homeobox_small_annot)


## ---- message=FALSE, warning=FALSE, out.width="45%",fig.show = "hold"---------
require(ggsci)
require(ggplot2)
##Scatter plot 
##Family
my_color1 <- pal_nejm(palette = "default")(5)

p1 <- cluster_scatter(my_Seurat,group = my_Seurat$Family,text_size = 3) + 
  
      ##set color of each group
      scale_color_manual(values = my_color1) + 
  
      ##set genepro custom theme
      theme_scatter(base_size = 7) + 
  
      ##set position of legend
      theme(legend.position = c(1,0),legend.justification = c(1,0),
            
            legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
      
      ##set number of column and point size of legends
      guides(color = guide_legend(ncol = 2,title = "Family",override.aes = list(size = 1.5)))

p1

##Cluster
my_color2 <- pal_jama(palette = "default")(7)

p2 <- cluster_scatter(my_Seurat,group = my_Seurat$seurat_clusters,text_size = 3) + 
  
      scale_color_manual(values = my_color2) + 
  
      theme_scatter(base_size = 7) + 
  
      theme(legend.position = c(1,0),legend.justification = c(1,0),
            
            legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
      
      guides(color = guide_legend(ncol = 2,title = "Family",override.aes = list(size = 1.5)))

p2

## ----fig.align="center", message=FALSE, warning=FALSE,fig.width = 6,fig.height = 4----

require(grDevices)
require(RColorBrewer)

compare_classification <- table(my_Seurat$Family,my_Seurat$seurat_clusters)

my_color3 <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(101)

compare_classification_color <- my_color3[round(compare_classification / 2) + 1]
  
pheatmap(compare_classification,color = "white",legend = F,
         
         cluster_rows = F,cluster_cols = F,angle_col = 0,
         
         display_numbers = T,number_format = "%d",
         
         number_color = compare_classification_color,
         
         fontsize = 12,fontsize_number = 12,fontface = "bold")
  

## ---- message=FALSE, warning=FALSE, out.width="45%",fig.show = "hold"---------

##Species classification

##Combine monocot, eudicot and Basal Magnoliophyta into angiosperms 
table(my_Seurat$Taxonomic1)

my_Seurat$new_Taxonomic <- as.factor(my_Seurat$Taxonomic1)

levels(my_Seurat$new_Taxonomic)

levels(my_Seurat$new_Taxonomic)[c(1,5,8)] <- "Angiospermae"

levels(my_Seurat$new_Taxonomic)

my_color3 <- scale_color_manual(values = c("#DDDDDD",pal_nejm(palette = "default")(3)))

my_color3 <- c("#DDDDDD",pal_ucscgb(palette = "default")(5))

my_size <- ifelse(my_Seurat$new_Taxonomic == "Angiospermae",1,2)

p3 <- cluster_scatter(my_Seurat,group = my_Seurat$new_Taxonomic,
                      
                      text_group = my_Seurat$Family,text_color = "#00000055",
                      
                      size = my_size,text_size = 3) +
  
      theme_scatter(base_size = 7) + scale_color_manual(values = my_color3) +
  
      theme(legend.position = c(1,0),legend.justification = c(1,0),
            
            legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
  
      guides(color = guide_legend(title = "Species",override.aes = list(size = 1.5)))

p3

##Combine Bryophyta, Chlorophytae, Coniferophyta, Lycopodiophyta, Marchantiophyta into Other plants

table(my_Seurat$Taxonomic2)

my_Seurat$new_Taxonomic2 <- as.factor(my_Seurat$Taxonomic2)

levels(my_Seurat$new_Taxonomic2)

levels(my_Seurat$new_Taxonomic2)[c(3,4,5,8,10)] <- "Other plants"

levels(my_Seurat$new_Taxonomic2)[4] <- "Other Eudicots"

my_color4 <- pal_jama(palette = "default")(7)

p4 <- cluster_scatter(my_Seurat,group = my_Seurat$new_Taxonomic2,
                      
                      text_group = my_Seurat$Family,
                      
                      text_color = "#00000055",text_size = 3) +
  
      theme_scatter(base_size = 7) + scale_color_manual(values = my_color4) +
  
      theme(legend.position = c(1,0),legend.justification = c(1,0),
            
            legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
  
      guides(color = guide_legend(title = "Species",override.aes = list(size = 1.5)))

p4

## ----fig.align="center", message=FALSE, warning=FALSE, out.width="55%"--------

##WOX Subfamily

table(my_Seurat$WOX_subgroup)

my_color5 <- c("#DDDDDD",pal_npg(palette = "nrc")(3))

my_color5 <- my_color5[c(3,2,1,4)]

my_Seurat$WOX_subgroup[is.na(my_Seurat$WOX_subgroup)] <- "Other family"

p5 <- cluster_scatter(my_Seurat,group = my_Seurat$WOX_subgroup,
                      
                      text_group = my_Seurat$Family,text_color = "#00000055",
                      
                      text_size = 3,size = .8) +
  
      theme_scatter(base_size = 7) + scale_color_manual(values = my_color5) +
  
      theme(legend.position = c(1,0),legend.justification = c(1,0),
            
            legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
  
      guides(color = guide_legend(title = "Species",override.aes = list(size = 1.5)))

p5

## ---- message=FALSE, warning=FALSE, out.width="45%",fig.show = "hold"---------
##Gene length

##Calculate the length of each sequence
tmp_list <- strsplit(gsub(".*/","",rownames(Homeobox_small_annot)),split = "-")

my_Seurat$seq_length <- as.numeric(unlist(lapply(tmp_list,"[",2))) -
  
              as.numeric(unlist(lapply(tmp_list,"[",1))) + 1

midpoint <- sum(range(my_Seurat$seq_length)) / 2

hist(my_Seurat$seq_length,main = "Sequence length")

p6 <- cluster_scatter(my_Seurat,group = my_Seurat$seq_length,text = T,
                      text_group = my_Seurat$Family,text_size = 3) + 
  
      scale_color_gradient2(low = "#4575B4",high = "#D73027",
                            mid = "#FEFEC0",midpoint = midpoint) +
  
      theme_scatter(base_size = 7) +
  
      theme(legend.position = c(1,0),legend.justification = c(1,0))+ 
  
      guides(color = guide_colourbar(title = "Length\n"))

p6


## ----eval=T, echo=T, warning=F, message=F, fig.height=3, fig.width=3,fig.align="center"----
##Monocle workflow to get pseudotime and state of each sequence
require(monocle)

my_monocle <- monocle_workflow(my_Seurat)

table(my_monocle$State)

hist(my_monocle$Pseudotime,main = "Pseudotime")

## ---- message=FALSE, warning=FALSE, out.width="45%",fig.show = "hold"---------

##Colored by State
my_color7 <- pal_aaas(palette = "default")(4)

p7 <- cluster_scatter(my_monocle,group = my_monocle$State,
                      
                      text_group = my_monocle$Family,text_size = 3) +
  
  theme_scatter(base_size = 7) + 
  
  theme(legend.position = c(1,0),legend.justification = c(1,0),
            
        legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
  
  scale_color_manual(values = my_color7) + 
  
  guides(color = guide_legend(ncol = 3,title = "State",override.aes = list(size = 1.5))) +
  
  labs(x = "Component 1",y = "Component 2")

p7

##Colored by Family
p8 <- cluster_scatter(my_monocle,group = my_monocle$Family,
                      
                      text_group = my_monocle$Family,text_size = 3) +
  
  theme_scatter(base_size = 7) + 
  
  theme(legend.position = c(1,0),legend.justification = c(1,0),
            
        legend.key.height = unit(0.2,"cm"),legend.key.width = unit(0.2,"cm")) + 
  
  scale_color_manual(values = my_color1) + 
  
  guides(color = guide_legend(ncol = 2,title = "State",override.aes = list(size = 1.5))) +
  
  labs(x = "Component 1",y = "Component 2")

p8

## ----eval=T, echo=T, warning=F, message=F-------------------------------------
##Group must in the metadata of Seurat object
colnames(my_Seurat@meta.data)

##We chose 'Family' as costum group, identifying the specific site of each family
my_Specific_site <- find_Specific_site(my_Seurat,ident = "Family",logfc = 0,min.pct = 0)

head(my_Specific_site)

## ----fig.align="center", message=FALSE, warning=FALSE, fig.height = 2,fig.width = 6----
##Get seq conservation (bits) of each site in each family (group)
my_seq_conservation <- get_seq_conservation(seqs = Homeobox_small_all$combined_alignment,
                                            
                                            group = my_Seurat$Family,type = "AA")
##Conservative heatmap for each site
p9 <- Position_site_heatmap(my_seq_conservation,group_color = my_color1) + 
  
                            theme(plot.margin = margin(5,20,5,20))

p9

## ----fig.align="center", message=FALSE, warning=FALSE, fig.height = 3,fig.width = 7----

##Specific site of each family
p10 <- Position_site_plot(dat = my_Specific_site,seurat = my_Seurat,
                          
                          site_color = my_color1,plot_margin = margin(5,5,60,5))

p10

## ----fig.align="center",fig.height = 6,fig.width = 6, message=FALSE, warning=FALSE----

##The 20 most conserved specific sites of each family are selected to draw weblogo
p11 <- position_site_weblogo(seqs = Homeobox_small_all$combined_alignment,
                             
                             Specific_site = my_Specific_site,text_angle = 45,max_seq_length = 20,
                             
                             group = my_Seurat$Family,stack_width_threshold = 0.1)

p11

##Specify the same site for each family to draw a weblogo 
p12 <- position_site_weblogo(seqs = Homeobox_small_all$combined_alignment,range = 91:110,
                             
                             max_seq_length = 20,
                             
                             Specific_site = my_Specific_site,text_angle = 45,
                             
                             group = my_Seurat$Family,stack_width_threshold = 0.1)

p12

## ---- message=FALSE, warning=FALSE--------------------------------------------

require(igraph)

table(my_Seurat$Family,my_Seurat$Taxonomic1)

##Order of species evolution: Chlorophytae -> Charophyta -> Marchantiophyta -> Bryophyta -> 
##Lycopodiophyta -> Coniferophyta -> Basal Magnoliophyta -> Monocots ->Eudicots
##We chose HB-PHD as root node since only HB-PHD family contain Chlorophytae genes.

Family_lineage<- constructingNetwork(alignment_score,group = my_Seurat$Family,rootNode = 2)


## ----fig.align="center",fig.height = 6,fig.width = 3, message=FALSE, warning=FALSE----

## Visualization of Gene lineage relationships
require(igraph)

MDST_plot(Family_lineage$MDST,node.color = my_color1,
          
          group = as.character(Family_lineage$group),
          
          edge.label = round(Family_lineage$MDST[,3],2),
          
          edge.label.adjust = 0.15,weight = F,rescale = T)


## ----fig.align="center",fig.height = 6,fig.width = 3, message=FALSE, warning=FALSE----

Species_lineage <- constructingNetwork(alignment_score,group = my_Seurat$new_Taxonomic2,rootNode = 3)
## Visualization of species lineage relationships
MDST_plot(Species_lineage$MDST,node.color = my_color4,
          
          group = as.character(Species_lineage$group),
          
          edge.label = round(Species_lineage$MDST[,3],2),
          
          edge.label.adjust = 0.15,weight = F,rescale = F)


