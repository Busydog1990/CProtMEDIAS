
#===============================================================================
#' Scatter plot for genes
#'
#' Scatter plot for alignment result after dimensionality reduction.
#'
#' @param dat Seurat/monocle object.
#' @param mapping Name of reduction to pull cell embeddings for. Only for Seurat object.
#' @param group Factor. Grouping information of each gene. Must contain the same length of genes.
#' @param size Numeric. The point size of scatter plot. Default:1.
#' @param text Logical. If TRUE, gene family will be labeled at scatter plot. Default:T.
#' @param text_group Factor. Grouping information of labeled gene family. Must contain the same length of genes.
#' @param text_color Character. Colors of labeled gene family. Default:"black".
#' @param text_size Numeric. The gene family text size of scatter plot. Default:5.
#' @param segment_size Numeric. The segment size of scatter plot. Default:0.75.
#' @param text_face Character. fontface of labeled gene family.
#' @import ggplot2
#' @import ggrepel
#' @import stats
#' @return Scatter plot for genes.
#' @export
##Scatter plot for genes
#-------------------------------------------------------------------------------
cluster_scatter <- function(dat,mapping = "umap",group = NULL,size = 1,text = T,
                            text_group = NULL,text_color = "black",
                            text_size = 5,segment_size = 0.75,text_face = "bold"){
  if ("Seurat" %in% class(dat) & attr(class(dat),"package") == "SeuratObject"){
    plot_dat <- data.frame(SeuratObject::Embeddings(dat,reduction = mapping))
    colnames(plot_dat) <- c("UMAP_1","UMAP_2")
    if (!is.null(group)){
      plot_dat$group <- group
      p1 <- ggplot(plot_dat,aes(x = UMAP_1,y = UMAP_2,color = group)) + geom_point(size = size)
    } else {
      p1 <- ggplot(plot_dat,aes(x = UMAP_1,y = UMAP_2)) + geom_point(size = size)
      return(p1)
    }
    if (text){
      if (is.null(text_group)){
        text_dat <- merge(stats::aggregate(data = plot_dat,UMAP_1~group,mean),
                          stats::aggregate(data = plot_dat,UMAP_2~group,mean))
        p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = UMAP_1,y = UMAP_2,label = group),
                                            color = text_color,fontface = text_face,size = text_size)
      } else {
        plot_dat$text_group <- text_group
        text_dat <- merge(aggregate(data = plot_dat,UMAP_1~text_group,mean),
                          aggregate(data = plot_dat,UMAP_2~text_group,mean))

      p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = UMAP_1,y = UMAP_2,label = text_group),
                                          color = text_color,fontface = text_face,size = text_size)
      }
    }
    return(p1)
  } else if ("CellDataSet" %in% class(dat) & attr(class(dat),"package") == "monocle"){
      plot_dat <- ggplot_build(plot_cell_trajectory(dat))
      segment_dat <- plot_dat$data[[1]]
      point_dat <- plot_dat$data[[2]]
      node_dat <- plot_dat$data[[3]]
      node_text_dat <- plot_dat$data[[4]]

      if (!is.null(group)){
        point_dat$group <- group
        p1 <- ggplot(point_dat,aes(x = x,y = y,color = group)) + geom_point(size = size)
        p1 <- p1 + geom_segment(data = segment_dat,aes(x = x,y = y,xend = xend,yend = yend),
                                size = segment_size,color = "black") +
          geom_point(data = node_dat,aes(x = x,y = y),color = I("black"),size = I(5)) +
          geom_text(data = node_text_dat,aes(x = x,y = y,label = label),color = I("white"),size = I(4))
      } else {
        p1 <- ggplot(point_dat,aes(x = x,y = y)) + geom_point(size = size)
        p1 <- p1 + geom_segment(data = segment_dat,aes(x = x,y = y,xend = xend,yend = yend),size = segment_size) +
          geom_point(data = node_dat,aes(x = x,y = y),color = I("black"),size = I(5)) +
          geom_text(data = node_text_dat,aes(x = x,y = y,label = label),color = I("white"),size = I(4))
        return(p1)
      }
      if (text){
        if (is.null(text_group)){
        text_dat <- merge(aggregate(data = point_dat,x~group,mean),
                          aggregate(data = point_dat,y~group,mean))
        p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = x,y = y,label = group),
                                            color = text_color,fontface = text_face,
                                            size = text_size)

        } else {
          point_dat$text_group <- text_group
          text_dat <- merge(aggregate(data = point_dat,x~text_group,mean),
                            aggregate(data = point_dat,y~text_group,mean))
          p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = x,y = y,label = text_group),
                                            color = text_color,fontface = text_face,
                                            size = text_size)
        }
      }
      return(p1)
  }

}
#===============================================================================
#' Theme of Scatter plot
#' genepro custom theme
#'
#' @param base_size font base size
#' @param base_family font base family
#'
#' @export
theme_scatter <- function(base_size=12, base_family=''){
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(axis.line = element_blank(),
          panel.background = element_rect(color = "black"),
          #legend.key.height = unit(0.5,"cm"),
          #legend.position = c(1,0),
          #legend.justification = c(1,0),
          legend.title = element_text(face = "bold"),
          legend.background = element_rect(fill = NA))
}
#===============================================================================
#' Scatter plot for genes belonging to multiple gene families
#'
#' Scatter plot for genes belonging to multiple gene families,
#' Genes belonging to multiple families will be given multiple colors
#'
#' @param dat Seurat/monocle object.
#' @param my_color Gene families colors. Only gene belonging to multiple gene families will be colored.
#' @param mapping Name of reduction to pull gene embedding for. Only for Seurat object.
#' @param group Factor. Grouping information of each gene. Must contain the same length of genes.
#' @param gene Character. Must contain duplicate gene names.
#' @param background_color Colors for gene NOT belonging to multiple gene families.
#' @param text Logical. If TRUE, gene family will be labeled at scatter plot. Default:T.
#' @param text_size Numeric. The gene family text size of scatter plot. Default:5.
#' @param text_group Factor. Grouping information of labeled gene family. Must contain the same length of genes.
#' @param text_face Character. fontface of labeled gene family.
#' @param text_color Character. Colors of labeled gene family. Default:"black".
#' @param point_size Numeric. The background point size of scatter plot. Default:1.
#' @param segment_size Numeric. The segment size of scatter plot. Only for monocle object. Default:1.
#' @param legend Logical. If TRUE, add a legend to the picture. Default:TRUE.
#' @param multiple_point_size Numeric. Point size of gene belonging to multiple gene families. Default:0.3.
#' @import ggplot2
#' @import ggrepel
#' @return Scatter plot for genes belonging to multiple gene families.
#' @export
#-------------------------------------------------------------------------------
multiple_group_scatter <- function(dat,my_color,mapping = "umap",group = NULL,gene = NULL,
                                   background_color = "grey80",text = T,
                                   text_group = NULL,text_face = "bold",
                                   text_color = "grey60",text_size = 5,legend = T,
                                   point_size = 1,segment_size = 1,multiple_point_size = 0.3){
  if ("Seurat" %in% class(dat) & attr(class(dat),"package") == "SeuratObject"){
    plot_dat <- data.frame(SeuratObject::Embeddings(dat,reduction = mapping))
    plot_dat$gene <- gene
    colnames(plot_dat) <- c("UMAP_1","UMAP_2")
    if (is.null(group)){
      plot_dat$group <- dat@active.ident
    } else {
      plot_dat$group <- group
    }
    plot_dat$group <- as.factor(plot_dat$group)
    my_group_levels <- levels(plot_dat$group)
    my_dup <- duplicated(plot_dat$gene)
    stopifnot("No genes belonging to multiple gene families" =
              sum(my_dup) > 0)
    plot_dat_dupplicate <- plot_dat[plot_dat$gene %in% plot_dat$gene[duplicated(plot_dat$gene)],]
    plot_dat_dupplicate_list <- split(plot_dat_dupplicate,plot_dat_dupplicate$gene)
    if (length(background_color) == 1){
      p1 <- ggplot(plot_dat,aes(x = UMAP_1,y = UMAP_2,color = I(background_color))) +
                   geom_point(size = I(point_size))
    } else if(length(background_color) == length(my_group_levels)){
      p1 <- ggplot(plot_dat,aes(x = UMAP_1,y = UMAP_2,color = group)) +
                   geom_point(size = I(point_size)) +
                   scale_color_manual(values = background_color)
    } else {stop("wrong background color length")}
    if (text){
      if (is.null(text_group)){
        text_dat <- merge(aggregate(data = plot_dat,UMAP_1~group,mean),
                          aggregate(data = plot_dat,UMAP_2~group,mean))
        if (length(text_color) == 1){
          p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = UMAP_1,y = UMAP_2,label = group),
                                              color = text_color,fontface = text_face,
                                              size = text_size)
        } else if(length(text_color) == length(my_group_levels)){
          p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = UMAP_1,y = UMAP_2,label = group),
                                              fontface = text_face,inherit.aes = F,color = text_color,
                                              size = text_size)
        }
      } else {
        text_dat <- merge(aggregate(data = plot_dat,UMAP_1~text_group,mean),
                          aggregate(data = plot_dat,UMAP_2~text_group,mean))
        if (length(text_color) == 1){
          p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = UMAP_1,y = UMAP_2,label = text_group),
                                              color = text_color,fontface = text_face,size = text_size)
        } else if(length(text_color) == length(my_group_levels)){
          p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = UMAP_1,y = UMAP_2,label = text_group),
                                              fontface = text_face,inherit.aes = F,color = text_color,
                                              size = text_size)
        }
      }
    }
    if (legend){
      data2 <- data.frame(x = plot_dat[1,1],y = plot_dat[1,2],group = my_group_levels,color = my_color)
      p1 <- p1 + geom_point(data = data2,aes(x = x,y = y,color = color),size = -1) +
        scale_color_manual(values = data2$color) +
        guides(color = guide_legend(override.aes = list(size = 3)))
    }
  } else if ("CellDataSet" %in% class(dat) & attr(class(dat),"package") == "monocle"){
    plot_dat <- ggplot_build(plot_cell_trajectory(dat))
    segment_dat <- plot_dat$data[[1]]
    point_dat <- plot_dat$data[[2]]
    point_dat <- point_dat[,-1]
    point_dat$gene <- gene
    node_dat <- plot_dat$data[[3]]
    node_text_dat <- plot_dat$data[[4]]
    if (is.null(group)){
      point_dat$group <- dat$State
    } else {
      point_dat$group <- group
    }
    point_dat$group <- as.factor(point_dat$group)
    my_group_levels <- levels(point_dat$group)
    my_dup <- duplicated(point_dat$gene)
    stopifnot("No genes belonging to multiple gene families" =
                sum(my_dup) > 0)
    plot_dat_dupplicate <- point_dat[point_dat$gene %in% point_dat$gene[duplicated(point_dat$gene)],]
    plot_dat_dupplicate_list <- split(plot_dat_dupplicate,plot_dat_dupplicate$gene)

    p1 <- ggplot(point_dat,aes(x = x,y = y)) + geom_point(size = I(point_size),color = I(background_color))
    p1 <- p1 + geom_segment(data = segment_dat,aes(x = x,y = y,xend = xend,yend = yend),size = segment_size) +
      geom_point(data = node_dat,aes(x = x,y = y),color = I("black"),size = I(5)) +
      geom_text(data = node_text_dat,aes(x = x,y = y,label = label),color = I("white"),size = I(4))
    if (text){
      if (is.null(text_group)){
        text_dat <- merge(aggregate(data = point_dat,x~group,mean),
                          aggregate(data = point_dat,y~group,mean))
        p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = x,y = y,label = group),
                                            color = text_color,fontface = text_face,
                                            size = text_size)

      } else {
        point_dat$text_group <- text_group
        text_dat <- merge(aggregate(data = point_dat,x~text_group,mean),
                          aggregate(data = point_dat,y~text_group,mean))
        p1 <- p1 + ggrepel::geom_text_repel(data = text_dat,aes(x = x,y = y,label = text_group),
                                            color = text_color,fontface = text_face,
                                            size = text_size)
      }
    }
    if (legend){

      data2 <- data.frame(x = point_dat$x[1],y = point_dat$y[1],group = my_group_levels,color = my_color)

      p1 <- p1 + geom_point(data = data2,aes(x = x,y = y,color = group),size = -1) +
        scale_color_manual(values = data2$color) +
        guides(color = guide_legend(override.aes = list(size = 3)))
    }
    plot_dat <- point_dat
  }

  df <- data.frame(group = levels(plot_dat$group),y = 1,
                   color = my_color)
  pb <- txtProgressBar(style=3)
  for (i in 1:length(plot_dat_dupplicate_list)){
    setTxtProgressBar(pb, i/length(plot_dat_dupplicate_list))
    my_list <- plot_dat_dupplicate_list[[i]]
    a <- nrow(my_list)
    dat <- df[df$group %in% my_list$group,]
    p_tmp <-  ggplot() + geom_bar(data=dat,aes(x = "",y = y,fill = group),stat="identity") +
      coord_polar("y",start=0) +
      scale_fill_manual(values = dat$color) +
      theme_bw()+
      theme(axis.title = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            panel.border = element_blank(),
            panel.background = element_blank(),
            panel.grid = element_blank(),
            plot.margin = margin(0,0,0,0),
            plot.background = element_blank())
    g<-ggplotGrob(p_tmp)
    p1 <- p1 + annotation_custom(g,xmin=my_list[1,1]-multiple_point_size,xmax=my_list[1,1]+multiple_point_size,
                                 ymin=my_list[1,2]-multiple_point_size,ymax=my_list[1,2]+multiple_point_size)
  }
  close(pb)
  return(p1)
}
#===============================================================================
#' Draw a rectangle to reflect the specific site of each family.
#'
#'
#' @param dat Data frame. Specific site data.
#' @param site Character. Site column in dat. Default:"gene".
#' @param cluster Character. Cluster column in dat. Default:"cluster".
#' @param seurat Seurat object.
#' @param site_color Character. Color of each group.
#' @param legend Logical. If TRUE, add a legend to the picture. Default:TRUE.
#' @param plot_margin margin around entire plot (unit with the sizes of the
#'                    top, right, bottom, and left margins).Default:margin(5,5,50,5)
#' @import ggplot2
#' @import ggsci
#' @import ggrepel
#' @return A ggplot of a rectangle to reflect the specific site of each family.
#' @export
#-------------------------------------------------------------------------------
Position_site_plot <- function(dat,site = "gene",
                               cluster = "cluster",seurat,
                               site_color,legend = T,
                               plot_margin = margin(5,5,50,5)){
  site_number <- nrow(seurat)
  my_site <- dat[,c(site,cluster)]
  colnames(my_site) <- c("site","cluster")
  cluster_num <- length(unique(my_site$cluster))
  if (class(my_site$site) != "numeric"){
    my_site$site <- as.numeric(gsub("\\D","",my_site$site))
  }
  stopifnot("special site beyond site number" = max(my_site$site) <= site_number)
  my_site_table <- table(my_site$cluster,my_site$site)
  my_site_table1 <- table(my_site$cluster)
  my_site_table2 <- table(my_site$site)
  my_site_all <- as.numeric(rep(names(my_site_table2),as.numeric(my_site_table2)))
  my_site_start <- c()
  my_site_end <- c()
  for (i in my_site_table2){
    start <- seq(0,1,1/i)
    length(start) <- i
    end <- seq(1,0,-1/i)
    length(end) <- i
    end <- rev(end)
    my_site_start <- c(my_site_start,start)
    my_site_end <- c(my_site_end,end)
  }
  my_site_color <- site_color[1:nrow(my_site_table)]
  my_fill <- rep(my_site_color,ncol(my_site_table))[as.logical(c(my_site_table))]
  my_rect <- data.frame(Family = names(my_site_table1),xmin = 0,xmax = 0,ymin = 0,ymax = 0)
  p1 <- ggplot() + annotate("rect",xmin = 0,xmax = site_number,ymin = 0,ymax = 1,size = 1,
                                                       color = "black",fill = "white") +
    annotate("rect",xmin = my_site_all - 1,xmax = my_site_all,
             ymin = my_site_start,ymax = my_site_end,fill = my_fill)

  if (legend){
    p1 <- p1 + geom_rect(data = my_rect,aes(xmin = xmin,xmax = xmax,ymin = ymin,
                                            ymax =  ymax,fill = Family)) +
      scale_fill_manual(values = site_color)
  }

  p1 <- p1 + scale_y_continuous(expand = c(0,0),limits = c(0,5)) + scale_x_continuous(expand = c(0,0)) +
    theme_void() %+replace%
    theme(axis.text.x = element_text(colour = "black",face = "bold"),plot.margin = plot_margin,
          legend.position = c(0.5,0.3),legend.direction = "horizontal")
  p1
}
#===============================================================================
#' Heatmap for specific site of each family
#'
#' @param seq_conservation Conservation matrix of each site.
#' @param group_color Gene families colors.
#' @param labs_x Character. Title of X lab.
#' @param labs_y Character. Title of Y lab.
#' @importFrom reshape2 melt
#' @importFrom grDevices colorRampPalette
#' @return Heatmap for specific site of each family.
#' @export
#-------------------------------------------------------------------------------
Position_site_heatmap <- function(seq_conservation,group_color = NULL,
                                  labs_x = "",labs_y = "Gene Family"){
  dat <- t(seq_conservation)
  dat <- melt(dat)
  dat$Var3 <- as.numeric(as.factor(dat$Var2))
  if (length(group_color) == length(unique(dat$Var2))){
    my_color_matrix <- as.matrix(stats::aggregate(data = dat,value~Var2,
                                           function(x)cut((x - min(x)) / (max(x) - min(x)),
                                                          breaks = seq(0,1,length.out = 1000),
                                                          labels = F,right = T,include.lowest = T)))

    my_group <- my_color_matrix[,1]
    my_color_matrix <- as.matrix(sapply(data.frame(my_color_matrix[,-1]),as.numeric))
    my_color_matrix <- my_color_matrix[nrow(my_color_matrix):1,]
    rownames(my_color_matrix) <- rev(my_group)
    my_color_database <- do.call(rbind,lapply(group_color,
                                              function(x)colorRampPalette(c("white",x))(1000)))
    my_color_database <- my_color_database[nrow(my_color_database):1,]
    my_color_matrix2 <- matrix(0,nrow(my_color_matrix),ncol(my_color_matrix))
    for (i in 1:nrow(my_color_matrix)){
      my_color_matrix2[i,] <- my_color_database[i,my_color_matrix[i,]]
    }
    my_color_matrix2 <- c(t(my_color_matrix2))
    p1 <- ggplot(dat,aes(x = Var1,y = Var2)) +
      geom_tile(fill = my_color_matrix2) +
      scale_y_discrete(labels = rev(unique(dat$Var2))) + theme_classic() +
      theme(axis.text.y = element_text(color = rev(group_color),face = "bold"))
   } else if (length(group_color) == 1){
    p1 <- ggplot(dat,aes(x = Var1,y = Var2)) + geom_tile() + scale_y_discrete(labels = rev(unique(dat$Var2))) +
      scale_fill_gradient(low = "white",high = group_color) + theme_classic() +
      theme(axis.text.y = element_text(color = group_color,face = "bold"))
    } else if (is.null(group_color)){
    p1 <- ggplot(dat,aes(x = Var1,y = Var2)) + geom_tile() +
      scale_y_discrete(labels = rev(unique(dat$Var2))) + theme_classic() +
      theme(axis.text.y = element_text(color = "black",face = "bold"))
  } else {
    stop("wrong length of group_color")
  }
  p1 <- p1 + labs(x = labs_x,y = labs_y) +
    scale_x_continuous(expand = c(0,0),breaks = c(1,ncol(seq_conservation)),limits = c(-1,ncol(seq_conservation)))+
    theme(axis.ticks.y = element_blank(),
          axis.line = element_blank(),
          axis.text.x = element_text(color = "black",face = "bold"),
          panel.background = element_rect(color = "black")
    )
  p1
}
#===============================================================================
#' weblogo for specific site of each family
#'
#' @param seqs Aligned sequences.
#' @param alignment_file Path of aligned sequences.
#' @param group Factor. Grouping information of each gene. Must contain the same length of genes.
#' @param type Character. Sequences type, one of DNA, RNA and AA.
#' @param max_seq_length Numeric. The maximum number of site plotted in the weblogo.
#' @param ncol Numeric. The number of columns in Weblogo
#' @param range If range = "conservation", Automatic selection of the most conservative sites.(default)
#'              If range = "site", Automatic selection of the most conservative sites from Specific_site data.
#'              If class(range) is data.frame, each column in the range represents the locus selected by each family.
#'              If class(range) is list, each element in the range represents the locus selected by each family.
#'              If is.null(range), first max_seq_length site will be selected by each family.
#'              if class(range) is integer/numeric, corresponding site will be selected by each family.
#' @param stack_width_threshold Numeric. Stack width lower than stack_width_threshold will be changed to 0.
#' @param Specific_site Specific_site data. Only available when range = "site"
#' @param legend Logical. If FALSE, legend of weblogo will be deleted. Default: FALSE.
#' @param stack_width Letter width in Weblogo. If is.null(stack_width), the letter width will be calculated automatically according to the conservatism of the locus.
#'                    If class(stack_width) is data.frame/list, The letter width will be extracted from stack_width.
#'                    If stack_width is a length 1 numeric/integer vector, width of all letters are set to stack_width.
#'                    If stack_width is a numeric/integer vector with length > 1. The letter width will be recycled.
#' @param position_sort Logical. If TRUE, position of sites will be sorted. Default: TRUE.
#' @param relative Logical. If FALSE, ylim of weblogo will be fixed.
#' @param mark Logical. If TRUE, Mark on Weblogo. Default: FALSE.
#' @param mark_data Data.frame or List. Only available when mark = TRUE.
#' @param title Logical. If TRUE, title each weblogo.
#' @param title.color The color of title of each weblogo
#' @param text_angle Numeric. Site number angle.Default:0.
#' @param site_size Numeric. Site size.Default:10.
#' @importFrom BiocGenerics width
#' @importFrom cowplot plot_grid
#' @import ggplot2
#' @import ggseqlogo
#' @importFrom grDevices colorRampPalette
#' @return weblogo for specific site of each family.
#' @export
#-------------------------------------------------------------------------------
##position site weblogo

position_site_weblogo <- function(seqs = NULL,alignment_file = NULL,group = NULL,type = "AA",
         max_seq_length = 30,ncol = 1,range = "conservation",
         stack_width_threshold = 0.2,Specific_site = NULL,legend = F,
         stack_width = NULL,position_sort = T,relative = F,mark = F,
         mark_data = NULL,title = T,title.color = NULL,text_angle = 0,
         site_size = 10){
  if (is.null(seqs)){
    my_seq_alignment <- switch(type,DNA = Biostrings::readDNAStringSet(alignment_file),
                               RNA = Biostrings::readRNAStringSet(alignment_file),
                               AA = Biostrings::readAAStringSet(alignment_file))
  } else if (is.null(alignment_file)){
    if ('character' %in% class(seqs)){
      my_seq_alignment <- switch(type,DNA = Biostrings::DNAStringSet(seqs),
                                 RNA = Biostrings::RNAStringSet(seqs),
                                 AA = Biostrings::AAStringSet(seqs))
    } else if(any(c('AAStringSet','DNAStringSet','RNAStringSet') %in% class(seqs))){
      my_seq_alignment <- seqs
    } else {stop('Wrong imput sequence type')}
  } else {stop("No imput sequences")}
  if (is.null(group)){
    my_seq_conservation <- get_seq_conservation(seqs = my_seq_alignment,group = NULL,type = type)
    seq_list <- list(my_seq_alignment)
  } else {
    my_seq_conservation <- get_seq_conservation(seqs = my_seq_alignment,group = group,type = type)
    seq_list <- split(my_seq_alignment,group)
    group_name <- names(seq_list)
  }
  if (is.null(title.color)){
    title.color <- rep("black",length(seq_list))
  } else {
    title.color <- rep(title.color,length(seq_list))
    length(title.color) <- length(seq_list)
  }
  my_seq_conservation <- switch(type,DNA = my_seq_conservation / 2,
                                RNA = my_seq_conservation / 2,
                                AA = my_seq_conservation / log2(20))
  for (i in 1:length(seq_list)){
    if ('matrix' %in% class(my_seq_conservation)){
      seq_conservation <- my_seq_conservation[i,]
    } else if (any(c('numeric','integer') %in% class(my_seq_conservation))){
      seq_conservation <- my_seq_conservation
    }
    seq_list_mat <- as.matrix(seq_list[[i]])
    if ('list' %in% class(range)){
      my_range <- range[[i]]
      seq_list_mat_avail <- seq_list_mat[,range[[i]]]
    } else if ('data.frame' %in% class(range)){
      my_range <- my_range[,i]
      seq_list_mat_avail <- seq_list_mat[,range[,i]]
    } else if ('integer' %in% class(range) | 'numeric' %in% class(range)){
      seq_list_mat_avail <- seq_list_mat[,range]
      my_range <- range
    } else if (is.null(range)){
      my_range <- 1:max(max_seq_length,Biostrings::width(my_seq_alignment[1]))
    } else if (range == "site"){
      if (is.null(group_name)){
        stop("Unable to automatically generate range, because no group_name")
      }
      if ('data.frame' %in% class(Specific_site)){
        my_Specific_site <- Specific_site[,c("cluster","gene")]
        my_Specific_site_group <- my_Specific_site[my_Specific_site$cluster == group_name[i],]
        if (any(c('numeric','integer') %in% class(my_Specific_site_group$gene))){
          my_range <- as.numeric(my_Specific_site_group$gene)
        } else if ('character' %in% class(my_Specific_site_group$gene)){
          my_range <- as.numeric(gsub(".* ","",my_Specific_site_group$gene))
        } else {
          stop("Unable to automatically generate range, because no group site")
        }
      }
    } else if (range == "conservation"){
      my_range <- order(seq_conservation,decreasing = T)
      if (length(my_range) > max_seq_length){
        length(my_range) <- max_seq_length
      }
    } else {
      stop("Please imput correct range")
    }
    if (position_sort){
      my_range <- sort(my_range)
    }
    seq_list_mat_avail <- seq_list_mat[,my_range]
    if (is.null(stack_width)){
      my_stack_width <- seq_conservation[my_range]
    } else if('integer' %in% class(stack_width) | 'numeric' %in% class(stack_width)){
      if (length(stack_width) == 1){
        my_stack_width <- stack_width
      } else if (length(stack_width) > 1){
        my_stack_width <- rep(stack_width,length(my_range))
        length(my_stack_width) <- length(my_range)
      }
    } else if ('list' %in% class(stack_width)){
      my_stack_width <- stack_width[[i]]
      length(my_stack_width) <- length(my_range)
    } else if ('data.frame' %in% class(stack_width)){
      my_stack_width <- stack_width[,i]
      my_stack_width <- rep(my_stack_width,length(my_range))
      length(my_stack_width) <- length(my_range)
    } else {
      stop("Wrong stack_width class!")
    }
    my_stack_width[is.na(my_stack_width)] <- 0.01
    my_stack_width[my_stack_width <= stack_width_threshold] <- 0.01
    my_plot_seq <- apply(seq_list_mat_avail,1,function(x)paste0(x,collapse = ""))
    p <- ggplot() + my_geom_logo(my_plot_seq,stack_width = my_stack_width) +
      ggseqlogo::theme_logo() +
      scale_x_continuous(labels = my_range,breaks = 1:length(my_range),
                         expand = c(0,0),limits = c(-1,max_seq_length + 1))
    if (!relative){
      if (type %in% c("DNA","RNA")){
        p <- p + scale_y_continuous(limits = c(-0.1,2.1))
      } else if (type %in% c("AA")){
        p <- p + scale_y_continuous(limits = c(-0.1,log2(20)+0.1))
      }
    }
    if (!legend){
      p <- p + theme(legend.position = "none")
    }
    if (mark){
      if ("list" %in% class(mark_data)){
        my_mark_data <- mark_data[[i]]
        if ('data.frame' %in% class(my_mark_data)){
          if (type %in% c("DNA","RNA")){
            ymin = -0.05
            ymax = 2 + 0.05
          } else if (type %in% c("AA")){
            ymin = -0.05
            ymax = log2(20) + 0.05
          }
          for (j in 1:nrow(my_mark_data)){
            #print(my_mark_data$xmin[j])
            p <- p + annotate(geom = as.character(my_mark_data$geom[j]),
                              xmin = my_mark_data$xmin[j],
                              xmax = my_mark_data$xmax[j],
                              ymin = ymin,
                              ymax = ymax,
                              col = my_mark_data$col[j],
                              fill = my_mark_data$fill[j],
                              alpha = my_mark_data$alpha[j])
          }
        }
      }else if ("data.frame" %in% class(mark_data)){
        if (type %in% c("DNA","RNA")){
          ymin = -0.05
          ymax = 2 + 0.05
        } else if (type %in% c("AA")){
          ymin = -0.05
          ymax = log2(20) + 0.05
        }
        my_mark_data <- mark_data
        for (j in 1:nrow(my_mark_data)){
          p <- p + annotate(geom = as.character(my_mark_data$geom[j]),
                            xmin = my_mark_data$xmin[j],
                            xmax = my_mark_data$xmax[j],
                            ymin = ymin,
                            ymax = ymax,
                            col = my_mark_data$col[j],
                            fill = my_mark_data$fill[j],
                            alpha = my_mark_data$alpha[j])
        }
      } else {
        warning("No mark data!")
      }
    }
    if (title){
      p <- p + labs(title = names(seq_list)[i]) +
        theme(plot.title = element_text(hjust = 0.5,face = "bold",color = title.color[i]))
    }
    p <- p + theme(plot.margin = margin(0,0,0,0),

                   axis.text.x = element_text(angle = text_angle,size = site_size))
    assign(paste0("p",i),value = p)
  }
  my_pic <- paste0("cowplot::plot_grid(",paste(paste0("p",1:length(seq_list)),collapse = ","),',ncol = ncol, align = "v")')
  final_plot <- eval(parse(text = my_pic))
  final_plot
}

#===============================================================================
##This function is in ggseqlogo package
my_geom_logo <- function(data = NULL, method='bits', seq_type='auto', namespace=NULL,
                      font='roboto_medium', stack_width=0.95, rev_stack_order=F, col_scheme = 'auto',
                      low_col='black', high_col='yellow', na_col='grey20',
                      plot=T, ...) {

  if(stack_width > 1 | stack_width < 0) stop('"stack_width" must be between 0 and 1')
  if(is.null(data)) stop('Missing "data" parameter!')
  if(!is.null(namespace)) seq_type = 'other'

  # Validate method
  all_methods = c('bits', 'probability','custom')#, 'tsl')
  pind = pmatch(method, all_methods)
  method = all_methods[pind]
  if(is.na(method)) stop("method must be one of 'bits' or 'probability', or 'custom'")

  # Convert character seqs to list
  if(is.character(data) | is.matrix(data)) data = list("1"=data)

  if(is.list(data)){
    # Set names for list if they dont exist
    if(is.null(names(data))) names(data) = seq_along(data)

    lvls = names(data)

    # We have list of sequences - loop and rbind
    data_sp = lapply(names(data), function(n){
      curr_seqs = data[[n]]
      my_logo_data(seqs = curr_seqs, method = method, stack_width = stack_width,
                rev_stack_order = rev_stack_order, seq_group = n, seq_type = seq_type,
                font = font, namespace=namespace)
    })
    data = do.call(rbind, data_sp)
    # Set factor for order of facet
    data$seq_group = factor(data$seq_group, levels = lvls)
  }

  if(!plot) return(data)

  # Get sequence type
  seq_type = attr(data, 'seq_type')
  cs = get_col_scheme( col_scheme, seq_type )

  legend_title = attr(cs, 'cs_label')

  data = merge(data, cs, by='letter', all.x=T)

  # Make sure you retain order after merge
  data = data[order(data$order),]

  # Do we have a gradient colscale
  colscale_gradient = is.numeric( cs$group )

  colscale_opts = NULL
  if(colscale_gradient){
    # Set gradient colours
    colscale_opts = scale_fill_gradient(name=legend_title, low = low_col,
                                        high = high_col, na.value = na_col)
  }else{
    # Make group -> colour map
    tmp = cs[!duplicated(cs$group) & !is.na(cs$group),]
    col_map = unlist( split(tmp$col, tmp$group) )

    # Set colour scale options
    colscale_opts = scale_fill_manual(values=col_map, name=legend_title, na.value=na_col)
  }

  # If letters and group are the same, don't draw legend
  guides_opts = NULL
  if(identical(cs$letter, cs$group)) guides_opts = guides(fill=F)

  y_lim = NULL
  extra_opts = NULL
  if(method == 'tsl'){
    y_lab = 'Depleted    Enriched'
    tmp = max(abs(data$y))
    #y_lim = c(-tmp, tmp)
    row_a = row_b = data[1,]
    row_a$y = -tmp
    row_b$y = tmp
    data = rbind(data, row_a, row_b)
    data$facet = factor(data$y > 0, c(T, F), c('Enriched', 'Depleted'))
    extra_opts = NULL#facet_grid(facet~., scales='free')
  }else if(method == 'custom'){
    y_lab = ''
  }else{
    y_lab = method
    substr(y_lab, 1, 1) = toupper(substr(y_lab, 1, 1))
  }

  # Group data
  data$group_by = with(data, interaction(seq_group, letter, position))

  data$x = data$x
  # Create layer
  logo_layer = layer(
    stat = 'identity', data = data,
    mapping = aes_string(x = 'x', y = 'y', fill='group', group='group_by'),
    geom = 'polygon',
    position = 'identity', show.legend = NA, inherit.aes = F,
    params = list(na.rm = T, ...)
  )


  breaks_fun = function(lim){
    # account for multiplicatuce expansion factor of 0.05
    1: floor( lim[2] / 1.05 )
  }

  # Expand 0.05 addidtive
  list(logo_layer, scale_x_continuous(breaks = breaks_fun, labels = identity),
       ylab(y_lab), xlab(''), colscale_opts, guides_opts, coord_cartesian(ylim=y_lim),
       extra_opts)
}
#===============================================================================
##This function is in ggseqlogo package
##We modified this function to allow different stacks in the weblogo.
##In some cases, the width of the stack is proportional to the
##fraction of valid symbols in that position. (Positions with many gaps have thin stacks.)
my_logo_data <- function( seqs, method='bits', stack_width=0.95,
                       rev_stack_order=F, font, seq_group=1,
                       seq_type = 'auto', namespace=NULL ){

  seq_width <- nchar(seqs)[1]
  # Get font
  font_df = get_font(font)

  # TODO
  # hh = twosamplelogo_method(seqs, seqs_bg, pval_thresh=0.05, seq_type = seq_type, namespace = namespace)

  # Generate heights based on method
  if(method == 'bits'){
    hh = bits_method(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace)
  }else if(method == 'probability'){
    hh = probability_method(seqs, decreasing = rev_stack_order, seq_type = seq_type, namespace = namespace)
  }else if(method == 'custom'){
    if(seq_type == 'auto') seq_type = guessSeqType(rownames(seqs))
    hh = matrix_to_heights(seqs, seq_type, decreasing = rev_stack_order)
  }else{
    stop('Invalid method!')
  }

  # Merge font df and heights
  ff = merge(font_df, hh, by = 'letter')
  # Scale x and ys to new positions
  if (length(stack_width) == 1){
    x_pad = stack_width/2
    ff$x = newRange(ff$x, ff$position - x_pad, ff$position + x_pad)
    ff$y = newRange(ff$y, ff$y0, ff$y1)
  } else if (length(stack_width) == seq_width){
    ff <- ff[order(ff$position),]
    ff$position.factor <- as.factor(ff$position)
    levels(ff$position.factor) <- 1:seq_width
    my_position <- table(ff$position.factor)
    x_pad = rep(stack_width/2,c(my_position))
    ff$x = newRange(ff$x, ff$position - x_pad, ff$position + x_pad)
    ff$y = newRange(ff$y, ff$y0, ff$y1)
    ff <- ff[order(ff$letter),]
  }
  # Rename columns
  ff = as.data.frame(ff)[,c('x', 'y', 'letter', 'position', 'order')]
  ff$seq_group = seq_group

  # Set sequence type as attribute, to be used downstream
  attr(ff, 'seq_type') = attr(hh, 'seq_type')

  # Return data table
  ff
}
#===============================================================================
#' MDST Network visualization with igraph
#'
#' @param MDST MDST Network.
#' @param node.color The fill color of the vertex.
#' @param group The vertex labels.
#' @param weight Logical. Whether to determine the width of the edges. according to the weight.Default:TRUE
#' @param node.label.color The color of the node labels.
#' @param edge.color The color of the edges.
#' @param arrow.size Size of arrows.
#' @param node.size Size of nodes.
#' @param frame.color The color of the frame of the vertices.
#' @param edge.label The edge labels.
#' @param edge.width Width of each edge.
#' @param edge.label.color The color of the edge labels.
#' @param edge.label.loc The location of edge labels, if edge.label.loc = 0, edge label locate at the same position of start node,
#'                       if edge.label.loc = 1, edge label locate at the same position of end node,
#'                       if edge.label.loc = 0.5, edge label locate at the middle of start and end node.Default:0.5.
#' @param edge.label.adjust edge label will be adjusted in the x-axis direction.Default:0.1.
#' @param edge.label.cex The font size for the edge labels.Default:0.8.
#' @param arrow.width The width of the arrows. Default:layout_as_tree
#' @param my_layout Either a function or a numeric matrix. It specifies how the vertices will be placed on the plot.
#' @param rescale Logical. Whether to rescale the coordinates to the [-1,1]x[-1,1](x[-1,1]) interval.
#'                This parameter is not implemented for tkplot.Default:FALSE, the layout will NOT be rescaled.
#' @param node.self Logical. Whether to draw an arrow pointing to the vertice itself.Default:FALSE.
#' @param ... Other parameter of plot.igraph.
#' @importFrom igraph add_edges
#' @importFrom igraph E<-
#' @importFrom igraph V<-
#' @importFrom igraph E
#' @importFrom igraph V
#' @importFrom igraph make_empty_graph
#' @importFrom igraph layout_as_tree
#' @importFrom grDevices colorRampPalette
#' @seealso \code{\link{plot.igraph},\link{igraph.plotting}}
#' @export
#'
##MDST Network visualization with igraph
#-------------------------------------------------------------------------------
MDST_plot <- function(MDST,node.color = NULL,group = NULL,weight = T,
                      node.label.color = "black",edge.color = "black",
                      arrow.size = 0.4,node.size = NULL,frame.color = NULL,
                      edge.label = NULL,edge.width = 2,edge.label.color = "black",
                      edge.label.loc = 0.5,edge.label.adjust = 0.1,
                      edge.label.cex = 0.8,arrow.width = 1.5,
                      my_layout = layout_as_tree,rescale = F,
                      node.self = F,...
){
  options(warn = 1)
  numCluster <- max(MDST[,1:2])
  if (!node.self){MDST <- MDST[MDST[,1] != MDST[,2],]}
  net <- add_edges(make_empty_graph(numCluster), t(MDST[,1:2]), weight=MDST[,3])
  if (!is.null(node.color)){V(net)$color <- node.color}
  if (!is.null(frame.color)){V(net)$frame.color <- frame.color
  } else {V(net)$frame.color <- NA}
  if (!is.null(node.label.color)){V(net)$label.color <- node.label.color}
  if (!is.null(edge.label)){
    E(net)$label <- edge.label
    E(net)$label.color <- edge.label.color
    E(net)$label.cex <- edge.label.cex
  }

  if (!is.null(group)){V(net)$label <- group}
  if (is.null(node.size)){V(net)$size <- 30
  } else {
    V(net)$size <- node.size / max(node.size) * 50
  }
  if (weight){E(net)$width <- (E(net)$weight / max(E(net)$weight)) * 5
  } else {E(net)$width <- edge.width}
  E(net)$arrow.size <- arrow.size
  E(net)$arrow.width <- arrow.width
  E(net)$color <- edge.color
  if ('function' %in% class(my_layout)){
    my_layout <- my_layout(net)
  } else if (any(c('data.frame',"matrix") %in% class(my_layout))){
    my_layout <- as.matrix(my_layout[,1:2])
  }
  #print(my_layout)
  #graph_attr(net, "layout") <- my_layout
  if (rescale){
    plot(net,layout = my_layout,rescale = rescale,...)
  } else {
    if (!is.null(edge.label.loc)){
      edge.label.pos.x <-  edge.label.loc * (-my_layout[,1][MDST[,1]] + my_layout[,1][MDST[,2]]) +
        my_layout[,1][MDST[,1]]
      edge.label.pos.y <-  edge.label.loc * (-my_layout[,2][MDST[,1]] + my_layout[,2][MDST[,2]]) +
        my_layout[,2][MDST[,1]]
      if (!is.null(edge.label.adjust)){
        for (i in 1:nrow(MDST)){
          if (my_layout[,1][MDST[,1]][i] < my_layout[,1][MDST[,2]][i]){
            edge.label.pos.x[i] <- edge.label.pos.x[i] + edge.label.adjust
          } else if (my_layout[,1][MDST[,1]][i] > my_layout[,1][MDST[,2]][i]){
            edge.label.pos.x[i] <- edge.label.pos.x[i] - edge.label.adjust
          } else if (my_layout[,1][MDST[,1]][i] == my_layout[,1][MDST[,2]][i]){
            if (my_layout[,1][MDST[,1]][i] <= mean(my_layout[,1])){
              edge.label.pos.x[i] <- edge.label.pos.x[i] - edge.label.adjust
            } else {
              edge.label.pos.x[i] <- edge.label.pos.x[i] + edge.label.adjust
            }
          }
        }
      }
      #print(edge.label.pos.x);print(edge.label.pos.y)
      E(net)$label.x <- edge.label.pos.x
      E(net)$label.y <- edge.label.pos.y
    }
    #print(my_layout)
    #print(net)
    plot(net,layout = my_layout,rescale = rescale,xlim = range(my_layout[,1]),ylim = range(my_layout[,2]),...)
  }
}
#===============================================================================
##This function is in ggseqlogo package
get_col_scheme = function(col_scheme, seq_type='auto'){

  # Check if user-defined color scheme
  if(is.data.frame(col_scheme)){
    if(!'ggseqlogo_cs' %in% class(col_scheme))
      stop('Colour scheme must be generated using "make_col_scheme" function')
    return(col_scheme)
  }

  # Get ambigious colour scheme
  col_scheme = match.arg(col_scheme, list_col_schemes(F))

  # Get default color scheme for sequence type
  if(col_scheme == 'auto'){
    if(seq_type == 'auto') stop('"col_scheme" and "seq_type" cannot both be "auto"')

    col_scheme = switch(tolower(seq_type), aa = 'chemistry',
                        dna = 'nucleotide', rna = 'nucleotide',
                        other='nucleotide')

  }


  # Pick from default color schemes
  cs = switch(col_scheme,
              # Color scheme based on chemistry of amino acids
              chemistry2 = data.frame(
                letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
                group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
                col = c(rep('#058644', 5), rep('#720091', 2), rep('#0046C5', 3), rep('#C5003E', 2), rep('#2E2E2E', 8)),
                stringsAsFactors = F
              ),

              # Color scheme based on chemistry of amino acids
              chemistry = data.frame(
                letter = c('G', 'S', 'T', 'Y', 'C', 'N', 'Q', 'K', 'R', 'H', 'D', 'E', 'P', 'A', 'W', 'F', 'L', 'I', 'M', 'V'),
                group = c(rep('Polar', 5), rep('Neutral', 2), rep('Basic', 3), rep('Acidic', 2), rep('Hydrophobic', 8)),
                col = c(rep('#109648', 5), rep('#5E239D', 2), rep('#255C99', 3), rep('#D62839', 2), rep('#221E22', 8)),
                stringsAsFactors = F
              ),

              # Hydrophobicity index (PMID: 7108955) from -4.5 to 4.5
              hydrophobicity = data.frame(
                letter = c('I', 'V', 'L', 'F', 'C', 'M', 'A', 'G', 'T', 'W',
                           'S', 'Y', 'P', 'H', 'D', 'E', 'N', 'Q', 'K', 'R'),
                group = c(4.5, 4.2, 3.8, 2.8, 2.5, 1.9, 1.8, -0.4, -0.7, -0.9, -0.8,
                          -1.3, -1.6, -3.2, -3.5, -3.5, -3.5, -3.5, -3.9, -4.5),
                stringsAsFactors=F
              ),

              # Colour based on nucleotide
              nucleotide2 = data.frame(
                letter = c('A', 'C', 'G', 'T', 'U'),
                col = c('darkgreen', 'blue', 'orange', 'red', 'red'),
                stringsAsFactors = F
              ),

              #alt red BA1200
              nucleotide = data.frame(
                letter = c('A', 'C', 'G', 'T', 'U'),
                col = c('#109648', '#255C99', '#F7B32B', '#D62839', '#D62839'),
                stringsAsFactors = F
              ),

              base_pairing = data.frame(
                letter = c('A', 'T', 'U', 'G', 'C'),
                group = c(rep('Weak bonds', 3), rep('Strong bonds', 2)),
                col = c(rep('darkorange', 3), rep('blue', 2)),
                stringsAsFactors = F
              ),

              # ClustalX color scheme:
              # http://www.jalview.org/help/html/colourSchemes/clustal.html
              clustalx = data.frame(
                letter = c('W', 'L', 'V', 'I', 'M', 'F', 'A', 'R', 'K', 'T', 'S', 'N', 'Q', 'D', 'E', 'H', 'Y', 'C', 'G', 'P'),
                col = c(rep('#197FE5', 7), rep('#E53319', 2), rep('#19CC19', 4), rep('#CC4CCC', 2),
                        rep('#19B2B2', 2), '#E57F7F', '#E5994C', '#B0B000'),
                stringsAsFactors = F
              ),

              # Taylor color scheme (PMID: 9342138)
              taylor = data.frame(
                letter = c('D','S','T','G','P','C','A','V','I','L','M','F','Y','W','H','R','K','N','Q','E'),
                col = c('#FF0000','#FF3300','#FF6600','#FF9900','#FFCC00','#FFFF00','#CCFF00','#99FF00',
                        '#66FF00','#33FF00','#00FF00','#00FF66','#00FFCC','#00CCFF','#0066FF','#0000FF',
                        '#6600FF','#CC00FF','#FF00CC','#FF0066'),
                stringsAsFactors = F
              )
  )

  if(!'group' %in% names(cs)) cs$group = cs$letter

  # Set attributes
  attr(cs, 'cs_label') = col_scheme
  class(cs) = c('data.frame','ggseqlogo_cs')

  return(cs)
}
#===============================================================================
##This function is in ggseqlogo package
bits_method <- function(seqs, decreasing, ...){
  # Make PFM
  pfm = makePFM(seqs, ...)

  # Get ic
  ic = attr(pfm, 'bits')
  if(all(ic == 0)){
    warning('All positions have zero information content perhaps due to too few input sequences. Setting all information content to 2.')
    ic = (ic * 0)+2
  }
  heights = t(t(pfm) * ic)

  seq_type = attr(pfm, 'seq_type')
  matrix_to_heights(heights, seq_type, decreasing)
}
#===============================================================================
##This function is in ggseqlogo package
probability_method <- function(seqs, decreasing, ...){
  # Make PFM
  pfm = makePFM(seqs, ...)
  seq_type = attr(pfm, 'seq_type')
  matrix_to_heights(pfm, seq_type, decreasing)
}
#===============================================================================
##This function is in ggseqlogo package
matrix_to_heights <- function(mat, seq_type, decreasing=T){

  mat[is.infinite(mat)] = 0

  if(any(duplicated(rownames(mat)))) stop('Matrix input must have unique row names')

  dat = lapply(1:ncol(mat), function(i){
    vals = mat[,i]

    pos = sort( vals[vals >= 0], decreasing = decreasing)
    neg = sort(vals[vals < 0], decreasing = !decreasing)
    #vals = sort(vals, decreasing = T)
    cs_pos = cumsum( pos )
    cs_neg = cumsum( neg )

    df_pos = df_neg = NULL

    if(length(pos) > 0)
      df_pos = data.frame(letter=names(pos), position=i,  y0=c(0, cs_pos[-length(cs_pos)]),
                          y1=cs_pos, stringsAsFactors = F)

    if(length(neg) > 0)
      df_neg = data.frame(letter=names(neg), position=i, y0=cs_neg, y1=c(0, cs_neg[-length(cs_neg)]),
                          stringsAsFactors = F)

    rbind(df_pos, df_neg)
  })

  dat = do.call(rbind, dat)

  # Adjust y spacing
  space_factor = 0.004
  y_pad = max(dat$y1) * space_factor
  dat$y0 = dat$y0 + y_pad
  dat = subset(dat, dat$y1 > dat$y0)

  # Dummy points to make sure full plot is drawn
  # Make sure position 1 and n have a dummy empty letter missing
  dummy = data.frame(letter=dat$letter[1], position=NA, y0=0, y1=0)

  # Missing first position
  if(dat$position[1] != 1){
    dummy$position = 1
    dat = rbind( dummy, dat )
  }

  # Missing last position
  if(dat$position[nrow(dat)] != ncol(mat)){
    dummy$position = ncol(mat)
    dat = rbind( dat, dummy )
  }

  rownames(dat) = NULL

  attr(dat, 'seq_type') = seq_type

  dat
}
#===============================================================================
##This function is in ggseqlogo package
get_font <- function(font){

  GGSEQLOGO_FONT_BASE = getOption('GGSEQLOGO_FONT_BASE')
  if(is.null(GGSEQLOGO_FONT_BASE)){
    GGSEQLOGO_FONT_BASE=system.file("extdata", "", package = "ggseqlogo")
    options(GGSEQLOGO_FONT_BASE=GGSEQLOGO_FONT_BASE)
  }

  #all_fonts = c('sf_bold', 'sf_regular', 'ms_bold', 'ms_regular', 'xkcd_regular')
  font = match.arg(tolower(font), ggseqlogo::list_fonts(F))
  font_filename = paste0(font, '.font')
  font_obj_name = sprintf('.ggseqlogo_font_%s', font)

  font_obj = getOption(font_obj_name)
  if(is.null(font_obj)){
    # Not loaded into global env yet - load it into options
    font_path = file.path(GGSEQLOGO_FONT_BASE, font_filename)
    font_obj_list = list( tmp=readRDS(font_path) )
    names(font_obj_list) = font_obj_name
    options(font_obj_list)
    font_obj = font_obj_list[[1]]
  }

  # Return font data
  font_obj
}
#===============================================================================
##This function is in ggseqlogo package
guessSeqType <- function(sp){
  # Ensure we have something
  if(length( intersect(sp, c(.AA_NAMESPACE(), .DNA_NAMESPACE(),.RNA_NAMESPACE())) ) == 0)
    stop('Could not get guess seq_type. Please explicitly define sequence type or use "other" with custom namespaces.')

  dat = setdiff(intersect(sp, .AA_NAMESPACE()), c(.DNA_NAMESPACE(),.RNA_NAMESPACE()))
  if(length(dat) > 0){
    return('AA')
  }else if('U' %in% sp){
    return('RNA')
  }
  return('DNA')
}
#===============================================================================
##This function is in ggseqlogo package
newRange <- function(old_vals, new_min=0, new_max=1){
  old_min = min(old_vals)
  old_max = max(old_vals)

  new_vals = (((old_vals - old_min) * (new_max - new_min)) / (old_max - old_min)) + new_min
  new_vals
}

#===============================================================================
