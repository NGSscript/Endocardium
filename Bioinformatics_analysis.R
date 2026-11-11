# load R packages --------------------------------------------------------------------
library(openxlsx)
library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(ggrepel)
library(tibble)
library(ComplexHeatmap)
library(patchwork)
library(pheatmap)
# GpTn-GpTp Differential expression analysis --------------------------------------------------------------------
GpTn_sample<-data.frame(Sample=c("GpTn41", "GpTn42", "GpTn61", "GpTn62", "GpTn63",
                                 "GpTp41", "GpTp42", "GpTp61", "GpTp62", "GpTp63"),
                        Groups=c(rep("GpTn",5),rep("GpTp",5)))
GpTn_counts <- read.table('RNAseq/GpTn_GpTp.count.txt',header = T, sep = "\t" ,row.names = 1)
gene_id_name_GpTn<-GpTn_counts[,c("gene_id","gene_name")]
GpTn_counts <- GpTn_counts[,GpTn_sample$Sample]
rownames(GpTn_sample)<-GpTn_sample$Sample
GpTn_sample$Groups=as.factor(GpTn_sample$Groups)

GpTn_dds <- DESeqDataSetFromMatrix(countData = GpTn_counts,
                                   colData=GpTn_sample, design = ~Groups) 
keep <- rowSums(counts(GpTn_dds))>=10 
GpTn_dds <- GpTn_dds[keep,]
dim(GpTn_dds)
GpTn_dds$Groups<-relevel(GpTn_dds$Groups, ref = 'GpTn')
GpTn_ddsDE <- DESeq(GpTn_dds)
GpTn_vsdata <- vst(GpTn_ddsDE, blind = FALSE)
res<- results(GpTn_ddsDE, alpha = 0.05)
resOrdered <- res[order(res$padj),]
resOrdered <- as.data.frame(resOrdered)
resOrdered$symbol<- mapIds(org.Mm.eg.db, keys=rownames(resOrdered),keytype = "ENSEMBL", column = "SYMBOL")
resOrdered<-na.omit(resOrdered)
resOrdered[["DifferentialExpression"]] <- ifelse(resOrdered$padj < 0.05 &abs(resOrdered$log2FoldChange) > 1,
                                                ifelse(resOrdered$log2FoldChange > 0, "Up-regulated", "Down-regulated"), "Non-significant")
table(resOrdered$DifferentialExpression)
up_gene<-resOrdered%>%filter(DifferentialExpression %in% c("Up-regulated"))%>%pull(symbol)
down_gene<-resOrdered%>%filter(DifferentialExpression %in% c("Down-regulated"))%>%pull(symbol)
up_GO <- enrichGO(gene = up_gene, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
down_GO <- enrichGO(gene = down_gene, OrgDb = org.Mm.eg.db,keyType = 'SYMBOL', ont = 'BP',pAdjustMethod = 'fdr', pvalueCutoff = 0.05,qvalueCutoff = 0.05,readable = T)
up_GO1<-as.data.frame(up_GO)
down_GO1<-as.data.frame(down_GO)


write.xlsx(resOrdered, file = 'RNAseq/result/DEG.xlsx')
write.xlsx(up_GO, file = 'RNAseq/result/downDEG_GOenrichment.xlsx')
write.xlsx(down_GO, file = 'RNAseq/result/upDEG_GOenrichment.xlsx')

# combined_plot  --------------------------------------------------------------

terms_up <- c("mesenchyme development",
"atrial septum development",
"endocardial cushion development",
"mesenchymal cell differentiation",
"secondary heart field specification",
"Wnt signaling pathway",
"negative regulation of cell adhesion",
"epithelial to mesenchymal transition",
"transforming growth factor beta receptor signaling pathway")

terms_down<-c( "blood coagulation",
"regulation of hemopoiesis",
"endothelium development",
"positive regulation of cell-cell adhesion",
"establishment of endothelial barrier")

df_up<-up_GO1%>%filter(Description %in% terms_up)
df_down<-down_GO1%>%filter(Description %in% terms_down)
df_up$Cluster<-"UP"
df_down$Cluster<-"UP"

df_up <- df_up %>%
  mutate(
    GeneRatio_num = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    Cluster = factor(Cluster)
  ) %>%
  group_by(Cluster) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()
df_up$log10_P <- -log10(df_up$pvalue)
df_up <- df_up %>%
  arrange(log10_P) %>%
  mutate(Term_Ordered = factor(Description, levels = Description))


df_down <- df_down %>%
  mutate(
    GeneRatio_num = sapply(strsplit(as.character(GeneRatio), "/"), function(x) as.numeric(x[1]) / as.numeric(x[2])),
    Cluster = factor(Cluster)
  ) %>%
  group_by(Cluster) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()
df_down$log10_P <- -log10(df_down$pvalue)
df_down <- df_down %>%
  arrange(log10_P) %>%
  mutate(Term_Ordered = factor(Description, levels = Description))



plot_up <- ggplot(df_up, aes(x = log10_P, y = Term_Ordered)) +
  geom_col(aes(fill = log10_P), width = 0.75) +
  geom_col(color = "black", fill = NA, width = 0.75) +
  scale_fill_gradientn(
    colours = c("#FFFBC8", "#FFC864", "#992200", "#260901"),
    name = "-log10(P)",
    limits = c(0, 12),oob = scales::squish
  ) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 3))+
  labs(
    x = "-log10(P-value)",
    y = "",
    title = "UP-regulated genes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 15, hjust = 1, colour = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
plot_down <- ggplot(df_down, aes(x = log10_P, y = Term_Ordered)) +
  geom_col(aes(fill = log10_P), width = 0.75) +
  geom_col(color = "black", fill = NA, width = 0.75) +
  scale_fill_gradientn(
    colours = c("#A3F0BE", "#83D6C5", "#020B6E", "#072A45"),
    limits = c(0, 12),oob = scales::squish, 
    name = "-log10(P)"
  ) +
  scale_x_continuous(limits = c(0, 12), breaks = seq(0, 12, 3))+
  labs(
    x = "-log10(P-value)",
    y = "",
    title = "DOWN-regulated genes"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 15, hjust = 1, colour = "black"), 
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )
combined_plot <- plot_up / plot_down +
  plot_layout(
    heights = c(nrow(df_up), nrow(df_down)), 
    guides = "collect"
  ) +
  plot_annotation(theme = theme(plot.margin = margin(5, 5, 5, 5)))
combined_plot





# volcano_plot -------------------------------------------------------------------
DEG<-read.xlsx('RNAseq/result/DEG.xlsx')
rownames(DEG)<-DEG$gene_id
volcano_plot_enhanced <- function(df, num_symbol = 5, logFC="log2FoldChange", FDR="padj", logFC_Value=1,
                                  Symbol="symbol", y_increased=10, labs="symbol",title="volcano") {
  df <- df %>% mutate(DEG = case_when(
    !!sym(logFC) > logFC_Value & !!sym(FDR) < 0.05 ~ "Up",
    abs(!!sym(logFC)) < logFC_Value | !!sym(FDR) > 0.05 ~ "None",
    !!sym(logFC) < -logFC_Value & !!sym(FDR) < 0.05 ~ "Down"
  ))  
  if (!logFC %in% colnames(df)) {
    cat(paste0("The '", logFC, "' column is not found in the dataset.\n"))
    logFC <- readline("Please enter the name of the 'logFC' column: ")
  }
  if (!FDR %in% colnames(df)) {
    cat(paste0("The '", FDR, "' column is not found in the dataset.\n"))
    FDR <- readline("Please enter the name of the 'FDR' column: ")
  }
  if (!Symbol %in% colnames(df)) {
    cat(paste0("The '", Symbol, "' column is not found in the dataset.\n"))
    Symbol <- readline("Please enter the name of the 'Symbol' column: ")
  }
  print(paste0("Column name of logFC is: ", logFC))
  print(paste0("Column name of FDR is: ", FDR))
  print(paste0("Column name of Symbol is: ", Symbol))
  upregulated_genes <- sum(df[[logFC]] > logFC_Value & df[[FDR]] < 0.05)
  downregulated_genes <- sum(df[[logFC]] < -logFC_Value & df[[FDR]] < 0.05)
  y_upper <- ceiling(max(-log10(df[[sym(FDR)]]), na.rm = TRUE))
  x_abs <- max(abs(c(
    floor(min(df[[sym(logFC)]], na.rm = TRUE)),
    ceiling(max(df[[sym(logFC)]], na.rm = TRUE))
  )))
  plot <- ggplot(data = df, aes(x = !!sym(logFC), y = -log10(!!sym(FDR)))) + 
    geom_point(aes(color = (!!sym(logFC)), size = -log10(!!sym(FDR))), alpha = 0.9) +
    geom_point(data = df %>% filter(DEG == "None"), aes(color = DEG, size = -log10(!!sym(FDR))), 
               color = "gray", alpha = 0.8) +
    scale_color_gradientn(colours = c("#39489f","#39bbec","#f9ed36","#f38466","#b81f25"), 
                          values = seq(0, 1, 0.2)) +
    scale_size_continuous(range = c(1, 6)) +
    geom_label_repel(data = df %>%
                       tidyr::drop_na() %>%
                       dplyr::filter(!!sym(Symbol) %in% c("Nkx2-5")),
                     aes(x = !!sym(logFC), y = -log10(!!sym(FDR)), label = !!sym(Symbol)),
                     size = 5,color = "#790000",
                     box.padding = unit(0.1, "lines"),  
                     point.padding = unit(0.6, "lines"), 
                     segment.color = "black", nudge_x = 1.5, nudge_y = 6,segment.curvature = 0.1, 
                     fontface="bold", family = "serif", 
                     max.overlaps = 10000)+
    geom_label_repel(data = df %>%tidyr::drop_na() %>%plyr::filter(!!sym(Symbol) %in% c("Sox9")),box.padding = unit(0.1, "lines"),  point.padding = unit(0.6, "lines"), segment.color = "black", nudge_x = 1.6, nudge_y = 3,segment.curvature = 0.1,fontface="bold", family = "serif",max.overlaps = 10000)+
    geom_label_repel(data = df %>%tidyr::drop_na() %>%dplyr::filter(!!sym(Symbol) %in% c("Twist1")),aes(x = !!sym(logFC), y = -log10(!!sym(FDR)), label = !!sym(Symbol)),size = 5,color = "#790000",box.padding = unit(0.1, "lines"),  point.padding = unit(0.6, "lines"), segment.color = "black", nudge_x = 1.8, nudge_y = -3,segment.curvature = -0.1,fontface="bold", family = "serif",max.overlaps = 10000)+
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = 4) + 
    scale_size(range = c(1,6)) + 
    xlim(-x_abs, x_abs) + 
    ylim(c(-1, y_upper + y_increased)) +
    theme_bw() + 
    theme(panel.grid = element_blank(),
          legend.background = element_rect(color = "#808080", linetype = 1),
          axis.text = element_text(size = 12, color = "#000000"),
          axis.title = element_text(size = 15),
          plot.title = element_text(hjust = 0.5, vjust = 8,size=18,face ="bold"),
          plot.margin = margin(t = 15, r = 5, b = 5, l = 5, unit = "mm"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 14)) + 
    coord_cartesian(clip = "off") +
    labs(title = title) +
    theme(text=element_text(family="serif"))+
    labs(title = title)+
    annotation_custom(grob = grid::segmentsGrob(
      y0 = unit(-10, "pt"),
      y1 = unit(-10, "pt"),
      arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "first"),
      gp = grid::gpar(lwd = 3, col = "#3288bd")), 
      xmin = -x_abs, xmax = -1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Down",
        gp = grid::gpar(col = "#3288bd",fontfamily = "serif",fontsize=15)),
      xmin = -x_abs, xmax = -1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::segmentsGrob(
        y0 = unit(-10, "pt"),
        y1 = unit(-10, "pt"),
        arrow = arrow(angle = 45, length = unit(.2, "cm"), ends = "last"),
        gp = grid::gpar(lwd = 3, col = "#d73027")), 
      xmin = x_abs, xmax = 1, ymin = y_upper + y_increased, ymax = y_upper + y_increased) +
    annotation_custom(
      grob = grid::textGrob(
        label = "Up",
        gp = grid::gpar(col = "#d73027",fontfamily  = "serif",fontsize=15)),
      xmin = x_abs, xmax = 1, ymin = y_upper + y_increased, ymax = y_upper + y_increased)
  
  
  return(plot)
}
volcano_plot <- volcano_plot_enhanced(DEG, title = "Nkx2-5+ lineage endocardium gene expression") 




# heatmap --------------------------------------------------------------------
GpTn_genes <- c("Emcn", "Cldn5", "Kdr", "Egfl7", "Robo4", "Zbtb46", "Gng11", "F8",
                "Wnt4", "Bmper", "Tgfb2", "Twist1", "Twist2", "Sox9", "Msx1",
                "Has2", "Vcan", "Postn", "Fn1", "Tnc", "Cd44", "Itgav", "Col3a1", "Loxl2")
GpTn_id<-mapIds(org.Mm.eg.db, keys=GpTn_genes,keytype = "SYMBOL", column = "ENSEMBL")
GpTn_id<-data.frame(gene_id=GpTn_id,gene_name=names(GpTn_id))
Gptn_vst_mat <- assay(GpTn_vsdata)

GpTn_mat <- merge(GpTn_id,Gptn_vst_mat,by.x="gene_id",by.y="row.names")
rownames(GpTn_mat)<-GpTn_mat$gene_name
GpTn_mat<-GpTn_mat[,-c(1,2)]
GpTn_mat<-GpTn_mat[GpTn_genes,]
GpTn_mat<-na.omit(GpTn_mat)

GpTn_anno_df <- as.data.frame(colData(GpTn_vsdata))[,c("Groups"),drop=FALSE]
GpTn_anno_df$Groups<-gsub("GpTn","Nkx2-5 negative endocardium",GpTn_anno_df$Groups)
GpTn_anno_df$Groups<-gsub("GpTp","Nkx2-5 linage endocardium",GpTn_anno_df$Groups)
colnames(GpTn_anno_df)<-c("Group")
GpTn_mat<-GpTn_mat[GpTn_genes,]
GpTn_col_col <- list(
  Group = c(
    "Nkx2-5 negative endocardium" = "#2983B1",
    "Nkx2-5 linage endocardium" = "#D37538"
  )
)


heatmap<-ComplexHeatmap::pheatmap(GpTn_mat,
         scale = 'row',       
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         annotation_col = GpTn_anno_df,
         show_rownames = TRUE,
         annotation_colors = c(GpTn_col_col),
         color = colorRampPalette(c("#004488", "#4488BB", "#F4F6F5", "#CC5533", "#942003"))(100),
         border_color = NA,
         show_colnames = FALSE,
         heatmap_legend_param = list(title = "Zcore"),fontsize = 12)



