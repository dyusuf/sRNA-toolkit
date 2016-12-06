setwd("/home/dilmurat/Work/combined-experiments/")
library(plyr)
library(ggplot2)

mirpath <- read.csv("miRPath_pathway_miRNA.csv", header=T, sep=",")
mirpath <- subset(mirpath, p_value<=0.001)
for (i in names(mirpath)) {mirpath[i] <- droplevels(mirpath[i])}
mirpath$pathway <- gsub("_", " ", mirpath$pathway) 

library(RColorBrewer)
getPalette = colorRampPalette(brewer.pal(8, "Accent"))
colourCount = length(unique(mirpath$miRNA))
sum_pathway <- ddply(mirpath, c("pathway"), summarise,
                     total = sum(gene_number))
mirpath$pathway <- factor(mirpath$pathway, levels = sum_pathway[order(sum_pathway$total, decreasing=F), "pathway"])
mirpath$miRNA <- factor(mirpath$miRNA, levels=as.vector(arrange(count(as.vector(mirpath$miRNA)), desc(freq))$x))
# mirpath$miRNA <- factor(mirpath$miRNA, levels= 
#                           c('hsa-miR-192-5p',
#                            'hsa-miR-25-3p',
#                            'hsa-miR-142-5p',
#                            'hsa-miR-30e-5p',
#                            'hsa-miR-22-3p',
#                            'hsa-miR-423-5p',
#                            'hsa-miR-92a-3p',
#                            'hsa-miR-103a-3p',
#                            'hsa-miR-486-5p',
#                            'hsa-miR-182-5p',
#                            'hsa-miR-16-5p',
#                            'hsa-miR-15a-5p',
#                            'hsa-miR-451a',
#                            'hsa-miR-27b-3p',
#                            'hsa-miR-126-5p',
#                            'hsa-miR-191-5p',
#                            'hsa-miR-10a-5p',
#                            'hsa-miR-143-3p',
#                            'hsa-miR-21-5p',
#                            'hsa-miR-10b-5p'))
mirpath$plasma_abundance <- "NA"
for (i in levels(mirpath$miRNA)){mirpath[mirpath["miRNA"] == i, ]$plasma_abundance <- 
                                   sort(subset(cor_exo_plasma, grepl(strsplit(i, "-")[[1]][3], rna))$plasma_abundance, decreasing=T)[1]} 
# assuming there would not be miRNAs in mirpath$miRNA which are derived from the same locus 

mirpath_class <- read.csv("mirpath_annotation-20141127.csv", header=T, sep="\t")
mirpath_class$pathway <- gsub("_"," ", mirpath_class$pathway)
mirpath$class <- "NA"
for (i in levels(mirpath$pathway)){mirpath[mirpath["pathway"] == i, ]$class <- 
                                     as.character(subset(mirpath_class, pathway == i)$class)} 

notable_pathway <- c(
  "PI3K-Akt_signaling_pathway",                  
  "MAPK_signaling_pathway",                  
  "Ubiquitin_mediated_proteolysis",             
  "Wnt_signaling_pathway",   
  "mTOR_signaling_pathway",  
  "ErbB_signaling_pathway",
  "Focal_adhesion", 
  "Regulation_of_actin_cytoskeleton",
  "T_cell_receptor_signaling_pathway",           
  "Neurotrophin_signaling_pathway", 
  "Long-term_potentiation",  
  "Cholinergic_synapse",                        
  "Glutamatergic_synapse",    
  "Endocytosis",                             
  "Insulin_signaling_pathway"                                         
)
notable_pathway <- gsub("_", " ", notable_pathway)
mirpath$notable_pathway <- "NA"
for (i in levels(mirpath$pathway)){if (i %in% notable_pathway) 
{mirpath[mirpath["pathway"] == i, ]$notable_pathway <- i} 
                                   else
                                   {mirpath[mirpath["pathway"] == i, ]$notable_pathway <- "others"}
} 

mirpath$miRNA <- factor(gsub("hsa-", "", mirpath$miRNA))

cluster1 <- c("miR-22-3p", "miR-25-3p",
              "miR-142-5p", "miR-192-5p",
              "miR-30e-5p")
cluster2 <- c("miR-486-5p", "miR-92a-3p",
              "miR-16-5p", "miR-451a",
              "miR-182-5p", "miR-103a-3p",
              "miR-15a-5p")
cluster3 <- c("miR-191-5p", "miR-126-5p",
              "miR-10a-5p", "miR-27b-3p",
              "miR-21-5p")
mirpath$cluster <- "NA"
for (i in levels(mirpath$miRNA)){
if (i %in% cluster1) 
  {mirpath[mirpath["miRNA"] == i, ]$cluster <- "cluster 1"} 
else if (i %in% cluster2) 
  {mirpath[mirpath["miRNA"] == i, ]$cluster <- "cluster 2"} 
else if (i %in% cluster3) 
{mirpath[mirpath["miRNA"] == i, ]$cluster <- "cluster 3"} 
else
  {mirpath[mirpath["miRNA"] == i, ]$cluster <- "others"}
} 

mirpath_summery <- ddply(mirpath, c("class", "notable_pathway", "cluster"), summarise,
                         gene_number = sum(gene_number)
)
class_summery <- ddply(mirpath, c("class", "cluster"), summarise,
                       gene_number = sum(gene_number)
)
mirpath_summery$class <- factor(mirpath_summery$class, levels = class_summery[order(class_summery$gene_number, decreasing=F), "class"])
mirpath_summery$notable_pathway <- factor(mirpath_summery$notable_pathway, levels=c(notable_pathway,"others"))

mirpath_summery_2_3 <- subset(mirpath_summery, cluster %in% c("cluster 2", "cluster 3"))

mirpath_4 <- subset(mirpath, 
                    miRNA %in% c("miR-486-5p", "miR-92a-3p", "miR-126-5p", "miR-27b-3p")
                      )
mirpath_summery_4 <- ddply(mirpath_4, c("class", "notable_pathway", "cluster"), summarise,
                         gene_number = sum(gene_number)
)

pdf("mirpathway_class.pdf", width=9, height=10)
ggplot(mirpath_summery_2_3, aes(x=class, y=gene_number, fill=notable_pathway)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=getPalette(colourCount), name="Pathway")  +
  theme(axis.text.y = element_text(size=7, colour="#000000")) +
  theme(axis.text.x = element_text(size=7, colour="#000000")) +
  theme(legend.text=element_text(size=7), legend.title = element_text(size=8)) +
  coord_flip() +
  ylab("Number of target genes") +
  xlab("Pathway") +
  facet_wrap(~cluster, ncol=3)
dev.off()  

mirpath$plasma_abundance <- as.numeric(mirpath$plasma_abundance)
mirpath$miRNA <- factor(mirpath$miRNA, levels=unique(as.vector(arrange(mirpath, desc(plasma_abundance))$miRNA)))
ordered_RNA <- gsub("hsa-", "", levels(mirpath$miRNA))
mirpath$miRNA <- factor(mirpath$miRNA, levels=ordered_RNA) 
mirpath_2_3 <- subset(mirpath, cluster %in% c("cluster 2", "cluster 3"))
miRNA_levels <- levels(mirpath_2_3$miRNA)

# write.table(mirpath_2_3, file = "mirpath_2_3.csv", quote=F, sep="\t") 
#order data with libraoffice
mirpath_2_3 <- read.csv("mirpath_2_3.csv", header=T, sep="\t")
mirpath_2_3$pathway <- factor(mirpath_2_3$pathway, levels = rev(unique(as.vector(mirpath_2_3$pathway))))
mirpath_2_3$miRNA <- factor(mirpath_2_3$miRNA, levels = miRNA_levels)

p_pathway <- ggplot(mirpath_2_3, aes(x=miRNA, y=pathway)) +
  geom_point(aes(color=log10(plasma_abundance), size=gene_number)) +
  scale_size_area(max_size = 3.5, name="Number of\ntarget genes") +
  #   scale_colour_gradientn(colours=c("red","violet","blue")) +
  #   geom_text(aes(y=as.numeric(pathway), label= round(gene_number,0)), hjust=1.3,
  #             size=2) +
  scale_color_gradient( low = "#132B43",
                        high = "red", name="Relative\nnormalized\nread count\n(log10)") +
  theme_minimal() +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=90, vjust=0.7, hjust=0.7)) +
  theme(legend.text = element_text(size=6), legend.title=element_text(size=6.5)) +
  xlab("miRNA") +
  ylab("Pathway") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
#   ggtitle("A") +
  theme(plot.title = element_text(hjust = 0)) +
  facet_wrap(~cluster, scale="free_x")

pdf("fig5_pathway.pdf")
p_pathway
dev.off()

top5_pathway <- levels(mirpath$pathway)[54:58]
Neurotrophin_signaling_pathway_mirs <- as.vector(subset(mirpath, pathway=="Neurotrophin_signaling_pathway")$miRNA)
MAPK_signaling_pathway_mirs <- as.vector(subset(mirpath, pathway=="MAPK_signaling_pathway")$miRNA)
Focal_adhesion_mirs <- as.vector(subset(mirpath, pathway=="Focal_adhesion")$miRNA)
Pathways_in_cancer_mirs <- as.vector(subset(mirpath, pathway=="Pathways_in_cancer")$miRNA)
PI3K_Akt_signaling_pathway_mirs <- as.vector(subset(mirpath, pathway=="PI3K-Akt_signaling_pathway")$miRNA)
intersect(c(Neurotrophin_signaling_pathway_mirs, MAPK_signaling_pathway_mirs), 
          c(Focal_adhesion_mirs, Pathways_in_cancer_mirs,
            PI3K_Akt_signaling_pathway_mirs))

top5_mirpath <- subset(mirpath, pathway %in% top5_pathway) 
top5_mirs <- unique(top5_mirpath$miRNA)

key_seqs$individual_id <-  as.factor(toupper(key_seqs$individual_id))
key_seqs$condition <- as.factor(sapply(key_seqs$condition, gsub, pattern="exosome145", replacement="exosome"))

key_series_plot <- subset(key_seqs,
                          sequence %in% c(
                            "GGCTGGTCCGATGGTAGTGGGTTATCAGAACT",
                            "GTTTCCGTAGTGTAGTGGTTATCACGTTCGCCT",
                            "CGCGACCTCAGATCAGACGTGGCGACCCGCTGAAT",
                            "TATTGCACTTGTCCCGGCCTGT",
                            "TCCTGTACTGAGCTGCCCCGAG",
                            #                                  "TGTAAACATCCTTGACTGGAAGCT",
                            "TTCACAGTGGCTAAGTTCTG",
                            #                                  "CATTGCACTTGTCTCGGTCTGA",
                            #                                  "AAGCTGCCAGTTGAAGAACTGT",
                            "TAGCTTATCAGACTGATGTTGA",
                            "TTTGGCAATGGTAGAACTCACA",
                            "TAGCAGCACGTAAATATTGGCG",
                            "TAGCAGCACATAATGGTTTG",
                            "CATTATTACTTTTGGTACGCG",
                            "AAGCTGCCAGTTGAAGAACTGT"
                          ))

key_series_9_seqs <- subset(key_seqs,
                            rna %in% c("miR-486-5p", "yRF-Y4-5p",
                                       "miR-92a-3p", "miR-16-5p",
                                       "miR-21-5p", "miR-30e-5p-3'R2",
                                       "miR-126-5p", "yRF-Y4-3p",
                                       "tRF-tRNA-Val(CAC/AAC)", "miR-22-3p"
                                       ))

key_series_9_seqs$rna <- droplevels(key_series_9_seqs$rna)

library(scales)
time_dynamics <- function(name) {
ggplot(subset(key_series_9_seqs, rna==name), aes(x=time_point, y=log(abundance+1,10), group=condition, color=condition)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=time_point, y=log(abundance+1,10)), size=1) + 
  #   scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF"), name= "RNA") +
  scale_x_datetime(breaks = unique(key_seqs$time_point), 
                   labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_brewer(palette="Paired", name="Sample") +
  #   scale_color_brewer(palette="Set1", name="RNA") +
  xlab("Time point") + ylab("Normalized relative read count (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~individual_id, nrow=4, scale = "free_y") +
  ggtitle(name) +
  theme(plot.title = element_text(hjust = 0)) 
}

pdf("dynamics_time_9seqs.pdf")
time_dynamics("yRF-Y4-5p")
time_dynamics("miR-21-5p")
time_dynamics("miR-126-5p")
time_dynamics("yRF-Y4-3p")
time_dynamics("tRF-tRNA-Val(CAC/AAC)")
time_dynamics("miR-486-5p")
time_dynamics("miR-92a-3p")
time_dynamics("miR-16-5p")
time_dynamics("miR-30e-5p-3'R2")
time_dynamics("miR-22-3p")
dev.off()

key_series_y4_5p <- subset(key_series_plot,
                           rna == "yRF-Y4-5p")
# wide_form <- reshape(key_series_y4_5p, idvar = c("sequence", "rna",  "sample_id"), timevar = "condition", direction = "wide")
# cor(wide_form$abundance.plasma, wide_form$abundance.exosome145)
# plot(wide_form$abundance.plasma, wide_form$abundance.exosome145)

p_y4_5p <- ggplot(key_series_y4_5p, aes(x=time_point, y=log(abundance+1,10), group=condition, color=condition)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=time_point, y=log(abundance+1,10)), size=1) + 
  #   scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF"), name= "RNA") +
  scale_x_datetime(breaks = unique(key_seqs$time_point), 
                   labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_brewer(palette="Paired", name="Sample") +
  #   scale_color_brewer(palette="Set1", name="RNA") +
  xlab("Time point") + ylab("Normalized relative read count (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~individual_id, nrow=4, scale = "free_y") +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0)) 



key_series_mir_486_5p <- subset(key_series_plot,
                           rna == "miR-486-5p")

p_486_5p <- ggplot(key_series_mir_486_5p, aes(x=time_point, y=log(abundance+1,10), group=condition, color=condition)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=time_point, y=log(abundance+1,10)), size=1) + 
  scale_x_datetime(breaks = unique(key_seqs$time_point), 
                   labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_brewer(palette="Paired", name="Sample", guide=F) +
  xlab("Time point") + ylab("Normalized relative read count (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~individual_id, nrow=4, scale = "free_y") +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0)) 

#load 150319_qPCR-Anna-plots.Rdata for the experimental results
# write.table(all, file = "diurnal-dynamics-qPCR-Anna-20150323.csv", quote=F, sep="\t", row.names=F)

diurnal_qPCR <- read.csv("diurnal-dynamics-qPCR-Anna-20150323.csv",  header=T, sep="\t")
diurnal_qPCR$timeF <- as.POSIXct(diurnal_qPCR$timeF)
diurnal_qPCR$sampleType <- factor(diurnal_qPCR$sampleType, levels=c("plasma", "exosome fraction", "protease treated exosome fraction"))
  
yRF5p <- subset(diurnal_qPCR, RNA=="yRF5p" )
yRF5p_qPCR_fig <- ggplot(yRF5p, aes(x=timeF, y=log(normalized.abundance,10), group=sampleType, color=sampleType)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=timeF, y=log(normalized.abundance,10)), size=1) + 
  scale_x_datetime(breaks = unique(yRF5p$timeF), 
                   labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"), name="Sample type") +
  xlab("Time point") + ylab("Normalized relative RNA abundance (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~Individual, nrow=4, scale = "free_y") +
  ggtitle("D") +
  theme(plot.title = element_text(hjust = 0))

pdf("yRF5p-qPCR-Anna-20150323.pdf")
yRF5p_qPCR_fig
dev.off()

miR486_qPCR <- subset(diurnal_qPCR, RNA=="miR486" )
miR486_qPCR_fig <- ggplot(miR486_qPCR, aes(x=timeF, y=log(normalized.abundance,10), group=sampleType, color=sampleType)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=timeF, y=log(normalized.abundance,10)), size=1) + 
  scale_x_datetime(breaks = unique(yRF5p$timeF), 
                   labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"), 
                     name="Sample type") +
  xlab("Time point") + ylab("Normalized relative RNA abundance (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~Individual, nrow=4, scale = "free_y") +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0))

pdf("miR486-qPCR-Anna-20150323.pdf")
miR486_qPCR_fig 
dev.off()

library(reshape)
absolute_number_qPCR <- read.csv("absoluteNumbersExosomePlasma_anna.csv",  header=T, sep="\t")
absolute_number_qPCR_long <- melt(absolute_number_qPCR[, c(1, 2, 3, 7)])
absolute_number_qPCR_long$individual <- absolute_number_qPCR_long$sample
absolute_number_qPCR_long$individual <- substr(absolute_number_qPCR_long$individual, 1, 3)
absolute_number_qPCR_long$time_point <- absolute_number_qPCR_long$sample
absolute_number_qPCR_long$time_point <- substr(absolute_number_qPCR_long$time_point, 5, 5)
absolute_number_qPCR_long$time_point <- gsub( "4", "10:30", absolute_number_qPCR_long$time_point)
absolute_number_qPCR_long$time_point <- gsub( "5", "12:00", absolute_number_qPCR_long$time_point)
# absolute_number_qPCR_long$time_point <- as.POSIXct(absolute_number_qPCR_long$time_point)
absolute_number_qPCR_long$compartment <- gsub( "crude exosomes", "exosome fraction", absolute_number_qPCR_long$compartment)
absolute_number_qPCR_long$compartment <- gsub( "protease treated plasma exosomes", "protease treated exosome fraction", absolute_number_qPCR_long$compartment)

absolute_number_qPCR_long$compartment <- factor(absolute_number_qPCR_long$compartment, levels = c("plasma", "exosome fraction", 
                                                                                                  "protease treated exosome fraction"))

yRF5p_abs <- subset(absolute_number_qPCR_long, target=="yRF-Y4-5p")
yRF5p_abs$time_point <- as.factor(sapply(yRF5p_abs$time_point, gsub, pattern="10:30", replacement="10:30 AM"))
yRF5p_abs$time_point <- as.factor(sapply(yRF5p_abs$time_point, gsub, pattern="12:00", replacement="12:00 PM"))

yRF5p_abs_fig <- ggplot(yRF5p_abs, aes(x=time_point, y=log(value,10), fill=compartment)) +
  theme_minimal() +
  geom_boxplot() +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_fill_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"), name="Sample type") +
  xlab("Time point") + ylab("Number of molecules per milliliter plasma (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0))

miR486_abs <- subset(absolute_number_qPCR_long, target=="miR-486-5p")
miR486_abs$time_point <- as.factor(sapply(miR486_abs$time_point, gsub, pattern="10:30", replacement="10:30 AM"))
miR486_abs$time_point <- as.factor(sapply(miR486_abs$time_point, gsub, pattern="12:00", replacement="12:00 PM"))

miR486_abs_fig <- ggplot(miR486_abs, aes(x=time_point, y=log(value,10), fill=compartment)) +
  theme_minimal() +
#   geom_bar(position="dodge") +
  geom_boxplot() +                         
#   scale_x_datetime(breaks = unique(absolute_number_qPCR_long$time_point), 
#                    labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_fill_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"), name="Sample type") +
  xlab("Time point") + ylab("Number of molecules per milliliter plasma (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
#   facet_wrap(~individual, nrow=4, scale = "free_y") +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0))


library(grid)

# Define layout for the plots (2 rows, 2 columns)
layt <- grid.layout(nrow = 2, ncol = 3, heights = c(4/8, 4/8), widths = c(3/9, 
                                                                          3/9, 3/9), default.units = c("null", "null"))
# View the layout of plots
# grid.show.layout(layt)

tmp <- ggplotGrob(p_y4_5p)
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
legend$vp$x <- unit(0.5, 'npc')
legend$vp$y <- unit(0.2, 'npc')

ylab <- tmp$grobs[[19]]
ylab$vjust <- 8

pdf("fig3_Y4_mir486_dynamics.pdf", width=10, height=7)
# Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout = layt))
print(p_486_5p + guides(fill=FALSE), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_y4_5p + theme(legend.position="none"), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(miR486_qPCR_fig + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(yRF5p_qPCR_fig + theme(legend.position="none"), vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(miR486_abs_fig + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(yRF5p_abs_fig + theme(legend.position="none"), vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
dev.off()

#!!!---20150702

#data from qPCR to describe dynamics
#data from RNA-Seq to describe dynamics
#data from qPCR to describe the abundance distribution in the three compartments
#combine the results into one figure


cluster1 <- c(
  "miR-22-3p",
  "miR-25-3p",
  "miR-423-5p",
  "miR-142-5p-5'L2-3'L3",
  "miR-192-5p",
  "miR-30e-5p-3'R2" #qPCR
)
cluster3 <- c(
  'yRF-Y4-5p', #qPCR
  'miR-191-5p',
  "miR-27b-3p-3'L1",
  'miR-126-5p',
  'miR-21-5p', #qPCR
  "miR-143-3p-3'L1",
  "miR-10a-5p-3'L1",
  'tRF-tRNA-Val(CAC/AAC)', #qPCR
  'tRF-tRNA-Gly(GCC/CCC)'
)
cluster2 <- c(
  'miR-486-5p', #qPCR
  "miR-486-5p-3'L1",
  "miR-486-5p-3'~A",
  "miR-486-5p-3'~U",
  "miR-486-5p-3'L2",
  "miR-486-5p-3'L1~A",
  "miR-486-5p-3'~AA",
  "miR-486-5p-3'R1",
  'miR-92a-3p', #qPCR 
  "miR-92a-3p-3'L1",
  "miR-92a-3p-3'L1~A",
  'miR-16-5p', #qPCR
  "miR-451a-3'L1",
  'miR-451a',
  "miR-451a-3'R1",
  "miR-451a-3'L2",
  "miR-22-3p-3'L1",
  "miR-25-3p-3'L2",
  "miR-423-5p-3'L2",
  "miR-423-5p-3'L1",
  "miR-103a/107-3p-3'L4",
  "miR-15a-5p-3'L2",
  "miR-182-5p-3'L2"
)
clusterNA <- c(
  "yRF-Y4-5p-3'L1",
  "yRF-Y4-5p-3'R1",
  'yRF-Y4-3p', #qPCR
  "yRF-Y4-3p-3'L4",
  "miR-10b-5p-3'L1",
  'rRF-RNA28S5'
)


diurnal_qPCR <- read.csv("diurnal-dynamics-qPCR-Anna-20150323.csv",  header=T, sep="\t")
diurnal_qPCR$timeF <- as.POSIXct(diurnal_qPCR$timeF)
diurnal_qPCR$sampleType <- factor(diurnal_qPCR$sampleType, levels=c("plasma", "exosome fraction", "protease treated exosome fraction"))

diurnal_qPCR_plasma <- subset(diurnal_qPCR, sampleType=="plasma")
diurnal_qPCR_plasma$sampleType <- droplevels(diurnal_qPCR_plasma$sampleType)

diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="miR126", replacement="miR-126-5p"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="miR16", replacement="miR-16-5p"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="miR30e", replacement="miR-30e-5p-3'R2"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="miR486", replacement="miR-486-5p"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="miR92", replacement="miR-92a-3p"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="tRFval", replacement="tRF-tRNA-Val(CAC/AAC)"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="miR126", replacement="miR-126-5p"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="yRF3p", replacement="yRF-Y4-3p"))
diurnal_qPCR_plasma$RNA <- as.factor(sapply(diurnal_qPCR_plasma$RNA, gsub, pattern="yRF5p", replacement="yRF-Y4-5p"))

cluster3_plasma <- subset(diurnal_qPCR_plasma, RNA %in% cluster3)

cluster3_summary <- ddply(cluster3_plasma, c("RNA", "Time.Point", "timeF", "sampleType"), summarise, 
                       N  = sum(!is.na(normalized.abundance)),
                       mean=mean(normalized.abundance, na.rm=TRUE), 
                       sd=sd(normalized.abundance, na.rm=TRUE), 
                       se=sd/sqrt(N),
                       min=log10(mean-se),
                       max=log10(mean+se),
                       log10_mean=log(mean, 10)) 

cluster3_summary$RNA <- factor(cluster3_summary$RNA, levels=c("yRF-Y4-5p", "miR-126-5p", "tRF-tRNA-Val(CAC/AAC)"))


cluster3_dynamics <- ggplot(cluster3_summary, aes(x=timeF, y=log10_mean, group=sampleType, color=sampleType)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=timeF, y=log10_mean), size=1.5) + theme(legend.position="none") +
  geom_errorbar(aes(ymax = max, ymin=min), position="dodge", width=0.5) +
  scale_x_datetime(breaks = unique(cluster3_summary$timeF), 
                   labels = date_format("%I:%M %p")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=75, hjust=0.9)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"), name="Sample type", guide=FALSE) +
  xlab("Time point") + ylab("Normalized relative RNA abundance \n (log10(mean)  \u00B1 s.d.)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~RNA, nrow=4, scale = "free_y") +
  theme(strip.text = element_text(size=8)) +
  ggtitle("F") +
  theme(plot.title = element_text(hjust = 0))

cluster2_plasma <- subset(diurnal_qPCR_plasma, RNA %in% cluster2)

cluster2_summary <- ddply(cluster2_plasma, c("RNA", "Time.Point", "timeF", "sampleType"), summarise, 
                          N  = sum(!is.na(normalized.abundance)),
                          mean=mean(normalized.abundance, na.rm=TRUE), 
                          sd=sd(normalized.abundance, na.rm=TRUE), 
                          se=sd/sqrt(N),
                          min=log10(mean-se),
                          max=log10(mean+se),
                          log10_mean=log(mean, 10)) 

cluster2_summary$RNA <- factor(cluster2_summary$RNA, levels=c("miR-486-5p", "miR-92a-3p", "miR-16-5p"))


cluster2_dynamics <- ggplot(cluster2_summary, aes(x=timeF, y=log10_mean, group=sampleType, color=sampleType)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=timeF, y=log10_mean), size=1.5) + theme(legend.position="none") +
  geom_errorbar(aes(ymax = max, ymin=min), position="dodge", width=0.5) +
  scale_x_datetime(breaks = unique(cluster2_summary$timeF), 
                   labels = date_format("%I:%M %p")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=75, hjust=0.9)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"), name="Sample type", guide=FALSE) +
  xlab("Time point") + ylab("Normalized relative RNA abundance \n (log10(mean)  \u00B1 s.d.)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~RNA, nrow=4, scale = "free_y") +
  theme(strip.text = element_text(size=8)) +
  ggtitle("E") +
  theme(plot.title = element_text(hjust = 0))

key_series_9_seqs_plasma <- subset(key_seqs,
                            rna %in% c("miR-486-5p", "yRF-Y4-5p",
                                       "miR-92a-3p", "miR-16-5p",
                                       "miR-21-5p", "miR-30e-5p-3'R2",
                                       "miR-126-5p", "yRF-Y4-3p",
                                       "tRF-tRNA-Val(CAC/AAC)", "miR-22-3p"
                            ) &
                              condition == "plasma")

key_series_9_seqs_plasma$rna <- droplevels(key_series_9_seqs$rna)
key_series_9_seqs_plasma$condition <- droplevels(key_series_9_seqs_plasma$condition)

key_series_9_seqs_plasma_summary <- ddply(key_series_9_seqs_plasma, c("rna", "time_point",  "condition"), summarise, 
      N  = sum(!is.na(abundance)),
      mean=mean(abundance, na.rm=TRUE), 
      sd=sd(abundance, na.rm=TRUE), 
      se=sd/sqrt(N),
      min=log10(mean-se),
      max=log10(mean+se),
      log10_mean=log(mean, 10)) 

cluster2_RNAseq_summary <- subset(key_series_9_seqs_plasma_summary, rna %in% cluster2)
cluster2_RNAseq_summary$rna <- factor(cluster2_RNAseq_summary$rna, levels=c("miR-486-5p", "miR-92a-3p", "miR-16-5p"))
cluster2_dynamics_RNAseq <- ggplot(cluster2_RNAseq_summary, aes(x=time_point, y=log10_mean, group=condition, color=condition)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=time_point, y=log10_mean), size=1.5) + theme(legend.position="none") +
  geom_errorbar(aes(ymax = max, ymin=min), position="dodge", width=0.5) +
  scale_x_datetime(breaks = unique(cluster2_RNAseq_summary$time_point), 
                   labels = date_format("%I:%M %p")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=75,  hjust=0.9)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"),  guide=FALSE) +
  xlab("Time point") + ylab("Normalized relative RNA abundance \n (log10(mean)  \u00B1 s.d.)")+
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~rna, nrow=4, scale = "free_y") +
  theme(strip.text = element_text(size=8)) +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0))

cluster3_RNAseq_summary <- subset(key_series_9_seqs_plasma_summary, rna %in% c("yRF-Y4-5p", "miR-126-5p", "tRF-tRNA-Val(CAC/AAC)"))
cluster3_RNAseq_summary$rna <- factor(cluster3_RNAseq_summary$rna, levels=c("yRF-Y4-5p", "miR-126-5p", "tRF-tRNA-Val(CAC/AAC)"))
cluster3_dynamics_RNAseq <- ggplot(cluster3_RNAseq_summary, aes(x=time_point, y=log10_mean, group=condition, color=condition)) +
  theme_minimal() +
  geom_line(size=0.5) +  
  geom_point(aes(x=time_point, y=log10_mean), size=1.5) + theme(legend.position="none") +
  geom_errorbar(aes(ymax = max, ymin=min), position="dodge", width=0.5) +
  scale_x_datetime(breaks = unique(cluster3_RNAseq_summary$time_point), 
                   labels = date_format("%I:%M %p")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=75, hjust=0.9)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_manual(values=c("#1F78B4", "#A6CEE3", "#B2DF8A"),  guide=FALSE) +
  xlab("Time point") + ylab("Normalized relative RNA abundance \n (log10(mean)  \u00B1 s.d.)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~rna, nrow=4, scale = "free_y") +
  theme(strip.text = element_text(size=8)) +
  ggtitle("D") +
  theme(plot.title = element_text(hjust = 0))

library(grid)

# Define layout for the plots (2 rows, 2 columns)
layt <- grid.layout(nrow = 2, ncol = 3, heights = c(4/8, 4/8), widths = c(3/9, 
                                                                          3/9, 3/9), default.units = c("null", "null"))
# View the layout of plots
# grid.show.layout(layt)

tmp <- ggplotGrob(p_y4_5p)
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
legend$vp$x <- unit(0.5, 'npc')
legend$vp$y <- unit(0.2, 'npc')

ylab <- tmp$grobs[[19]]
ylab$vjust <- 8

pdf("fig3_cluster2-3_dynamics.pdf", width=10, height=7)
# Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout = layt))
print(miR486_abs_fig + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(yRF5p_abs_fig + theme(legend.position="none"), vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(cluster2_dynamics_RNAseq + guides(fill=FALSE), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(cluster3_dynamics_RNAseq + theme(legend.position="none"), vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(cluster2_dynamics + theme(legend.position="none"), vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(cluster3_dynamics + theme(legend.position="none"), vp = viewport(layout.pos.row = 2, layout.pos.col = 3))
dev.off()
