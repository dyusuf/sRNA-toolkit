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
  {mirpath[mirpath["miRNA"] == i, ]$cluster <- "cluster1"} 
else if (i %in% cluster2) 
  {mirpath[mirpath["miRNA"] == i, ]$cluster <- "cluster2"} 
else if (i %in% cluster3) 
{mirpath[mirpath["miRNA"] == i, ]$cluster <- "cluster3"} 
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

mirpath_4 <- subset(mirpath, 
                    miRNA %in% c("miR-486-5p", "miR-92a-3p", "miR-126-5p", "miR-27b-3p")
                      )
mirpath_summery_4 <- ddply(mirpath_4, c("class", "notable_pathway", "cluster"), summarise,
                         gene_number = sum(gene_number)
)

pdf("mirpathway_class.pdf", width=9, height=10)
ggplot(mirpath_summery_4, aes(x=class, y=gene_number, fill=notable_pathway)) +
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
p_pathway <- ggplot(mirpath, aes(x=miRNA, y=pathway)) +
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
  xlab("MiRNA") +
  ylab("Pathway") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0)) 

pdf("test.pdf")
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
                            "CATTATTACTTTTGGTACGCG"
                          ))

key_series_9_seqs <- subset(key_seqs,
                            rna %in% c("miR-486-5p", "yRF-Y4-5p",
                                       "miR-92a-3p", "miR-16-5p",
                                       "miR-21-5p", "miR-30e-5p-3'R2",
                                       "miR-126-5p", "yRF-Y4-3p",
                                       "tRF-tRNA-Val(CAC/AAC)"
                                       ))

key_series_9_seqs$rna <- droplevels(key_series_9_seqs$rna)

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
  xlab("Time point") + ylab("Relative normalized read count (log10)") +
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
  xlab("Time point") + ylab("Relative normalized read count (log10)") +
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
  #   scale_color_manual(values=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF"), name= "RNA") +
  scale_x_datetime(breaks = unique(key_seqs$time_point), 
                   labels = date_format("%H:%M")) +
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=45, vjust=0.4)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  scale_color_brewer(palette="Paired", name="Sample", guide=F) +
  #   scale_color_brewer(palette="Set1", name="RNA") +
  xlab("Time point") + ylab("Relative normalized read count (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~individual_id, nrow=4, scale = "free_y") +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0)) 

# Define layout for the plots (2 rows, 2 columns)
layt <- grid.layout(nrow = 2, ncol = 2, heights = c(4/8, 4/8), widths = c(4.5/8, 
                                                                          3.5/8), default.units = c("null", "null"))
# View the layout of plots
# grid.show.layout(layt)

tmp <- ggplotGrob(p_y4_5p)
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
legend$vp$x <- unit(0.5, 'npc')
legend$vp$y <- unit(0.2, 'npc')

ylab <- tmp$grobs[[19]]
ylab$vjust <- 8

pdf("pathyway_dynamics.pdf", width=10, height=7)
# Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout = layt))
print(p_pathway, vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1))
print(p_y4_5p+theme(legend.position="none")+ylab(NULL), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(p_486_5p+ylab(NULL), vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(grid.draw(legend))
print(grid.draw(ylab))
dev.off()

ggplot(key_series_plot, aes(x=time_point, y=log(abundance+1,10), group=condition, color=condition)) +
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
  xlab("Time point") + ylab("Relative normalized read count (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~individual_id, nrow=4, scale = "free_y") +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0)) +
  facet_wrap(~rna)