library(plyr)
library(ggplot2)

setwd("/home/dilmurat/Work/combined-experiments")
len_combined <- read.csv("2.sample_distri-read_len.csv", header=T, sep="\t")
# len_human <- subset(len_combined, rank1 =="Homo_sapiens" & condition=="p")
# len_human$condition <- droplevels(len_human$condition)
len_human <- subset(len_combined, rank1 =="Homo_sapiens")
len_human$rank1 <- droplevels(len_human$rank1)
len_human <- subset(len_human, condition == "p")
len_human$condition <- droplevels(len_human$condition)

len_human$RNA_family <- gsub("_", " ", len_human$RNA_family)
len_human[len_human$RNA_family!="miRNA",]$RNA_family <- paste(len_human[len_human$RNA_family!="miRNA",]$RNA_family, "fragment")


RNA_human_sum <- ddply(len_human, c("RNA_family", "sample_id", "condition"), summarise, total=sum(len_abundance))
RNA_human_sum$condition <- droplevels(RNA_human_sum$condition)

RNA_human_sd <- ddply(RNA_human_sum, c("RNA_family", "condition"), summarise, 
                       N    = sum(!is.na(total)),
                       mean=mean(total, na.rm=TRUE), 
                       sd=sd(total, na.rm=TRUE), 
                       se=sd/sqrt(N),
                       min=log10(mean-se),
                       max=log10(mean+se),
                       log10_mean=log(mean, 10),
                       cv=sd/mean)

RNA_human_sd$RNA_family <- factor(RNA_human_sd$RNA_family, levels=rev(unique(as.vector(arrange(RNA_human_sd, desc(mean))$RNA_family))))
RNA_human_sum$RNA_family <- factor(RNA_human_sum$RNA_family, levels=levels(RNA_human_sd$RNA_family))

RNA_composition <- ggplot(RNA_human_sum, aes(x=RNA_family, y=log10(total))) +
    geom_boxplot(outlier.size=1) +
#   geom_bar(stat='identity', position='dodge') +
#   geom_segment( aes(xend=RNA_family), yend=0, colour="grey50") +
#   geom_point(size=2, shape=1)  +
#   geom_errorbar(aes(ymin=min, ymax=max),
#                 width=.2,                    # Width of the error bars
#                 position=position_dodge(.9)) +
  theme_minimal() +
  theme(axis.text.x = element_text(size=6, colour="#000000")) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  ylab("Normalized relative read count (log10)") +
  xlab("ex-sRNA class") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  scale_y_continuous(breaks=seq(0, 15, 1)) +
  coord_flip() +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0)) 

pdf("RNA_composition-20150310.pdf")
RNA_composition
dev.off()

RNA_summery <- ddply(len_human, c("sequence_length", "condition", "RNA_family"), 
                     summarise, 
                     abundance_average=mean(len_abundance), 
                     top_percentage=mean(top_percentage),
                     sd=sd(len_abundance, na.rm=TRUE),
                     N  = sum(!is.na(len_abundance)),
                     se=sd/sqrt(N),
                     cv=sd/abundance_average)
RNA_fraction <- ddply(RNA_summery,  c("sequence_length"), transform,
                      percentage=abundance_average/sum(abundance_average),
                      se_percentage = se/sum(abundance_average)             
                      )


RNA_fraction_3RNAs <- subset(RNA_fraction,
                         RNA_family %in% c("Y RNA fragment", "miRNA", "tRNA fragment"))

RNA_fraction_3RNAs$RNA_family <- factor(RNA_fraction_3RNAs$RNA_family, levels= c("miRNA", "Y RNA fragment", "tRNA fragment")) 



selective_RNA <-ggplot(RNA_fraction_3RNAs, aes(x=sequence_length, y=log10(abundance_average), group=condition)) +
  #   geom_point(size = 4, position=position_dodge(width=0.1)) + 
  theme_minimal() +
  geom_line() +  
  geom_point(aes(x=sequence_length, y=log10(abundance_average),  size=top_percentage), alpha=0.5) + 
  scale_shape_manual(values= c(16, 22, 23, 24,   32),  name = "Most abundant\nsequence\nat given length") +
  scale_size_area(name="Average fraction of most\nabundant sequence\nat given length") +
  #   geom_text(x=22, y=-6, parse=TRUE, aes(label=label), data=data.frame(RNA_family = c("miRNA"), 
  #                                    label = c(paste("italic(CV) == ",round( 0.8070114, 2)))), size=3) +
  scale_x_continuous(breaks=seq(15, 100, 2)) +
  theme(axis.text.x = element_text(size=6, colour="#000000")) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  xlab("Sequence length (nucleotides)") + ylab("Normalized relative read count (log10)") +
  theme(axis.title = element_text(size=8, colour="#000000")) +
#   annotate("text", x=33, y=5,
#            label="yRF-Y4-5p", hjust=0.1, size=1.8)  +
#   annotate("text", x=29, y=4.25,
#            label="yRF-Y4-5p-3'L1", hjust=0.7, size=1.8)  +
#   annotate("text", x=25, y=2.7,
#            label="yRF-Y4-3p", hjust=0.6, size=1.8) +
  ggtitle("C") +
  theme(plot.title = element_text(hjust = 0)) +
  facet_wrap(~RNA_family, nrow=3)

pdf("selective_RNA_size_distri-20150311.pdf")
selective_RNA
dev.off()


RNA_fraction_rRNA <- subset(RNA_fraction,
                            RNA_family == "rRNA fragment" & sequence_length <= 47)
p_rRNA <- ggplot(RNA_fraction_rRNA, aes(x=sequence_length, y=log10(abundance_average), group=condition)) +
  #   geom_point(size = 4, position=position_dodge(width=0.1)) + 
  theme_minimal() +
  geom_line() +
  geom_point(aes(x=sequence_length, y=log10(abundance_average),  size=top_percentage), alpha=0.5) + 
  #   scale_shape_manual(values= c(16, 22, 23, 24,   32),  name = "Most abundant\nsequence\nat given length") +
  scale_size_area(name="Average fraction of most\nabundant sequence\nat given length") +
  #   geom_text(x=22, y=-6, parse=TRUE, aes(label=label), data=data.frame(RNA_family = c("miRNA"), 
  #                                    label = c(paste("italic(CV) == ",round( 0.8070114, 2)))), size=1) +
  scale_x_continuous(breaks=seq(15, 100, 2)) +
  theme(axis.text.x = element_text(size=6, colour="#000000")) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  xlab("Sequence length (nucleotides)") + ylab("Normalized relative read count (log10)") +
  theme(axis.title.y = element_text(size=8, colour="#000000")) +
  theme(axis.title.x = element_text(size=8, colour="#000000")) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0)) +
  facet_wrap(~RNA_family)

# Define layout for the plots (2 rows, 2 columns)
layt <- grid.layout(nrow = 2, ncol = 2, heights = c(3/6, 3/6), widths = c(3/6, 
                                                                          3/6), default.units = c("null", "null"))
# View the layout of plots
# grid.show.layout(layt)

pdf("fig1_exRNA_charateristics-20150311.pdf", width=10, height=7)
grid.newpage()
pushViewport(viewport(layout = layt))
print(RNA_composition, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_rRNA, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
print(selective_RNA, vp = viewport(layout.pos.row = 1:2, layout.pos.col = 2))
dev.off()


# len_human["variability"] <- "NA"
# for (i in 1:length(len_human$sample_id)) {if (len_human[i,]$condition != "p") 
#                                           len_human[i,]$variability <- toupper(substr(len_human[i,]$sample_id, 1, 3))
#                                           else len_human[i,]$variability <- "C10"}
# 
# RNA_human_sum <- ddply(len_human, c("RNA_family", "sample_id", "condition", "variability"), summarise, total=sum(len_abundance))
# RNA_human_sum$condition <- droplevels(RNA_human_sum$condition)
# 
# RNA_human_sd <- ddply(RNA_human_sum, c("RNA_family", "condition", "variability"), summarise, 
#                        N    = sum(!is.na(total)),
#                        mean=mean(total, na.rm=TRUE), 
#                        sd=sd(total, na.rm=TRUE), 
#                        se=sd/sqrt(N),
#                        min=log10(mean-se),
#                        max=log10(mean+se),
#                        log10_mean=log(mean, 10) ) 
# RNA_human_sd$RNA_family <- factor(RNA_human_sd$RNA_family, levels=rev(unique(as.vector(arrange(RNA_human_sd, desc(mean))$RNA_family))))
# RNA_human_sum$RNA_family <- factor(RNA_human_sum$RNA_family, levels=levels(RNA_human_sd$RNA_family))
# RNA_human_sum$variability <- factor(RNA_human_sum$variability, levels=c("C10", "P01", "P02", "P03", "P04"))
# RNA_composition <- ggplot(RNA_human_sum, aes(x=RNA_family, y=log10(total), fill=variability)) +
#     geom_boxplot(outlier.size=1) +
# #   geom_bar(stat='identity', position='dodge') +
# #   geom_segment( aes(xend=RNA_family), yend=0, colour="grey50") +
# #   geom_point(size=2, shape=1)  +
# #   geom_errorbar(aes(ymin=min, ymax=max),
# #                 width=.2,                    # Width of the error bars
# #                 position=position_dodge(.9)) +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size=6, colour="#000000")) +
#   theme(axis.text.y = element_text(size=6, colour="#000000")) +
#   scale_fill_brewer(palette="Set1", name="Dataset") +
#   ylab("Normalized relative read count (log10)") +
#   xlab("RNA class") +
#   theme(axis.title = element_text(size=8, colour="#000000")) +
#   scale_y_continuous(breaks=seq(0, 15, 1)) +
#   coord_flip() +
#   ggtitle("A") +
#   theme(plot.title = element_text(hjust = 0)) 
# 
# pdf("RNA_composition-20150309.pdf")
# RNA_composition
# dev.off()
# 
# RNA_summery <- ddply(len_human, c("sequence_length", "condition", "RNA_family", "variability"), 
#                      summarise, 
#                      abundance_average=mean(len_abundance), 
#                      top_percentage=mean(top_percentage),
#                      sd=sd(len_abundance, na.rm=TRUE),
#                      N  = sum(!is.na(len_abundance)),
#                      se=sd/sqrt(N),
#                      cv=sd/abundance_average)
# RNA_fraction <- ddply(RNA_summery,  c("variability"), transform,
#                       percentage=abundance_average/sum(abundance_average),
#                       se_percentage = se/sum(abundance_average)             
#                       )
# 
# 
# RNA_fraction_y <- subset(RNA_fraction,
#                          RNA_family == "Y RNA")
# p_yRNA <-ggplot(RNA_fraction_y, aes(x=sequence_length, y=log10(abundance_average), group=variability, color=variability)) +
#   #   geom_point(size = 4, position=position_dodge(width=0.1)) + 
#   theme_minimal() +
#   geom_line() +  
#   geom_point(aes(x=sequence_length, y=log10(abundance_average),  size=top_percentage), alpha=0.5) + 
#   scale_shape_manual(values= c(16, 22, 23, 24,   32),  name = "Most abundant\nsequence\nat given length") +
#   scale_size_area(name="Fraction of most\nabundant sequence\nat given length") +
#   #   geom_text(x=22, y=-6, parse=TRUE, aes(label=label), data=data.frame(RNA_family = c("miRNA"), 
#   #                                    label = c(paste("italic(CV) == ",round( 0.8070114, 2)))), size=3) +
#   scale_x_continuous(breaks=seq(15, 100, 2)) +
#   theme(axis.text.x = element_text(size=6, colour="#000000")) +
#   theme(axis.text.y = element_text(size=6, colour="#000000")) +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
#   xlab("Sequence length (nucleotides)") + ylab("Normalized relative read count (log10)") +
#   theme(axis.title = element_text(size=8, colour="#000000")) +
# #   annotate("text", x=33, y=5,
# #            label="yRF-Y4-5p", hjust=0.1, size=1.8)  +
# #   annotate("text", x=29, y=4.25,
# #            label="yRF-Y4-5p-3'L1", hjust=0.7, size=1.8)  +
# #   annotate("text", x=25, y=2.7,
# #            label="yRF-Y4-3p", hjust=0.6, size=1.8) +
#   ggtitle("B") +
#   theme(plot.title = element_text(hjust = 0)) +
#   scale_color_brewer(palette="Set1", name="Dataset")
# 
# pdf("yRF_size_distri.pdf")
# p_yRNA
# dev.off()
# 
# RNA_fraction_miRNA <- subset(RNA_fraction,
#                          RNA_family == "miRNA")
# 
# p_miRNA <- ggplot(RNA_fraction_miRNA, aes(x=sequence_length, y=log10(abundance_average), group=variability, color=variability)) +
#   #   geom_point(size = 4, position=position_dodge(width=0.1)) + 
#   theme_minimal() +
#   geom_line() +  
#   geom_point(aes(x=sequence_length, y=log10(abundance_average),  size=top_percentage), alpha=0.5) + 
#   #   scale_shape_manual(values= c(16, 22, 23, 24,   32),  name = "Most abundant\nsequence\nat given length") +
#   scale_size_area(name="Fraction of most\nabundant sequence\nat given length") +
#   #   geom_text(x=22, y=-6, parse=TRUE, aes(label=label), data=data.frame(RNA_family = c("miRNA"), 
#   #                                    label = c(paste("italic(CV) == ",round( 0.8070114, 2)))), size=3) +
#   scale_x_continuous(breaks=seq(15, 100, 2)) +
#   theme(axis.text.x = element_text(size=6, colour="#000000")) +
#   theme(axis.text.y = element_text(size=6, colour="#000000")) +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
#   xlab("Sequence length (nucleotides)") + ylab("Normalized relative read count (log10)") +
#   theme(axis.title = element_text(size=8, colour="#000000")) +
# #   annotate("text", x=21, y=6,
# #            label="miR-486-5p", hjust=-0.5, size=1.8)  +
# #   annotate("text", x=21, y=5.3,
# #            label="miR-486-5p-3'L1", hjust=1.13, size=1.8)  +
# #   annotate("text", x=23, y=5.2,
# #            label="miR-486-5p-3'~A", hjust=-0.1, size=1.8) +
#   ggtitle("C") +
#   theme(plot.title = element_text(hjust = 0)) +
#   scale_color_brewer(palette="Set1", name="Dataset")
# 
# pdf("miRNA_size_distri.pdf")
# p_miRNA
# dev.off()
# 
# RNA_fraction_tRNA <- subset(RNA_fraction,
#                                RNA_family == "tRNA")
# p_tRNA <- ggplot(RNA_fraction_tRNA, aes(x=sequence_length, y=log10(abundance_average), group=variability, color=variability)) +
#   #   geom_point(size = 4, position=position_dodge(width=0.1)) + 
#   theme_minimal() +
#   geom_line() +  
#   geom_point(aes(x=sequence_length, y=log10(abundance_average),  size=top_percentage), alpha=0.5) + 
#   #   scale_shape_manual(values= c(16, 22, 23, 24,   32),  name = "Most abundant\nsequence\nat given length") +
#   scale_size_area(name="Fraction of most\nabundant sequence\nat given length") +
#   #   geom_text(x=22, y=-6, parse=TRUE, aes(label=label), data=data.frame(RNA_family = c("miRNA"), 
#   #                                    label = c(paste("italic(CV) == ",round( 0.8070114, 2)))), size=3) +
#   scale_x_continuous(breaks=seq(15, 100, 2)) +
#   theme(axis.text.x = element_text(size=6, colour="#000000")) +
#   theme(axis.text.y = element_text(size=6, colour="#000000")) +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
#   xlab("Sequence length (nucleotides)") + ylab("Normalized relative read count (log10)") +
#   theme(axis.title = element_text(size=8, colour="#000000")) +
# #   annotate("text", x=33, y=2.45,
# #            label="tRF-tRNA-Val(CAC/AAC)", hjust=-0.1, size=1.8)  +
# # #   annotate("text", x=32, y=2.3,
# #            label="tRF-tRNA-Gly(GCC/CCC)", hjust=1.05, size=2.5)  +
#   ggtitle("E") +
#   theme(plot.title = element_text(hjust = 0))  +
#   scale_color_brewer(palette="Set1", name="Dataset")
# 
# pdf("tRF_size-distri.pdf")
# p_tRNA
# dev.off()
# 
# RNA_fraction_rRNA <- subset(RNA_fraction,
#                             RNA_family == "rRNA" & sequence_length <= 47)
# p_rRNA <- ggplot(RNA_fraction_rRNA, aes(x=sequence_length, y=log10(abundance_average), group=variability, color=variability)) +
#   #   geom_point(size = 4, position=position_dodge(width=0.1)) + 
#   theme_minimal() +
#   geom_line() +
#   geom_point(aes(x=sequence_length, y=log10(abundance_average),  size=top_percentage), alpha=0.5) + 
#   #   scale_shape_manual(values= c(16, 22, 23, 24,   32),  name = "Most abundant\nsequence\nat given length") +
#   scale_size_area(name="Fraction of most\nabundant sequence\nat given length") +
#   #   geom_text(x=22, y=-6, parse=TRUE, aes(label=label), data=data.frame(RNA_family = c("miRNA"), 
#   #                                    label = c(paste("italic(CV) == ",round( 0.8070114, 2)))), size=1) +
#   scale_x_continuous(breaks=seq(15, 100, 2)) +
#   theme(axis.text.x = element_text(size=6, colour="#000000")) +
#   theme(axis.text.y = element_text(size=6, colour="#000000")) +
#   theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
#   xlab("Sequence length (nucleotides)") + ylab("Normalized relative read count (log10)") +
#   theme(axis.title.y = element_text(size=8, colour="#000000")) +
#   theme(axis.title.x = element_text(size=8, colour="#000000")) +
#   ggtitle("D") +
#   theme(plot.title = element_text(hjust = 0)) +
#   scale_color_brewer(palette="Set1", name="Dataset")
# 
# pdf("rRF_size_distri.pdf")
# p_rRNA
# dev.off()

