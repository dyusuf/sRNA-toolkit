setwd("/home/dilmurat/Work/combined-experiments/")
library(plyr)
library(ggplot2)

all_players <- read.csv("2.sample_distri-all_players.csv", header=T, sep="\t")

plasma_players <- subset(all_players, condition=="p" & rank1 == "Homo_sapiens")

player_summery <- ddply(plasma_players, 
                        c("rna", "sequence", "sequence_length", "GC", "RNA", "condition", "key_player"),
                        summarise, 
                        abundance_average=mean(abundance), 
                        sd=sd(abundance, na.rm=TRUE),
                        cv=sd/abundance_average)
player_summery <- ddply(player_summery,  c("condition"), transform,
                        fraction=abundance_average/sum(abundance_average))

# player_summery <- ddply(plasma_players, c( "rna", "sequence", "sequence_length", "GC", "RNA", "condition", "key_player"), 
#                         summarise, 
#                         abundance_average=mean(abundance)
# )

levels(player_summery$RNA) <- c(levels(player_summery$RNA), "Other RNA classes")
for (i in 1:length(player_summery$RNA)) {
  if (!(player_summery[i,]$RNA %in% c("miRNA", "Y_RNA", "rRNA", "tRNA", "snRNA", "snoRNA", "lncRNA", "miscRNA", "piRNA", "mRNA"))) 
    player_summery[i,]$RNA <-  factor("Other RNA classes")
}

player_summery$RNA <- as.factor(sapply(player_summery$RNA, gsub, pattern="Y_RNA", replacement="Y RNA"))

player_summery$RNA <- factor(player_summery$RNA, levels=c(
  "miRNA", "Y RNA", "rRNA",   "tRNA", "snoRNA", "snRNA",  "lncRNA", "piRNA", "miscRNA", "mRNA", "Other RNA classes"
))

player_summery <- subset(player_summery, sequence_length>17)

top_mirs_cohort <- subset(player_summery, RNA == "miRNA" & abundance_average >= 1000 & sequence_length>=19)
mir_seqs <- unique(as.vector(top_mirs_cohort$sequence))

interesting_ones <- as.vector(arrange(subset(player_summery, sequence_length==21), desc(abundance_average))[1,]$sequence)
interesting_ones <- c(as.vector(arrange(subset(player_summery, sequence_length==22), desc(abundance_average))[1,]$sequence), interesting_ones)
interesting_ones <- c(as.vector(arrange(subset(player_summery, sequence_length==23), desc(abundance_average))[1,]$sequence), interesting_ones)
interesting_ones <- c(as.vector(arrange(subset(all_players, sequence_length==23 & RNA == "Y_RNA" & condition == "CD4H24h"), desc(abundance))[1,]$sequence), interesting_ones)
interesting_ones <- c(as.vector(arrange(subset(player_summery, sequence_length==27 & RNA == "Y RNA"), desc(abundance_average))[1,]$sequence), interesting_ones)
interesting_ones <- c(as.vector(arrange(subset(player_summery, sequence_length==31 & RNA == "Y RNA"), desc(abundance_average))[1,]$sequence),
                      as.vector(arrange(subset(player_summery, sequence_length==32 & RNA == "Y RNA"), desc(abundance_average))[1,]$sequence),
                      as.vector(arrange(subset(player_summery, sequence_length==33 & RNA == "Y RNA"), desc(abundance_average))[1,]$sequence),
                      interesting_ones)
interesting_ones <- c(as.vector(arrange(subset(player_summery, sequence_length==35 & RNA == "rRNA"), desc(abundance_average))[1,]$sequence),
                      interesting_ones)
interesting_ones <- c(as.vector(arrange(subset(player_summery, sequence_length==32 & RNA == "tRNA"), desc(abundance_average))[1,]$sequence),
                      as.vector(arrange(subset(player_summery, sequence_length==33 & RNA == "tRNA"), desc(abundance_average))[1,]$sequence),
                      interesting_ones)

interest_seqs <- unique(c(mir_seqs, interesting_ones))

top_series <- subset(all_players, rank1 =="Homo_sapiens" & condition %in% c("plasma", "exosome145"))
top_series$rank1 <- droplevels(top_series$rank1)
top_series$condition <- droplevels(top_series$condition)
top_series$RNA <- as.factor(sapply(top_series$RNA, gsub, pattern="Y_RNA", replacement="Y RNA"))

#select subset
key_seqs <- subset(top_series, sequence %in% interest_seqs)
#replace sample_ids with proper time points
key_seqs["individual_id"] <- "NA"
key_seqs["time_point"] <- "NA"
for (i in 1:length(key_seqs$sample_id)) {key_seqs[i,]$individual_id <- substr(key_seqs[i,]$sample_id, 1, 3)}
for (i in 1:length(key_seqs$sample_id)) {key_seqs[i,]$time_point <- substr(key_seqs[i,]$sample_id, 4, 5)}

T1 <- ISOdate(2012, 07, 09, 18, 05, tz = "")
T2 <- T1 + 2*3600 + 5*60
T3 <- T2 + 4*3600 + 7*3600 + 30*60
T4 <- T3 + 3*3600 - 10*60
T5 <- T4 + 1.5 * 3600
T6 <- T5 + 3*3600
T7 <- T6 + 3*3600
T8 <- T7 + 6*3600 + 8*3600

key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T1", replacement=T1))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T2", replacement=T2))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T3", replacement=T3))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T4", replacement=T4))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T5", replacement=T5))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T6", replacement=T6))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T7", replacement=T7))
key_seqs$time_point <- as.factor(sapply(key_seqs$time_point, gsub, pattern="T8", replacement=T8))
key_seqs$time_point <- as.POSIXct(key_seqs$time_point)

#the 11 sequences, associated with peak sequence, shown in expression pattern plot
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_1001", replacement="miR-486-5p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_2345", replacement="miR-486-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_950", replacement="miR-486-5p-3'~A"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="yRF_4412", replacement="yRF-Y4-3p-3'L4"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="yRF_118", replacement="yRF-Y4-3p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="yRF_177", replacement="yRF-Y4-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="yRF_3095", replacement="yRF-Y4-5p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="yRF_1632", replacement="yRF-Y4-5p-3'R1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="rRF_356", replacement="rRF-RNA28S5"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="tRF_167", replacement="tRF-tRNA-Gly(GCC/CCC)"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="tRF_1663", replacement="tRF-tRNA-Val(CAC/AAC)"))
#the 33 sequences selected with abundance >= 100 normalized read count
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-103-2~mir-107_4069", replacement="miR-103a/107-3p-3'L4"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-10a_2816", replacement="miR-10a-5p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-10b_1671", replacement="miR-10b-5p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-126_1910", replacement="miR-126-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-142_1000", replacement="miR-142-5p-5'L2-3'L3"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-143_3202", replacement="miR-143-3p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-15a_4194", replacement="miR-15a-5p-3'L2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-16-2_2374", replacement="miR-16-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-182_1810", replacement="miR-182-5p-3'L2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-191_1434", replacement="miR-191-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-192_3998", replacement="miR-192-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-21_117", replacement="miR-21-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-22_272", replacement="miR-22-3p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-22_3185", replacement="miR-22-3p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-25_2636", replacement="miR-25-3p-3'L2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-25_335", replacement="miR-25-3p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-27b_2362", replacement="miR-27b-3p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-30e_3214", replacement="miR-30e-5p-3'R2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-423_141", replacement="miR-423-5p-3'L2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-423_1960", replacement="miR-423-5p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-423_1996", replacement="miR-423-5p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-451_1655", replacement="miR-451a-3'L2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-451_2898", replacement="miR-451a-3'R1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-451_4357", replacement="miR-451a"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-451_525", replacement="miR-451a-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_1177", replacement="miR-486-5p-3'~U"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_2004", replacement="miR-486-5p-3'~AA"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_2191", replacement="miR-486-5p-3'L1~A"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_826", replacement="miR-486-5p-3'L2"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-486_945", replacement="miR-486-5p-3'R1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-92-2~mir-92-1_1449", replacement="miR-92a-3p-3'L1"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-92-2~mir-92-1_2042", replacement="miR-92a-3p"))
key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-92-2~mir-92-1_4484", replacement="miR-92a-3p-3'L1~A"))
# key_seqs$rna <- as.factor(sapply(key_seqs$rna, gsub, pattern="mir-181a-2_86", replacement="miR-181a-5p"))

key_seqs_series = subset(key_seqs, condition=="plasma")
long_form <- key_seqs_series[,c(4,5,8)]
wide_form <- reshape(long_form, idvar = c("sample_id"), timevar = "rna", direction = "wide")
wide_form <- wide_form[2:length(wide_form)]
wide_form[is.na(wide_form)] <- 0
names(wide_form) <- substring(names(wide_form), first=11, last=100)
key_cor <- cor(wide_form)

key_seqs$sample_id <-  as.factor(toupper(key_seqs$sample_id))
long_form <- key_seqs[,c(3,4,5,6,8)]
wide_form <- reshape(long_form, idvar = c("sequence", "rna",  "sample_id"), timevar = "condition", direction = "wide")
wide_form[is.na(wide_form)] <- 0
require(plyr)
func <- function(wide_form)
{
  return(data.frame(COR = cor(wide_form$abundance.plasma, wide_form$abundance.exosome145),
                    plasma_abundance=sum(wide_form$abundance.plasma)
  ))
}
cor_exo_plasma <- ddply(wide_form, .(rna, sequence), func)
cor_exo_plasma$condition <- "plasma vs exosome"
cor_exo_plasma <- ddply(cor_exo_plasma,  c("condition"), transform, 
                        plasma_percentage=plasma_abundance/sum(plasma_abundance)
)


library(reshape2)
library(grid)
library(ggplot2)

po.nopanel <- list(theme(panel.background = element_blank(), panel.grid.minor = element_blank(), 
                         panel.grid.major = element_blank()))
# plot XY quantative fill
                   
key_cor_long <- melt(key_cor)
hc <- hclust(dist(key_cor))
key_cor_long$Var1 <- factor(key_cor_long$Var1, levels=hc$labels[c(hc$order)])
key_cor_long$Var2 <- factor(key_cor_long$Var2, levels=rev(hc$labels[c(hc$order)]))

pxy <- ggplot(key_cor_long, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value), colour = "white") + 
#   scale_fill_gradient(low = "midnightblue", high = "green2",
#                                           name="Pearson\ncorrelation (r)",
#                             guide = guide_colorbar(direction = "horizontal",
#                                                  title.position = "top",
#                                                  label.position="bottom")) +
  scale_fill_gradient(low = "midnightblue", high = "green2", name="Pearson\ncorrelation\n(r)") +
  theme(axis.title = element_blank()) + 
  theme(axis.text.x = element_text(size=6, colour="#000000", angle=90,  hjust=1)) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(axis.ticks=element_blank())     +   
  theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
  theme(plot.margin=unit(c(0.5,-0.2,0.6,-0.3), "cm")) +                   
#                      scale_x_discrete(breaks = NULL) +
#                      scale_y_discrete(breaks = NULL) +                   
  po.nopanel        +
  ggtitle("A")   +
  theme(plot.title = element_text(hjust = 0)) 

pxy_x <- ggplot(key_cor_long, aes(Var1, y=-10000)) + 
                     geom_tile() + 
                     theme(axis.title = element_blank()) + 
                     theme(axis.ticks=element_blank())     +     
                     scale_y_discrete(breaks = NULL) + 
                     theme(rect = element_blank()) +
                     theme(line = element_blank()) +
                     theme(plot.margin=unit(c(0.2,-0.2,0.5,-0.3), "cm")) +
                     theme(axis.text.x = element_text(size=6, colour="#000000", angle=90, vjust=0.5,hjust=1)) 
pxy_y <- ggplot(key_cor_long, aes( x=-10000, Var2)) + 
  geom_tile() + 
  theme(axis.title = element_blank()) + 
  theme(axis.ticks=element_blank())     +     
  scale_x_discrete(breaks = NULL) + 
  theme(rect = element_blank()) +
  theme(line = element_blank()) +
  theme(plot.margin=unit(c(0.5,0,0.6,0.5), "cm")) +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  ggtitle("A")   +
  theme(plot.title = element_text(hjust = 1, colour="white")) 

cor_exo_plasma$rna <- factor(cor_exo_plasma$rna, levels=rev(hc$labels[c(hc$order)]))
ploty <- ggplot(cor_exo_plasma, aes( condition, rna)) + geom_tile(aes(fill = COR), 
          colour = "white") + scale_fill_gradient(low = "midnightblue", high = "green2", name="Pearson\ncorrelation\n(r)") +
          theme(legend.text=element_text(size=6), legend.title = element_text(size=6.5)) +
         theme( axis.title = element_blank()) + 
#   theme(axis.text.y = element_text(size=7, colour="#000000", hjust=1), axis.ticks.y=element_blank()) +
#                      scale_x_discrete(breaks = NULL) +
                     scale_y_discrete(breaks = NULL) +
                     theme(plot.margin=unit(c(0.5,0.5,0.6,-0.2), "cm")) + 
                     po.nopanel        +
                     ggtitle("B") +
                    theme(plot.title = element_text(hjust = 0)) +
                    theme(axis.text.x = element_text(size=7, colour="#000000")) +
                    theme(axis.ticks=element_blank()) 

grob_ploty <- ggplotGrob(ploty)
ploty_x <- grob_ploty$grobs[[5]]
ploty_x$vp$x <-  unit(0.81, 'npc')
ploty_x$vp$y <-  unit(0.2, 'npc')

# Define layout for the plots (2 rows, 2 columns)
layt <- grid.layout(nrow = 2, ncol = 4, heights = c(7/8, 1/8), widths = c(1/8, 5/8, 1/8, 
                                                                          1/8), default.units = c("null", "null"))
# View the layout of plots
# grid.show.layout(layt)

tmp <- ggplotGrob(pxy)
legend <- tmp$grobs[[8]]
legend$vp$x <-  unit(0.9, 'npc')
legend$vp$y <-  unit(0.6, 'npc')

pdf("key_seq_corM_corPE-20150310.pdf", width=10, height=7)
# Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout = layt))
print(ploty+theme(legend.position="none")+ scale_x_discrete(breaks = NULL), vp = viewport(layout.pos.row = 1, layout.pos.col = 3))
print(pxy+theme(legend.position="none", axis.text = element_blank()), vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
print(pxy_x, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
print(pxy_y, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
grid.draw(ploty_x)
grid.draw(legend)
dev.off()



###################################################################################
key_seqs$individual_id = toupper(key_seqs$individual_id)

key_seqs$condition <- as.factor(sapply(key_seqs$condition, gsub, pattern="exosome145", replacement="exosome fraction"))


indvid_cv <- ddply(key_seqs, 
                c("rna", "sequence",  "individual_id", "condition"),
                summarise, 
                abundance_average=mean(abundance), 
                sd=sd(abundance, na.rm=TRUE),
                cv=sd/abundance_average)

indvid_cv$condition <- factor(indvid_cv$condition, levels=c("plasma", "exosome fraction"))

p_indi <- ggplot(indvid_cv, aes(x=cv, color=individual_id)) + 
  stat_ecdf(size=1, alpha=0.8) + 
  scale_color_brewer(palette="Set1", name="Individual ID") +
  ylab("Empirical cumulative distribution function") +
  xlab("Coefficient of variation")+
  theme_minimal() +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(axis.text.x = element_text(size=6, colour="#000000")) +
  theme(legend.text = element_text(size=6), legend.title=element_text(size=6.5)) +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~condition, nrow=2) +
  ggtitle("B") +
  theme(plot.title = element_text(hjust = 0)) 


cluster1 <- c(
              "miR-22-3p",
              "miR-25-3p",
              "miR-423-5p",
              "miR-142-5p-5'L2-3'L3",
              "miR-192-5p",
              "miR-30e-5p-3'R2"
              )
cluster3 <- c(
              'yRF-Y4-5p',
              'miR-191-5p',
              "miR-27b-3p-3'L1",
              'miR-126-5p',
              'miR-21-5p',
              "miR-143-3p-3'L1",
              "miR-10a-5p-3'L1",
              'tRF-tRNA-Val(CAC/AAC)',
              'tRF-tRNA-Gly(GCC/CCC)'
              )
cluster2 <- c(
              'miR-486-5p',
              "miR-486-5p-3'L1",
              "miR-486-5p-3'~A",
              "miR-486-5p-3'~U",
              "miR-486-5p-3'L2",
              "miR-486-5p-3'L1~A",
              "miR-486-5p-3'~AA",
              "miR-486-5p-3'R1",
              'miR-92a-3p',
              "miR-92a-3p-3'L1",
              "miR-92a-3p-3'L1~A",
              'miR-16-5p',
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
              'yRF-Y4-3p',
              "yRF-Y4-3p-3'L4",
              "miR-10b-5p-3'L1",
              'rRF-RNA28S5'
              )

time_cv <- ddply(key_seqs, 
                 c("rna", "sequence",  "time_point", "condition"),
                 summarise, 
                 abundance_average=mean(abundance), 
                 sd=sd(abundance, na.rm=TRUE),
                 cv=sd/abundance_average)
time_cv$condition <- factor(time_cv$condition, levels=c("plasma", "exosome fraction"))

time_cv["cluster"] <- "NA"
for (i in 1:length(time_cv$rna)) {if (time_cv[i,]$rna %in% cluster1) 
                                    time_cv[i,]$cluster <- "cluster 1"
                                  else if (time_cv[i,]$rna %in% cluster2) 
                                    time_cv[i,]$cluster <- "cluster 2"
                                  else if (time_cv[i,]$rna %in% cluster3) 
                                      time_cv[i,]$cluster <- "cluster 3"
                                  else time_cv[i,]$cluster <- "cluster NA"
                                  }


time_cv$time_point <- as.character(time_cv$time_point)
two_points_cv <- subset(time_cv, grepl("10:30", time_point) | grepl("12:00", time_point))
two_points_cv$time_point <- gsub("2012-07-10 ", "", two_points_cv$time_point) 
day_cv <- subset(time_cv, grepl("07-10", time_point))
day_cv$time_point <- gsub("2012-07-10 ", "", day_cv$time_point) 

day_cv_cluster1 <- subset(day_cv, cluster=="cluster 1")
day_cv_cluster2 <- subset(day_cv, cluster=="cluster 2")
day_cv_cluster3 <- subset(day_cv, cluster=="cluster 3")
day_cv_cluster1_3 <- subset(day_cv, cluster %in% c("cluster 3", "cluster 1"))

time_cv$time_point <- gsub("2012-07-09", "day 1", time_cv$time_point) 
time_cv$time_point <- gsub("2012-07-10", "day 2", time_cv$time_point) 
time_cv$time_point <- gsub("2012-07-11", "day 3", time_cv$time_point) 

time_cv$time_point <- substr(time_cv$time_point, 1, 11)

p_time <- ggplot(time_cv, aes(x=cv, color=time_point)) + 
  stat_ecdf(size=1, alpha=0.8) + 
  scale_color_brewer(palette="Paired", name="Time point") +
  ylab("Empirical cumulative distribution function") +
  xlab("Coefficient of variation")+
  theme_minimal() +
  theme(axis.text.y = element_text(size=6, colour="#000000")) +
  theme(axis.text.x = element_text(size=6, colour="#000000")) +
  theme(legend.text = element_text(size=6), legend.title=element_text(size=6.5)) +
  theme(axis.title = element_text(size=8, colour="#000000")) +
  facet_wrap(~condition, nrow=2) +
  ggtitle("A") +
  theme(plot.title = element_text(hjust = 0)) 

# Define layout for the plots (2 rows, 2 columns)
layt <- grid.layout(nrow = 1, ncol = 2, heights = c(4/8, 4/8), widths = c(4/8, 
                                                                          4/8), default.units = c("null", "null"))
# View the layout of plots
# grid.show.layout(layt)

pdf("fig4_inter-intra-var.pdf", width=10, height=7)
# Draw plots one by one in their positions
grid.newpage()
pushViewport(viewport(layout = layt))
print(p_time, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p_indi, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

