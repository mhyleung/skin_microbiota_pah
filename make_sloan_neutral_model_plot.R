library(readr)
library(dplyr)
library(ggplot2)

skin <- read_tsv("predictions/prediction_dalian_scalp_10452.txt")
names(skin)[1] <- c("OTU")
skin <- skin %>% mutate(Partition=ifelse(freq < pred.lwr, "Below", "Neutral"))
skin <- skin %>% mutate(Partition=ifelse(freq > pred.upr, "Above", Partition))
write_tsv(skin,"prediction_dalian_scalp_10452.txt")

#spp <- read_tsv("rarefied_otu_table_individual_1086_cast/OTU_table_clean_1086_WKS_3Y.cast.txt")
spp <- read_tsv("OTU_table_clean5_wt_control_400_490_65.tidy_fixed_cast_10452_Dalian_Scalp.txt")
N <- mean(apply(spp, 1, sum))
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N
p <- as.data.frame(p)
colnames(p) <- c("mean_abundance")
write.table(p, "OTU_mean_abundance_dalian_scalp_10452.txt", sep="\t", row.names=TRUE)
#write.table(p, "rarefied_otu_table_indivdiual_otu_mean_abundance/OTU_mean_abundance_WKS_3Y.txt", sep="\t", row.names=TRUE)
#Please change the colnames in abundance manually
#merge mean_abundance with sloan prediciton results
abundance <- read_tsv("OTU_mean_abundance_dalian_scalp_10452.txt")
names(p)[1] <- c("OTU")
names(p)[2] <- c("Abundance")
merge <- left_join(skin, abundance)
write_tsv(merge, "prediction_dalian_scalp_10452.txt")

#Join OTU with taxonomy
Table <- read.tidy("prediction_dalian_scalp_10452.txt")
Tax <- read.tidy("otu_taxonomy.txt")
Merge <- merge(Table,Tax,by="OTU",all.X=TRUE)
write.tidy(Merge,"predictions/prediction_dalian_scalp_10452_wtax.txt")

#Plot the Sloan neutral model for all individuals
Table <- read.tidy("predictions/prediction_dalian_cheek_10452_wtax.txt")
Plot <- ggplot(Table,aes(x=log(Abundance)))
Plot <- Plot + geom_point(aes(y=freq, colour=Partition), size=0.5)
Plot <- Plot + geom_line(aes(y=freq.pred))
Plot <- Plot + geom_line(aes(y=pred.lwr), linetype="dotted")
Plot <- Plot + geom_line(aes(y=pred.upr), linetype="dotted")
Plot <- Plot + scale_color_brewer(palette = "Dark2")
Plot <- Plot + theme_classic()
Plot <- Plot + theme(legend.title = element_blank())
Plot <- Plot + theme(plot.title = element_text(hjust = 0.5, face="bold"))
Plot <- Plot + ggtitle("Dalian Scalp")
Plot <- Plot + xlab(paste0("log (Mean Relative Abundance)")) +ylab(paste0("Occurrence Frequency"))
ggsave("sloan_neutral_model_dalian_cheek_10452.pdf",width=4.23,height=3.25,units="in",useDingbats=FALSE)
