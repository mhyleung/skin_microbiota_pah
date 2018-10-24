##Significant relationships between PAH exposure and species abundance were plotted. Input file contains species abundance in %
##and PAH exposure for each sample.

Table <- read.tidy("skin_bacteria_pah.txt")

Table$Species <- factor(Table$Species, levels = c("S. epidermidis","Alicycliphilus sp.", "Neisseria sp.", "K. sedentarius", "M. luteus"))

Plot <- ggplot(Table,aes(x=PAH_conc, y=Species_conc, colour=Species))
Plot <- Plot + geom_point()
Plot <- Plot + facet_wrap(~PAH,scales="free", nrow=2)
Plot <- Plot + theme_classic()
Plot <- Plot + xlab(paste0("[PAH]")) + ylab(paste0("Species Relative Abundance (%)"))
ggsave("skin_species_pah_no_trend_line.pdf",width=9.37,height=4.86,units="in",useDingbats=FALSE)

