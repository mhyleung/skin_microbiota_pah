##FindFungi results were arranged in cast format as in metaphlan2 output, and was used as input for script below

Table <- read.tidy("findfungi_output_plot.txt")
Table <- melt(Table)
Table$Species<- factor(Table$Species,levels=c("Malassezia restricta","Malassezia globosa","Melampsora pinitorqua","Malassezia yamatoensis","Aspergillus sp.","Malassezia slooffiae","Minor/Unclassified"))
Table$variable<- factor(Table$variable,levels=c("4C","50C","59C","75C","77C","4S","50S","59S","75S","77S","205C","213C","266C","299C","205S","213S"))

Plot <- ggplot(Table,aes(x=variable, y=value,fill=Species))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + scale_fill_brewer(palette="Set3")
Plot <- Plot + theme_classic()
Plot <- Plot + ylab(paste0("Within-Domain Relative Abundance (%)"))
Plot <- Plot + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
Plot <- Plot + scale_y_continuous(expand = c(0,0))
Plot <- Plot + theme(legend.position="bottom",legend.title=element_blank(),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5))
ggsave("findfungi_relab_fungus.pdf",width=9.6,height=4.49,units="in",useDingbats=FALSE)
