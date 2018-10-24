##Taxonomic analysis and heatmap was performed and created using metaphlan and associated hclust.py script. taxonomic bar plot as presented on
##manuscript was created using script below. Merged metaphlan table of species-level abundnace data was used as input:

Table <- read.tidy("merged_metaphlan_table.txt")
Table$Species<- factor(Table$Species,levels=c("Propionibacterium acnes","Betapapillomavirus 3","Propionibacterium phage P101A","Propionibacterium phage P101D","Dasheen mosaic virus","Paracoccus denitrificans","Enhydrobacter aerosaccus","Corynebacterium lipophiloflavum","Paracoccus sp.","Minor/Unclassified"))
Table$Sample<- factor(Table$Sample,levels=c("4C","50C","59C","75C","77C","4S","50S","59S","75S","77S","205C","213C","266C","299C","205S","213S"))

Plot <- ggplot(Table,aes(x=Sample, y=Normalized,fill=Species))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + scale_fill_brewer(palette="Set3")
Plot <- Plot + theme_classic()
Plot <- Plot + ylab(paste0("Within-Domain Relative Abundance (%)")) + facet_wrap(~Domain)
Plot <- Plot + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1))
Plot <- Plot + scale_y_continuous(expand = c(0,0))
Plot <- Plot + theme(legend.position="bottom",legend.title=element_blank(),axis.title.x=element_blank(),plot.title = element_text(hjust = 0.5))
ggsave("metaphlan_relab_new.pdf",width=9.6,height=4.49,units="in",useDingbats=FALSE)
