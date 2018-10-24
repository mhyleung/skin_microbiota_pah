##AIC model data generated from sloan script. AIC data for neutral/binomial/poisson models presented. Data to be used below:
Table <- read.tidy("aic_models.txt")
Table$Model <- factor(Table$Model,levels=c("Neutral","Binomial","Poisson"))
Plot <- ggplot(Table,aes(x=Sample,y=AIC,fill=Model))
Plot <- Plot + geom_bar(stat="identity",position="dodge", width=0.5)
Plot <- Plot + xlab(paste0("Sample Group")) + ylab(paste0("Akaike Information Criterion Score"))
Plot <- Plot + theme_classic() + theme(legend.title=element_blank())
Plot <- Plot + scale_fill_manual(values=c("#cd8c95","#6e7b8b","#cdb38b"))
ggsave("aig.pdf",width=6.48,height=5,units="in")


##Group OTUs into genera, and tabulate whether OTUs are within/above/below Sloan model
Table <- read.tidy("sloan_by_common_genus.txt")
Table$Genus <- factor(Table$Genus, levels=c("Overall Bacteria","Propionibacterium","Staphylococcus","Enhydrobacter","Paracoccus","Corynebacterium", "Overall Fungi", "Malassezia"))
Plot <- ggplot(Table,aes(x=Genus, y=Percentage, fill=Partition))
Plot <- Plot + geom_bar(stat="identity") + facet_grid(Site~City,scales="free")
Plot <- Plot + scale_fill_brewer(palette = "Dark2")
Plot <- Plot + ylab(paste0("% of OTUs")) + theme(axis.title.x = element_blank())
Plot <- Plot + theme_classic() + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1)) + theme(legend.position="bottom")
Plot <- Plot + theme(legend.title = element_blank(),axis.title.x=element_blank()) + scale_y_continuous(expand = c(0, 0))
ggsave("sloan_by_genus.pdf",width=5,height=5, units="in",useDingbats=FALSE)

