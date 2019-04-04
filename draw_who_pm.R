Table <- read.tidy("worldwide_PAH.txt")
Table$Country <- factor(Table$Country, levels = c("Baoding, China","Dalian, China","Australia","Canada","China","France","Japan","UK","USA"))

Plot <- ggplot(Table,aes(x=Ten,y=Two,colour=Country))
Plot <- Plot + geom_point(size=1) + theme_classic() + scale_fill_brewer(palette="Set3")
Plot <- Plot + xlab(paste0("[PM10] (ug/m3)")) + ylab(paste0("[PM10] (ug/m3)")) + theme(legend.title=element_blank())
ggsave("WHO_data.pdf",useDingbats=FALSE)