##Natural connectivity reduction presented as the percentage reduction of the original natural connectivity upon removal of each node. Input data
##from spiec_easi directed attack scripts generate absolute natural connecfivity figures, which were manually converted to percentages and saves as
##a txt file to be used as input below:

#Open attack distribution table
Table <- read.tidy("Attack_Combined_Percentage.txt")
Table$Group <- factor(Table$Group, levels=c("Baoding Cheek Healthy","Baoding Cheek Acne","Dalian Cheek Healthy","Dalian Cheek Acne"))

Plot <- ggplot(Table, aes(x=Position_Perc,y=Percentage,colour=Group))
Plot <- Plot + geom_line()
Plot <- Plot + scale_x_continuous(limits = c(0,75))
Plot <- Plot + scale_color_manual(values=c("blue","green","red","black"))
Plot <- Plot + theme_classic()
Plot <- Plot + ylab(paste0("% Reduction in Natural Connectivity")) + xlab(paste0("% of Nodes Removed"))
Plot <- Plot + theme(legend.title = element_blank(), axis.text = element_text(size=12), axis.title=element_text(size=14))
ggsave("reduction_connectivity_cheek.pdf",width=7, height=7,units="in")
