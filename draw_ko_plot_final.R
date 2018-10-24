##Humann2 and associated scripts were used to generate functional potential data and KO conversion. Maaslin was used to identify KOs significantly
##associated with PAH exposure. The list of significant KOs were manually compiled and grouped by KO functional classes. The counts of KOs from each
##subclass were tabulated and presented as a barplot using script below:

Table <- read.tidy("pah_ko_sig_plot.txt")
Table$Class <- factor(Table$Class, levels=c("Environmental Information Processing","Genetic Information Processing","Metabolism","Mobile Genetic Element","Signalling and Cellular Processes","Unclassified"))
Table$Subclass <- factor(Table$Subclass, levels=c("Signal Transduction","DNA Repair/Recombination Protein","Folding/Sorting/Degradation","Membrane Trafficking","Replication and Repair","Transcription","Transfer RNA Biogenesis","Amino Acid/Protein Metabolism","Carbohydrate Metabolism","Cofactor/Vitamin Metabolism","Energy Metabolism","Glycan Biosynthesis/Metabolism","Lipid Metabolism","Lipopolysaccharide Metabolism","Terpenoid/Polyketide Metabolism","Xenobiotics Biodegradation","Transposase","Defense System","Mobility Protein","Secretion System","Structural Protein","Transporter Protein","Unclassified"))

colourset <- c(
"#8b8378", #Environmental Info Processing
"midnightblue","blue3","blue","cornflowerblue","DarkSlateGray3","DarkSlateGray1", #Genetic Information Processing
"#660000","#990000","#CC3333","#FF6666","#FF3333","#FF3366","#FF6699","#FF66CC","#FF33CC", #Metabolism
"#FFFF00", #Mobile Genetic Element
"#330066","#660099","#6633CC","#9933CC","#CC00FF", #Signalling and Cellular Processes
"#666666") #Unclassified

#c71585","#ff4500","#8b2500","#ff0000","#d02090","#8b2252","#8b2252","#cd2626","#f08080","#8b814c","#8b4789","#8470ff","#5d478b","#912cee","#8b8682","black")

Plot <- ggplot(Table,aes(x=PAH,fill=Subclass))
Plot <- Plot + geom_bar()
Plot <- Plot + scale_fill_manual(values=colourset)
Plot <- Plot + theme(legend.position="none")
Plot <- Plot + theme_classic() + theme(legend.position = "none",axis.text.x = element_text(angle=45,hjust=1, vjust=1))
Plot <- Plot + scale_y_continuous(expand = c(0,0))
Plot <- Plot + ylab(paste0("Significant KOs")) + xlab(paste0("PAH"))
ggsave("pah_ko_count_plot.png", width=7,height=7.units="in")


#Create separate subclass plots to add legend information to main plot generated above
#First add all colours to shades
#Environmental Information Processing
colourset_environmental <- c("#8b8378")
Plot_Environmental <- ggplot(subset(Table,Class == "Environmental Information Processing"), aes(PAH,Count,fill=Subclass)) + geom_bar(aes(fill=Subclass),stat="identity") + scale_fill_manual("Environmental Information Processing", values=colourset_environmental)
Plot_Environmental <- Plot_Environmental + theme(legend.text = element_text(size=6),legend.title = element_text(size=10))
ggsave("plot_environmental.png")

#Genetic Information Processing
colourset_genetic <- c("midnightblue","blue3","blue","cornflowerblue","DarkSlateGray3","DarkSlateGray1")
Plot_Genetic <- ggplot(subset(Table,Class == "Genetic Information Processing"), aes(PAH,Count,fill=Subclass)) + geom_bar(aes(fill=Subclass),stat="identity") + scale_fill_manual("Genetic Information Processing", values=colourset_genetic)
Plot_Genetic <- Plot_Genetic + theme(legend.text = element_text(size=6),legend.title = element_text(size=10))
ggsave("plot_genetic.png")

#Metabolism
colourset_metabolism <- c("#660000","#990000","#CC3333","#FF6666","#FF3333","#FF3366","#FF6699","#FF66CC","#FF33CC")
Plot_Metabolism <- ggplot(subset(Table,Class == "Metabolism"), aes(PAH,Count,fill=Subclass)) + geom_bar(aes(fill=Subclass),stat="identity") + scale_fill_manual("Metabolism", values=colourset_metabolism)
Plot_Metabolism <- Plot_Metabolism + theme(legend.text = element_text(size=6),legend.title = element_text(size=10))
ggsave("plot_metabolism.png")

#Mobile Genetic Element
colourset_mge <- c("#FFFF00")
Plot_MGE <- ggplot(subset(Table,Class == "Mobile Genetic Element"), aes(PAH,Count,fill=Subclass)) + geom_bar(aes(fill=Subclass),stat="identity") + scale_fill_manual("Mobile Genetic Element", values=colourset_mge)
Plot_MGE <- Plot_MGE + theme(legend.text = element_text(size=6),legend.title = element_text(size=10))
ggsave("plot_mge.png")

#Signalling
colourset_signal <- c("#330066","#660099","#6633CC","#9933CC","#CC00FF")
Plot_Signal <- ggplot(subset(Table,Class == "Signalling and Cellular Processes"), aes(PAH,Count,fill=Subclass)) + geom_bar(aes(fill=Subclass),stat="identity") + scale_fill_manual("Signalling and Cellular Processes", values=colourset_signal)
Plot_Signal <- Plot_Signal + theme(legend.text = element_text(size=6),legend.title = element_text(size=10))
ggsave("plot_signal.png")

#Unclassified
Plot_Unclassified <- ggplot(subset(Table,Class == "Unclassified"), aes(PAH,Count,fill=Subclass)) + geom_bar(aes(fill=Subclass),stat="identity") + scale_fill_manual("Unclassified", values="#666666")
Plot_Unclassified <- Plot_Unclassified + theme(legend.text = element_text(size=6),legend.title = element_text(size=10))
ggsave("plot_unclassified.png")

#Build plots
environmental <- ggplot_gtable(ggplot_build(Plot_Environmental))
genetic <- ggplot_gtable(ggplot_build(Plot_Genetic))
metabolism <- ggplot_gtable(ggplot_build(Plot_Metabolism))
mge <- ggplot_gtable(ggplot_build(Plot_MGE))
signalling <- ggplot_gtable(ggplot_build(Plot_Signal))
unclassified <- ggplot_gtable(ggplot_build(Plot_Unclassified))

grob_environmental<-environmental$grobs[[15]]
grob_genetic<-genetic$grobs[[15]]
grob_metabolism <- metabolism$grobs[[15]]
grob_mge <- mge$grobs[[15]]
grob_signalling <- signalling$grobs[[15]]
grob_unclassified <- unclassified$grobs[[15]]

grid.arrange(arrangeGrob(Plot,arrangeGrob(grob_environmental,grob_genetic,grob_metabolism,grob_mge,grob_signalling,grob_unclassified), ncol = 2,widths=c(3/4,1/4)))
