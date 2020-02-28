#Plot Coefficient bar chart
Table <- read.tidy("humann2_subclass_pathways_sign_count.txt")
Table$PAH <- factor(Table$PAH, levels=c("Benzo[k]fluoranthene","Acenaphthylene","Fluorene","Dibenz[a,h]anthracene","Benzo[a]pyrene","Phenanthrene","Pyrene","Benzo[a]anthracene","Chrysene","Indeno[1,2,3-c,d]pyrene",
                                          "Anthracene","Benzo[g,h,i]perylene","Acenaphthene","Benzo[b]fluoranthene","Fluoranthene"))
Table$Sign <- factor(Table$Sign, levels=c("Positive","Negative"))
Plot <- ggplot(Table,aes(x=PAH,y=Count,fill=Sign))
Plot <- Plot + geom_bar(stat="identity")
Plot <- Plot + coord_flip() + theme_classic() + scale_fill_brewer(palette="Set2") + theme(legend.position="none")
ggsave("pah_ko_number_new.pdf",width=6,height=8,units="in")
