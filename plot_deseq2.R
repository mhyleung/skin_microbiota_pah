##DeSeq2 was performed using QIIME1's "compare_categories.py" script. The output was used as input on this script to generate plot

#Merge with taxonomy (do separately for bacteria and fungi)
Table <- read.tidy("Deseq2_Combined.txt")
Tax <- read.tidy("../spiec/taxatableall.txt")
Merge <- merge(Table,Tax,by="OTU",all.X=TRUE)
write.tidy(Merge,"Deseq2_Combined_w_tax.txt")


#Open data with a combined table of both bacterial and fungal data
Table <- read.tidy("Deseq2_Combined_w_tax_cross_domain.txt")
Table$Phenotype <- factor(Table$Phenotype,levels=c("Healthy","Acne","Dandruff"))
Plot <- ggplot(Table,aes(x=reorder(Name,log2FoldChange),y=log2FoldChange, fill=City))
Plot <- Plot + geom_bar(stat="identity",width=0.4)
Plot <- Plot + facet_wrap(Site~Phenotype,nrow=1)
Plot <- Plot + theme_classic() + coord_flip() + scale_fill_manual(values=c("#9243e7","#73f54a"))
Plot <- Plot + xlab(paste0("OTU and Genus")) + ylab("DeSeq2 Log-Fold Differential")
Plot <- Plot + theme(legend.psosition = "none",axis.text.y = element_text(size=6))
Plot <- Plot + geom_hline(yintercept = 0, linetype="dotted")
ggsave("deseq_both_sites.pdf",width=6,height=6)
