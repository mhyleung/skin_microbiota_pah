#require("breakaway")
library(breakaway)
library(MASS)

#Use cast OTU table as input, and open metadata
data <- read.tidy("../otu_SILVA_revise/OTU_table_clean5_wt_control_400_490_65.tidy_fixed_cast.txt")
metadata <- read.tidy("../meta/meta_w_dandruff_wt_400_490_65_9C_209C.txt")

#Create frequency table
frequencytablelist <- lapply(apply(data,2,table),as.data.frame)
frequencytablelist <- lapply(frequencytablelist,function(x) x[x[,1]!=0,])

#Test breakaway on frequency count tables to see frequency table has been compiled successfully
breakaway(frequencytablelist[[2]])

#estimate shannon diversity
shannon(frequencytablelist[[2]])
#briefly look at variability of repeated testings
resample_estimate(data[,2], shannon)
resample_estimate(data[,2], shannon)
resample_estimate(data[,2], shannon)
#Compute and compare Shannon evenness estimates by iterating Shannon diversity estimates on each sample of dataset (999 iterations per sample to account for variability), include in data the 25% and 75% quantiles

ns <- unlist(lapply(frequencytablelist, function(x) sum(x[,2])))
estimates_shannon <- matrix(NA,nrow=dim(data)[2],ncol=4)
rownames(estimates_shannon) <- colnames(data)
colnames(estimates_shannon) <- c("shannon_est","shannon_seest","shannon_lcb","shannon_ucb")
for (i in 2:dim(data)[2]) {
  resample_estimate(data[,i], shannon, my_sample_size = ns)
  samples <- replicate(999, resample_estimate(data[,i], shannon, my_sample_size = ns))
  estimates_shannon[i,1] <- mean(samples)
  estimates_shannon[i,2] <- sd(samples)
  estimates_shannon[i,3:4] <- quantile(samples, c(0.025, 0.975))
}

write.tidy(estimates_shannon,"estimates_shannon.txt",row.names=TRUE)
