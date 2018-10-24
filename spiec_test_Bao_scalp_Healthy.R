#libraries
library(readr)
library(devtools)
library(SpiecEasi)
library(phyloseq)
library(tibble)
library(igraph)
library(stringr)
library(tidyr)
library(ggplot2)
library(Matrix)
library(orca)
library(ggnetwork)
#install_github('wilkox/forcats')
library(wilkoxmisc)
library(forcats)
library(intergraph)
library(network)
library(plyr)
library(dplyr)
library(tidyverse)

#Filter OTUs from table where prevalence is below 25% of dataset (bacteria)
OTUTableBCH <- read_tsv("OTU_Bacterial_Baoding_Scalp_Healthy.txt")
# cast OTU table
OTUTableBCH <- dcast(OTUTableBCH, Sample ~ OTU, value.var = "Count", fill = 0)
OTUTableBCH[1] <- NULL
# select most prevalent OTUs by occurrence frequency
OTUTableBCH <- 1*(OTUTableBCH>0)
freqBCH <- apply(OTUTableBCH, 2, mean)
freqBCH <- freqBCH[freqBCH >=0.25]
freqBCH <- as.matrix(freqBCH)
abundOTUsBCH <- rownames(freqBCH)
prevalenceBCH <- freqBCH[,1]


#calculate the abundance of prevalent OTUs (bacteria)
OTUTableBCH <- read_tsv("OTU_Bacterial_Baoding_Scalp_Healthy.txt")
# cast OTU table
OTUTableBCH <- dcast(OTUTableBCH, Sample ~ OTU, value.var = "Count", fill = 0)
OTUTableBCH[1] <- NULL
N <- apply(OTUTableBCH,1,sum)
OTUTableBCH <- read_tsv("OTU_Bacterial_Baoding_Scalp_Healthy.txt") %>% filter(OTU %in% abundOTUsBCH)
OTUTableBCH <- dcast(OTUTableBCH, Sample ~ OTU, value.var = "Count", fill = 0)
OTUTableBCH[1] <- NULL
abundBCH <- OTUTableBCH/N
abundBCH <- apply(abundBCH,2,mean)
abundanceBCH <- abundBCH

#Store taxonomic information of prevalent OTUs (bacteria)
OTUTableBCH <- read_tsv("OTU_Bacterial_Baoding_Scalp_Healthy.txt") %>%
filter(OTU %in% abundOTUsBCH) %>%
select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>%
unique()
GenusBCH <- OTUTableBCH$Genus
ClassBCH <- OTUTableBCH$Class
PhylumBCH <- OTUTableBCH$Phylum
KingdomBCH <- OTUTableBCH$Kingdom

#Prepare phyloseq OTU table with prevalent OTUs (bacteria)
OTUTableBCH <- "OTU_Bacterial_Baoding_Scalp_Healthy.txt" %>%
read_tsv %>%
filter(OTU %in% abundOTUsBCH) %>%
select(OTU, Sample, Count) %>%
mutate(OTU = word(OTU, 1, sep = ";")) %>%
# Select only most abundant OTUs
# Delete the marked lines below to include all OTUs
group_by(OTU) %>%                   # <- Delete
mutate(TotalCount = sum(Count)) %>% # <- Delete
ungroup %>%                         # <- Delete
filter(TotalCount >= 100) %>%      # <- Delete
select(-TotalCount) %>%             # <- Delete
spread(OTU, Count, fill = 0) %>%
as.data.frame %>%
remove_rownames %>%
column_to_rownames("Sample") %>%
as.matrix %>%
otu_table(taxa_are_rows = FALSE)

#Prepare phyloseq taxa table (bacteria)
TaxaTableBCH <- "OTU_Bacterial_Baoding_Scalp_Healthy.txt" %>%
read_tsv %>%
filter(OTU %in% abundOTUsBCH) %>%
select(OTU, Domain = Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
mutate(OTU = word(OTU, 1, sep = ";")) %>%
unique %>%
as.data.frame %>%
remove_rownames %>%
column_to_rownames("OTU") %>%
as.matrix %>%
tax_table

PhyseqBCH <- phyloseq(OTUTableBCH, TaxaTableBCH)


#Step 1 of SPIEC-EASI: perform normalization, pseudo-count, and CLR-transformation of count data (just bacteria data)
data_clr_bch <- t(clr(OTUTableBCH+1, 1))



##Optional: subselect first 100 columns (OTUs) to test script
#data_clr_bac_small <- data_clr_bac[,1:100]


################################################
################################################
################################################
################################################
################################################


#Filter OTUs from table where prevalence is below 25% of dataset (fungi)
OTUTableFCH <- read_tsv("OTU_Fungal_Baoding_Scalp_Healthy.txt")
# cast OTU table
OTUTableFCH <- dcast(OTUTableFCH, Sample ~ OTU, value.var = "Count", fill = 0)
OTUTableFCH[1] <- NULL
# select most prevalent OTUs by occurrence frequency
OTUTableFCH <- 1*(OTUTableFCH>0)
freqFCH <- apply(OTUTableFCH, 2, mean)
freqFCH <- freqFCH[freqFCH >=0.25]
freqFCH <- as.matrix(freqFCH)
abundOTUsFCH <- rownames(freqFCH)
prevalenceFCH <- freqFCH[,1]


#calculate the abundance of prevalent OTUs (fungi)
OTUTableFCH <- read_tsv("OTU_Fungal_Baoding_Scalp_Healthy.txt")
# cast OTU table
OTUTableFCH <- dcast(OTUTableFCH, Sample ~ OTU, value.var = "Count", fill = 0)
OTUTableFCH[1] <- NULL
N <- apply(OTUTableFCH,1,sum)
OTUTableFCH <- read_tsv("OTU_Fungal_Baoding_Scalp_Healthy.txt") %>% filter(OTU %in% abundOTUsFCH)
OTUTableFCH <- dcast(OTUTableFCH, Sample ~ OTU, value.var = "Count", fill = 0)
OTUTableFCH[1] <- NULL
abundFCH <- OTUTableFCH/N
abundFCH <- apply(abundFCH,2,mean)
abundanceFCH <- abundFCH


#Store taxonomic information (fungi)
OTUTableFCH <- read_tsv("OTU_Fungal_Baoding_Scalp_Healthy.txt") %>%
filter(OTU %in% abundOTUsFCH) %>%
select(OTU, Kingdom, Phylum, Class, Order, Family, Genus) %>%
unique()
GenusFCH <- OTUTableFCH$Genus
ClassFCH <- OTUTableFCH$Class
PhylumFCH <- OTUTableFCH$Phylum
KingdomFCH <- OTUTableFCH$Kingdom

#Prepare phyloseq OTU table (fungi)
OTUTableFCH <- "OTU_Fungal_Baoding_Scalp_Healthy.txt" %>%
read_tsv %>%
filter(OTU %in% abundOTUsFCH) %>%
select(OTU, Sample, Count) %>%
mutate(OTU = word(OTU, 1, sep = ";")) %>%
# Select only most abundant OTUs
# Delete the marked lines below to include all OTUs
group_by(OTU) %>%                   # <- Delete
mutate(TotalCount = sum(Count)) %>% # <- Delete
ungroup %>%                         # <- Delete
filter(TotalCount >= 100) %>%      # <- Delete
select(-TotalCount) %>%             # <- Delete
spread(OTU, Count, fill = 0) %>%
as.data.frame %>%
remove_rownames %>%
column_to_rownames("Sample") %>%
as.matrix %>%
otu_table(taxa_are_rows = FALSE)

#Prepare phyloseq taxa table (fungi)
TaxaTableFCH <- "OTU_Fungal_Baoding_Scalp_Healthy.txt" %>%
read_tsv %>%
filter(OTU %in% abundOTUsFCH) %>%
select(OTU, Domain = Kingdom, Phylum, Class, Order, Family, Genus, Species) %>%
mutate(OTU = word(OTU, 1, sep = ";")) %>%
unique %>%
as.data.frame %>%
remove_rownames %>%
column_to_rownames("OTU") %>%
as.matrix %>%
tax_table

PhyseqFCH <- phyloseq(OTUTableFCH, TaxaTableFCH)

#Step 1 of SPIEC-EASI: perform normalization, pseudo-count, and CLR-transformation of count data (fungal)
data_clr_fch <- t(clr(OTUTableFCH+1, 1))


#Combine bacterial and fungal CLR-transformed data. This combined file will be used for steps 2-3 of SPIEC-EASI pipeline
data_clr_combined_bch <- cbind(data_clr_bch,data_clr_fch)

#Combine bacterial and fungal abundance and prevalence data by combining tables from separate domains. Want to retain abundance data as percentage of community within each domain, not total relative abundance of the two domains combined
freqBCH <- c(freqBCH,freqFCH)
abundBCH <- c(abundBCH,abundFCH)
prevalenceBCH <- freqBCH
abundanceBCH <- abundBCH
GenusCombinedBCH <- c(GenusBCH,GenusFCH)
ClassCombinedBCH <- c(ClassBCH,ClassFCH)
PhylumCombinedBCH <- c(PhylumBCH,PhylumFCH)
KingdomCombinedBCH <- c(KingdomBCH, KingdomFCH)

#Step 2 of SPIEC-EASI: Inverse covariance estimation (from SPIEC-EASI github)
sparseiCov <- function(data, method, npn=FALSE, verbose=FALSE, cov.output = TRUE, ...) {
    
    if (npn) data <- huge::huge.npn(data, verbose=verbose)
    
    args <- list(...)
    
    method <- switch(method, glasso = "glasso", mb = "mb", stop("Method not supported"))
    
    if (is.null(args$lambda.min.ratio)) args$lambda.min.ratio <- 1e-3
    
    if (method %in% c("glasso")) {
        do.call(huge::huge, c(args, list(x=data, method=method, verbose=verbose,
        cov.output = cov.output)))
        
    } else if (method %in% c('mb')) {
        est <- do.call(huge::huge.mb, c(args, list(x=data, verbose=verbose)))
        est$method <- 'mb'
        est$data <- data
        est$sym  <- ifelse(!is.null(args$sym), args$sym, 'or')
        return(est)
    }
}

sparseiCov_bch <- sparseiCov(data_clr_combined_bch, method="mb", lambda.min.ratio=1e-3, verbose=TRUE, nlambda=20)


#Step 3 of SPIEC-EASI: Model selection to pick right lambda (from SPIEC-EASI github)

#' Model selection for picking the right \code{lambda} penalty.
#' This is identical to huge::huge.stars except that the subsampling loop is replaced with an mclapply function to add parallelization capabilities.
#'
#' @param est an estimate/model as produced by the sparseiCov function
#' @param criterion character string specifying criterion/method for model selection accepts 'stars' [default], 'ric', 'ebic'
#' @param stars.thresh variability threshold for stars selection
#' @param ebic.gamma tuning parameter for ebic
#' @param stars.subsample.ratio The default value 'is 10*sqrt(n)/n' when 'n>144' and '0.8' when 'n<=144', where 'n' is the sample size.
#' @param rep.num number of subsamplings when \code{criterion} = stars.
#' @param ncores number of cores to use. Need multiple processers if \code{ncores > 1}
#' @param normfun normalize internally if data should be renormalized
#' @importFrom parallel mclapply
#' @export
icov.select <- function(est, criterion = 'stars', stars.thresh = 0.05, ebic.gamma = 0.5,
stars.subsample.ratio = NULL, rep.num = 20, ncores=1, normfun=function(x) x, verbose=FALSE) {
    gcinfo(FALSE)
    if (est$cov.input) {
        message("Model selection is not available when using the covariance matrix as input.")
        class(est) = "select"
        return(est)
    }
    if (!est$cov.input) {
        if (est$method == "mb" && is.null(criterion))
        criterion = "stars"
        if (est$method == "ct" && is.null(criterion))
        criterion = "ebic"
        n = nrow(est$data)
        d = ncol(est$data)
        nlambda = length(est$lambda)
        if (criterion == "ric") {
            if (verbose) {
                message("Conducting rotation information criterion (ric) selection....")
                #        flush.console()
            }
            if (n > rep.num) {
                nr = rep.num
                r = sample(n, rep.num)
            }
            if (n <= rep.num) {
                nr = n
                r = 1:n
            }
            out = .C("RIC", X = as.double(est$data), dd = as.integer(d),
            nn = as.integer(n), r = as.integer(r), nr = as.integer(nr),
            lambda_opt = as.double(0), PACKAGE = "huge")
            est$opt.lambda = out$lambda_opt/n
            rm(out)
            gc()
            if (verbose) {
                message("done\n")
                #        flush.console()
            }
            if (verbose) {
                message("Computing the optimal graph....")
                #        flush.console()
            }
            if (est$opt.lambda > max(cor(est$data)))
            est$refit = Matrix(0, d, d)
            else {
                if (est$method == "mb")
                est$refit = huge::huge.mb(est$data, lambda = est$opt.lambda,
                sym = est$sym, idx.mat = est$idx.mat, verbose = FALSE)$path[[1]]
                if (est$method == "glasso") {
                    if (!is.null(est$cov)) {
                        tmp = huge::huge.glasso(est$data, lambda = est$opt.lambda,
                        scr = est$scr, cov.output = TRUE, verbose = FALSE)
                        est$opt.cov = tmp$cov[[1]]
                    }
                    if (is.null(est$cov))
                    tmp = huge::huge.glasso(est$data, lambda = est$opt.lambda,
                    verbose = FALSE)
                    est$refit = tmp$path[[1]]
                    est$opt.icov = tmp$icov[[1]]
                    rm(tmp)
                    gc()
                }
                if (est$method == "ct")
                est$refit = huge::huge.ct(est$data, lambda = est$opt.lambda,
                verbose = FALSE)$path[[1]]
            }
            est$opt.sparsity = sum(est$refit)/d/(d - 1)
            if (verbose) {
                cat("done\n")
                #        flush.console()
            }
        }
        if (criterion == "ebic" && est$method == "glasso") {
            if (verbose) {
                cat("Conducting extended Bayesian information criterion (ebic) selection....")
                #        flush.console()
            }
            est$ebic.score = -n * est$loglik + log(n) * est$df + 4 * ebic.gamma * log(d) * est$df
            est$opt.index = which.min(est$ebic.score)
            est$refit = est$path[[est$opt.index]]
            est$opt.icov = est$icov[[est$opt.index]]
            if (est$cov.output)
            est$opt.cov = est$cov[[est$opt.index]]
            est$opt.lambda = est$lambda[est$opt.index]
            est$opt.sparsity = est$sparsity[est$opt.index]
            if (verbose) {
                message("done\n")
                #        flush.console()
            }
        }
        if (criterion == "stars") {
            if (is.null(stars.subsample.ratio)) {
                if (n > 144)
                stars.subsample.ratio = 10 * sqrt(n)/n
                if (n <= 144)
                stars.subsample.ratio = 0.8
            }
            
            #            for (i in 1:nlambda) merge[[i]] <- Matrix(0, d, d)
            
            if (verbose) {
                mes = "Conducting Subsampling....."
                message(mes, appendLF = FALSE)
                #        cat("\n")
                #        flush.console()
            }
            #    for (i in 1:rep.num) {
            premerge <- parallel::mclapply(1:rep.num, function(i) {
                #                if (verbose) {
                #                  mes <- paste(c("Conducting Subsampling....in progress:",
                #                    floor(100 * i/rep.num), "%"), collapse = "")
                #                  cat(mes, "\r")
                #                  flush.console()
                #                }
                #                merge <- replicate(nlambda, Matrix(0, d,d))
                ind.sample = sample(c(1:n), floor(n * stars.subsample.ratio),
                replace = FALSE)
                if (est$method == "mb")
                tmp = huge::huge.mb(normfun(est$data[ind.sample, ]), lambda = est$lambda,
                scr = est$scr, idx.mat = est$idx.mat, sym = est$sym,
                verbose = FALSE)$path
                if (est$method == "ct")
                tmp = huge::huge.ct(normfun(est$data[ind.sample, ]), lambda = est$lambda,
                verbose = FALSE)$path
                if (est$method == "glasso")
                tmp = huge::huge.glasso(normfun(est$data[ind.sample, ]), lambda = est$lambda,
                scr = est$scr, verbose = FALSE)$path
                #                for (j in 1:nlambda) merge[[j]] <- merge[[j]] + tmp[[j]]
                
                rm(ind.sample)
                gc()
                return(tmp)
            }, mc.cores=ncores)
            #  }
            # merge <- lapply(merge, as.matrix)
            #      merge <- lapply(merge, simplify2array)
            #      est$merge <- lapply(1:dim(merge)[3], function(i) merge[,,i]/rep.num)
            
            merge <- Reduce(function(l1, l2) lapply(1:length(l1),
            function(i) l1[[i]] + l2[[i]]), premerge, accumulate=FALSE)
            
            if (verbose) {
                message("done")
                #        cat("\n")
                #        flush.console()
            }
            est$variability = rep(0, nlambda)
            est$merge <- vector('list', nlambda)
            for (i in 1:nlambda) {
                est$merge[[i]]  <- merge[[i]]/rep.num
                est$variability[i] <- 4 * sum(est$merge[[i]] * (1 - est$merge[[i]]))/(d * (d - 1))
            }
            est$opt.index = max(which.max(est$variability >=
            stars.thresh)[1] - 1, 1)
            est$refit = est$path[[est$opt.index]]
            est$opt.lambda = est$lambda[est$opt.index]
            est$opt.sparsity = est$sparsity[est$opt.index]
            if (est$method == "glasso") {
                est$opt.icov = est$icov[[est$opt.index]]
                if (!is.null(est$cov))
                est$opt.cov = est$cov[[est$opt.index]]
            }
        }
        est$criterion = criterion
        class(est) = "select"
        return(est)
    }
}

spiec.easi.combined_bch <- icov.select(sparseiCov_bch, rep.num = 20, ncores = 16, verbose=TRUE)

#Merge bacterial and fungal phyloseq
PhyseqCombined_BC_Healthy <- merge_phyloseq(PhyseqBCH,PhyseqFCH)

# 6 Convert to igraph and network objects
ig.combined_bch <- spiec.easi.combined_bch %>%
.$refit %>%
adj2igraph(vertex.attr = list(name = taxa_names(PhyseqCombined_BC_Healthy)))

vertex_attr(ig.combined_bch)
edge_attr(ig.combined_bch)

betweenness_bch <- betweenness(ig.combined_bch)
ig.combined_bch <- set_vertex_attr(ig.combined_bch,"betweenness",value=betweenness_bch)
degree_bch <- degree(ig.combined_bch)
ig.combined_bch <- set_vertex_attr(ig.combined_bch,"degree", value=degree_bch)


# set OTUs vertex attributes
ig.combined_bch <- set_vertex_attr(ig.combined_bch, "prevalence", value=prevalenceBCH)
ig.combined_bch <- set_vertex_attr(ig.combined_bch, "abundance", value=abundanceBCH)
GenusCombinedBCH[is.na(GenusCombinedBCH)] <- "unclassified_Genus"
ig.combined_bch <- set_vertex_attr(ig.combined_bch, "Genus", value=GenusCombinedBCH)
ClassCombinedBCH[is.na(ClassCombinedBCH)] <- "unclassified_Class"
ig.combined_bch <- set_vertex_attr(ig.combined_bch, "Class", value=ClassCombinedBCH)
ig.combined_bch <- set_vertex_attr(ig.combined_bch, "Phylum", value=PhylumCombinedBCH)
ig.combined_bch <- set_vertex_attr(ig.combined_bch, "Kingdom", value=KingdomCombinedBCH)


#Extrac OTU names for later use
OTUNamescombined_BC_Healthy <- PhyseqCombined_BC_Healthy %>% otu_table %>% colnames

# 7 Determine correlation weights and direction (positive or negative
# correlation), and save
elist.weighted.combined_bch <- spiec.easi.combined_bch %>%
getOptBeta %>%
symBeta(mode = 'maxabs') %>%
summary %>%
# Assuming OTU numbers in i and j correlate to the order of OTUs in
# se.summer3zw, restore OTU names
left_join(data.frame(OTU1 = OTUNamescombined_BC_Healthy, i = 1:length(OTUNamescombined_BC_Healthy))) %>%
left_join(data.frame(OTU2 = OTUNamescombined_BC_Healthy, j = 1:length(OTUNamescombined_BC_Healthy)))

elist.weighted.combined_bch <- elist.weighted.combined_bch %>%
mutate(CorSign = ifelse(x > 0, "Positive", "Negative")) %>%
mutate(CorMagnitude = abs(x))

#Save elist document
write.tidy(elist.weighted.combined_bch,"elist.weighted.combined_bch.txt")

ig.combined_bch <- set_edge_attr(ig.combined_bch, "CorrelationSign", value=elist.weighted.combined_bch$CorSign)
ig.combined_bch <- set_edge_attr(ig.combined_bch, "CorrelationMagnitude", value=elist.weighted.combined_bch$CorMagnitude)


# remove weak edges (< 0.15)
Weak = which(E(ig.combined_bch)$CorrelationMagnitude<0.15)
ig.combined_bch = delete_edges(ig.combined_bch, Weak) 

# remove vertices without strong edges (independent vertices)
Isolated = which(degree(ig.combined_bch)==0)
ig.combined_bch = delete_vertices(ig.combined_bch, Isolated)

# plot the network
library(RColorBrewer)
Cor <- E(ig.combined_bch)$CorrelationMagnitude


par(mfrow=c(1,2))
coords <- layout.fruchterman.reingold(ig.combined_bch)
E(ig.combined_bch)$color[E(ig.combined_bch)$CorrelationSign == "Negative"] <- 'red'
E(ig.combined_bch)$color[E(ig.combined_bch)$CorrelationSign == "Positive"] <- 'blue'


pal <- brewer.pal(length(unique(ig.combined_bca$Kingdom)), "Set3")

pal2 <- c("#ff7f00", "#c41185",  "#320956", "#b15928","black","purple","blue", "red","green","brown","orange")

fine = 2000 # this will adjust the resolving power.
plot(ig.combined_bch,
     vertex.label=NA,
     vertex.size= sqrt(V(ig.combined_bch)$abundance*2000)+2,
     vertex.color= pal2[as.numeric(as.factor(vertex_attr(ig.combined_bch, "Kingdom")))],
     edge.width=Cor*4,
     layout=coords)
legend("topleft", legend=levels(as.factor(V(ig.combined_bch)$Kingdom)), bg= "white",col = pal2, bty="o",pch=20, pt.cex = 1.3, cex = 0.6, text.col=pal2 , horiz = FALSE, inset = c(0.3, 0.83))
title("Baoding Scalp - Healthy", line=-1)

# Perform fast-greedy community detection on network graph
kc = fastgreedy.community(ig.combined_bch)

# Perform edge-betweenness community detection on network graph
gc = edge.betweenness.community(ig.combined_bch)

# Plot community networks determined by fast-greedy and edge-betweenness methods side-by-side
# color the vertex by class rank
#pal <- brewer.pal(length(unique(V(ig.summer)$Class)), "Paired")
pal <- c("#1f78b4","#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00", "#cab2d6",  "#ffff99", "#b15928","grey", "black")
my_color =pal[as.numeric(as.factor(V(ig.combined_bch)$ClassCombinedBCH))]
pal2 <- c("red", "grey")
my_color2 =pal2[as.numeric(as.factor(E(ig.combined_bch)$CorrelationSignBCH))]


par(mfrow = c(1, 2))
plot(kc, ig.combined_bch,
vertex.label = NA,
vertex.size = 5,
mark.groups=kc,
edge.color=my_color2,
col=my_color,
#membership=NULL,
#vertex.color= my_color,
main = "Fast-greedy \nmodularity=0.6 groups=8",
edge.width=Cor*5,
layout=coords)
plot(gc, ig.combined_bch,
vertex.label =NA,
vertex.size = 5,
edge.color=my_color2,
col=my_color,
#vertex.color= pal[as.numeric(as.factor(vertex_attr(ig.summer, "Class")))],
main = "Edge-betweenness \nmodularity=0.59 groups=9",
edge.width=Cor*5,
#layout_nicely(ig.bac.small),
layout=coords)




##############draw degree distribution

dd.BCH <- degree_distribution(ig.combined_bch, cumulative=FALSE)
sum(seq_along(dd.BCH)*dd.BCH)-1

dd.BCA <- degree_distribution(ig.combined_bca, cumulative=FALSE)
sum(seq_along(dd.BCA)*dd.BCA)-1

dd.DCH <- degree_distribution(ig.combined_dch, cumulative=FALSE)
sum(seq_along(dd.DCH)*dd.DCH)-1

dd.DCA <- degree_distribution(ig.combined_dca, cumulative=FALSE)
sum(seq_along(dd.DCA)*dd.DCA)-1

png("degree_distribution_cheek.png",width=4,height=4,units="in",res=1200)
par(mar=c(4,4,4,4))
plot(seq_along(dd.BCH)-1, dd.BCH, type='b', xlim=c(0,30), ylim=c(0,0.2),
ylab="Frequency", xlab="Degree", col='blue')
points(seq_along(dd.BCA)-1, dd.BCA, type='b', col='green')
points(seq_along(dd.DCH)-1, dd.DCH, type='b', col='red')
points(seq_along(dd.DCA)-1, dd.DCA, type='b', col='black')
legend("topright", c("Baoding Healthy","Baoding Acne","Dalian Healthy", "Dalian Acne"), col=c("blue","green","red","black"),pch=1, lty=1, cex=0.6)
title(main="Cheek")
dev.off()
      
#Network stability estimation by natural connectivity
natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
nc.attack <- function(ig) {
  hubord <- order(rank(betweenness(ig)), rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.95), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  })
}


natcon <- function(ig) {
  N <- vcount(ig)
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
  nc / (N - log(N))
}
nc.attack <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.95), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  })
}

nc.BCH <- nc.attack(ig.combined_bch)
nc.DCH <- nc.attack(ig.combined_dch)
nc.DCA <- nc.attack(ig.combined_dca)

##Plot Natural connectivity

png("connectivity_cheek.png",width=4,height=4,units="in",res=1200)
plot(seq(0,0.8,len=length(nc.Winter)),nc.Winter, type='l', col='blue', ylim=c(0,max(nc.Winter)), xlab="Proportion of Removed Nodes", ylab="Natural Connectivity",font.lab=2)
points(seq(0,0.8,len=length(nc.Spring)),nc.Spring,type='l', col='green')
points(seq(0,0.8,len=length(nc.Summer)),nc.Summer,type='l', col='red')
points(seq(0,0.8,len=length(nc.DCA)),nc.DCA,type='l')
dev.off()
