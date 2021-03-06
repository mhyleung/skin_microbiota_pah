##The following script will perform the Sloan's model of random assembly, including generating Akaike information criterion
##scores, and prevalence and abundance data for Sloan modelling figure as presented on manuscript.
library(devtools)
library(readr)
library(Hmisc)
library(bbmle)
library(wilkoxmisc)
library(reshape2)
library(minpack.lm)

spp <- read_tsv("OTU_table_clean5_wt_control_400_490_65.tidy_fixed_cast_10452_Dalian_Cheek.txt")

#If stats=FALSE the function will return a table of observed and predicted values for each otu.
sncm.fit_f <- function(spp, pool=NULL, stats=FALSE, taxon=NULL){
    require(minpack.lm)
    require(Hmisc)
    require(stats4)
    
    options(warn=-1)
    
    #Calculate the number of individuals per community
    N <- mean(apply(spp, 1, sum))
    
    
    #Calculate the average relative abundance of each taxa across communities
    if(is.null(pool)){
        p.m <- apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    } else {
        p.m <- apply(pool, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    }
    
    #Calculate the occurrence frequency of each taxa across communities
    spp.bi <- 1*(spp>0)
    freq <- apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]
    
    #Combine
    C <- merge(p, freq, by=0)
    C <- C[order(C[,2]),]
    C <- as.data.frame(C)
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
    p <- C.0[,2]
    freq <- C.0[,3]
    names(p) <- C.0[,1]
    names(freq) <- C.0[,1]
    
    #Calculate the limit of detection
    d = 1/N
    
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
    m.ci <- confint(m.fit, 'm', level=0.95)
    
    ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
    sncm.LL <- function(m, sigma){
        R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
        R = dnorm(R, 0, sigma)
        -sum(log(R))
    }
    m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
    
    ##Calculate Akaike's Information Criterion (AIC)
    aic.fit <- AIC(m.mle, k=2)
    bic.fit <- BIC(m.mle)
    
    ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
    
    pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for binomial model
    bino.LL <- function(mu, sigma){
        R = freq - pbinom(d, N, p, lower.tail=FALSE)
        R = dnorm(R, mu, sigma)
        -sum(log(R))
    }
    bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    
    aic.bino <- AIC(bino.mle, k=2)
    bic.bino <- BIC(bino.mle)
    
    ##Goodness of fit for binomial model
    bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
    Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
    
    bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma){
        R = freq - ppois(d, N*p, lower.tail=FALSE)
        R = dnorm(R, mu, sigma)
        -sum(log(R))
    }
    pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    
    aic.pois <- AIC(pois.mle, k=2)
    bic.pois <- BIC(pois.mle)
    
    ##Goodness of fit for Poisson model
    pois.pred <- ppois(d, N*p, lower.tail=FALSE)
    Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
    
    pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Results
    if(stats==TRUE){
        fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqr.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
        fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
        return(fitstats)
    } else {
        A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
        A <- as.data.frame(A)
        colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
        if(is.null(taxon)){
            B <- A[order(A[,1]),]
        } else {
            B <- merge(A, taxon, by=0, all=TRUE)
            row.names(B) <- B[,1]
            B <- B[,-1]
            B <- B[order(B[,1]),]
        }
        return(B)
    }
}

stats.otu <- sncm.fit_f(spp)
write.table(stats.otu, "prediction_dalian_cheek_10452.txt", sep="\t", col.names=NA)

#If stats=TRUE the function will return fitting statistics.
sncm.fit_t <- function(spp, pool=NULL, stats=TRUE, taxon=NULL){
    require(minpack.lm)
    require(Hmisc)
    require(stats4)
    
    options(warn=-1)
    
    #Calculate the number of individuals per community
    N <- mean(apply(spp, 1, sum))
    
    #Calculate the average relative abundance of each taxa across communities
    if(is.null(pool)){
        p.m <- apply(spp, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    } else {
        p.m <- apply(pool, 2, mean)
        p.m <- p.m[p.m != 0]
        p <- p.m/N
    }
    
    #Calculate the occurrence frequency of each taxa across communities
    spp.bi <- 1*(spp>0)
    freq <- apply(spp.bi, 2, mean)
    freq <- freq[freq != 0]
    
    #Combine
    C <- merge(p, freq, by=0)
    C <- C[order(C[,2]),]
    C <- as.data.frame(C)
    C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),] #Removes rows with any zero (absent in either source pool or local communities)
    p <- C.0[,2]
    freq <- C.0[,3]
    names(p) <- C.0[,1]
    names(freq) <- C.0[,1]
    
    #Calculate the limit of detection
    d = 1/N
    
    ##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
    m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE), start=list(m=0.1))
    m.ci <- confint(m.fit, 'm', level=0.95)
    
    ##Fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
    sncm.LL <- function(m, sigma){
        R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
        R = dnorm(R, 0, sigma)
        -sum(log(R))
    }
    m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
    
    ##Calculate Akaike's Information Criterion (AIC)
    aic.fit <- AIC(m.mle, k=2)
    bic.fit <- BIC(m.mle)
    
    ##Calculate goodness-of-fit (R-squared and Root Mean Squared Error)
    freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1-p), lower.tail=FALSE)
    Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
    
    pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for binomial model
    bino.LL <- function(mu, sigma){
        R = freq - pbinom(d, N, p, lower.tail=FALSE)
        R = dnorm(R, mu, sigma)
        -sum(log(R))
    }
    bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    
    aic.bino <- AIC(bino.mle, k=2)
    bic.bino <- BIC(bino.mle)
    
    ##Goodness of fit for binomial model
    bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
    Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
    
    bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Calculate AIC for Poisson model
    pois.LL <- function(mu, sigma){
        R = freq - ppois(d, N*p, lower.tail=FALSE)
        R = dnorm(R, mu, sigma)
        -sum(log(R))
    }
    pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
    
    aic.pois <- AIC(pois.mle, k=2)
    bic.pois <- BIC(pois.mle)
    
    ##Goodness of fit for Poisson model
    pois.pred <- ppois(d, N*p, lower.tail=FALSE)
    Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
    RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))
    
    pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
    
    ##Results
    if(stats==TRUE){
        fitstats <- data.frame(m=numeric(), m.ci=numeric(), m.mle=numeric(), maxLL=numeric(), binoLL=numeric(), poisLL=numeric(), Rsqr=numeric(), Rsqrf.bino=numeric(), Rsqr.pois=numeric(), RMSE=numeric(), RMSE.bino=numeric(), RMSE.pois=numeric(), AIC=numeric(), BIC=numeric(), AIC.bino=numeric(), BIC.bino=numeric(), AIC.pois=numeric(), BIC.pois=numeric(), N=numeric(), Samples=numeric(), Richness=numeric(), Detect=numeric())
        fitstats[1,] <- c(coef(m.fit), coef(m.fit)-m.ci[1], m.mle@coef['m'], m.mle@details$value, bino.mle@details$value, pois.mle@details$value, Rsqr, Rsqr.bino, Rsqr.pois, RMSE, RMSE.bino, RMSE.pois, aic.fit, bic.fit, aic.bino, bic.bino, aic.pois, bic.pois, N, nrow(spp), length(p), d)
        return(fitstats)
    } else {
        A <- cbind(p, freq, freq.pred, pred.ci[,2:3], bino.pred, bino.pred.ci[,2:3])
        A <- as.data.frame(A)
        colnames(A) <- c('p', 'freq', 'freq.pred', 'pred.lwr', 'pred.upr', 'bino.pred', 'bino.lwr', 'bino.upr')
        if(is.null(taxon)){
            B <- A[order(A[,1]),]
        } else {
            B <- merge(A, taxon, by=0, all=TRUE)
            row.names(B) <- B[,1]
            B <- B[,-1]
            B <- B[order(B[,1]),]
        }
        return(B)
    }
}

stats.all <- sncm.fit_t(spp)
write_tsv(stats.all, "stats_dalian_cheek_10452.txt")

