library(diptest)
library(mixtools)
library(sROC)



EL_statistic<-function(expression_pattern,ret.pvals=FALSE){
    ep<-expression_pattern
    
    CDF<-kCDF(ep,xgrid=ep,adjust=1)
    nep<-CDF$Fhat[match(ep,CDF$x)]
    
    el<-log(nep/(1-nep))
    
    if (ret.pvals){
        pvals<-rep(NA,length(nep))
        pvals[el>=0]<-(1-nep[el>=0])
        pvals[el<0]<-nep[el<0]
        
        return(list(el=el,pvals=pvals))
    }
    
    return(el)
}


detect_bimodal_signal<-function(signal,display=TRUE){
        dip_pvals<-dip.test(signal)
        positions<-names(dip_pvals)
        
        if (display){
                hist(signal,breaks=round(length(signal)/10))
        }
        
        
        if (dip_pvals$p.value<0.05){
            
            mixmdl = normalmixEM(signal,maxit = 5000)
        
            mixmdl$posterior<-mixmdl$posterior[,order(mixmdl$mu)]
            mixmdl$sigma<-mixmdl$sigma[order(mixmdl$mu)]
            mixmdl$lambda<-mixmdl$lambda[order(mixmdl$mu)]
            mixmdl$mu<-mixmdl$mu[order(mixmdl$mu)]
            
            if (display){
                
            tmp<-hist(signal,50)
        
            par(new=TRUE)
            plot(sort(signal),mixmdl$lambda[1]*dnorm(sort(signal),mean = mixmdl$mu[1],sd = mixmdl$sigma[1]),type = 'l',col='darkgreen',lwd=3,
                 ylim=c(0,max(tmp$density)),xlim=c(tmp$breaks[1],tmp$breaks[length(tmp$breaks)]),
                 frame.plot = FALSE,
                 xaxt='n',yaxt='n',xlab='',ylab='')
        
            par(new=TRUE)
            plot(sort(signal),mixmdl$lambda[2]*dnorm(sort(signal),mean = mixmdl$mu[2],sd = mixmdl$sigma[2]),type = 'l',col='darkred',lwd=3,
                 xlim=c(tmp$breaks[1],tmp$breaks[length(tmp$breaks)]),ylim=c(0,max(tmp$density)),
                 frame.plot = FALSE,
                xaxt='n',yaxt='n',xlab='',ylab='')
        
            
            }
            
            tmp<-sort(mixmdl$x)[which(log10(mixmdl$posterior[order(mixmdl$x),2]/mixmdl$posterior[order(mixmdl$x),1])>1)]
            discTH<-min(tmp[which(tmp>min(mixmdl$mu) & tmp<max(mixmdl$mu))])
            
        }else{
                
            discTH<-Inf
            
            }
        return(discTH)
        }
