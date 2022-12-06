#!/usr/bin/env Rscript

# Retrieve command line arguments
args<-commandArgs(trailingOnly=TRUE)
file.const<-args[1]
file.prof<-args[2]
file.sol<-args[3]
plot.name<-args[4]

# Set up functions
prof<-function(plot.name){
    
    # Read in constants
    const<-read.table(file.const,header=TRUE)

    # Set up vertical profile
    z<-c()
    for (k in seq(1,const$nk)){
        z<-c(z,(const$nk-k+1)*const$dz*exp(const$zexp*(const$nk-k)))
    }
    z<-rev(z)
    z.fine<-seq(min(z),max(z),length.out=1000)

    # Interpolate profile
    if (suppressWarnings(library(pracma,quietly=TRUE,logical=TRUE))){
        df.pchip<-pchip(z,df(z,const),z.fine)
        f.pchip<-c()
        idxloop<-seq(2,(length(z.fine)))
        if (const$dir < 0){idxloop<-rev(idxloop)}
        for (k in idxloop){
            if (const$dir*z.fine[k] < const$dir*const$zdep){
                f.pchip<-c(f.pchip,NA)
            } else {
                f.low<-f.pchip[length(f.pchip)]
                if (is.null(f.low) || is.na(f.low)){f.low<-0}
                delF<-const$dir*(z.fine[k]-z.fine[k-1])*(df.pchip[k]+df.pchip[k-1])/2
                f.pchip<-c(f.pchip,f.low+delF)
            }
        }
        f.pchip<-c(f.pchip,NA)
        if (const$dir < 0){f.pchip<-rev(f.pchip)}
    } else {
        df.pchip<-rep(NA,length(z.fine))
        f.pchip<-rep(NA,length(z.fine))
    }

    # Read in model profile
    f.model<-read.table(file.prof)
    sol.model<-scan(file.sol,quiet=TRUE)

    # Plot profiles
    postscript(file=paste(plot.name,".ps",sep=''),horizontal=FALSE,
               paper="special",width=7,height=7,bg='white')
    layout(matrix(c(1,2,3,3),ncol=2,byrow=TRUE),heights=c(0.8,0.2))
    par.default<-par(mar=c(2,3,3,2))
    plot(df(z.fine,const),z.fine,type='l',lwd=3,panel.first={
        abline(v=0,col='grey');
        abline(h=z,col='grey');
        lines(df.pchip,z.fine,type='l',lwd=3,col='orange')
    },
    main='Function',
    ylab='Height',
    yaxt='n')
    axis(2,at=z,labels=round(z),las=1)
    points(df(z,const),z,lwd=3,pch=19)
    plot(f(z.fine,const),z.fine,type='l',lwd=3,panel.first={
        abline(v=0,col='grey');
        abline(h=z,col='grey');
        abline(h=const$zdep,col='blue',lwd=2);
        abline(v=const$a,col='red',lwd=2);
        lines(f.pchip,z.fine,type='l',lwd=3,col='orange')
        lines(f.model[,2],f.model[,1],type='l',lwd=3,col='green')
        abline(h=sol.model,lwd=2,col='green',lty=2)
    },
    main='Integral',
    ylab='Height',
    yaxt='n')
    axis(2,at=z,labels=round(z),las=1)
    points(f(z,const),z,lwd=3,pch=19)
    par(mar=c(0,0,0,0))
    plot(1,type='n',axes=FALSE,xlab='',ylab='')
    legend('center',legend=c('Analytic','PCHIP','Departure Level','RHS','Model','Model Solution'),
           col=c('black','orange','blue','red','green','green'),ncol=2,
           lwd=c(3,3,2,2,3,2),lty=c(1,1,1,1,1,2),cex=1.2)
    par(par.default)
    dev.off()
}
df<-function(z,const){
    return(const$alpha*sin(2*pi/const$L*z) + const$beta*z)
}
f.indef<-function(z,const){
    return(-const$alpha*const$L/(2.*pi)*cos(2*pi/const$L*z) + (const$beta/2)*z^2)
}
f<-function(z,const){
    return(f.indef(z,const)-f.indef(const$zdep,const))
}

# Generate plots
prof(plot.name)
