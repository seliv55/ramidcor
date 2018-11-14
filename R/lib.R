binom<-function(n){ #calculation of binomial cefficients
 bf<-numeric(n); bf[1]<-1;
 if(n>1){ for(i in 2:n) bf[i]<-bf[i-1]*i;  bc<-numeric(n); 
   for(i in 1:(n-1)) bc[i]=bf[n]/bf[i]/bf[n-i];}
   else  bc<-1;
     bc }

carbons<-function(n, bc) {#natural distribution of 13C according to the number of carbons in whole fragment
  pc12<-0.989; pc13<-0.011;
    mc<-numeric(n+1); mc[1]<-pc12^n;
 if(n>1) {for(i in 1:(n-1)) mc[i+1]<-bc[i]*pc12^(n-i)*pc13^i; mc[n+1]<-pc13^n;}
   else {if(n==1) mc[2]<-pc13;}
       mc }
    
silicium<-function(mc){#correction of natural mass distribution accounting for Si
  pSi0<-0.9223; pSi1<-0.0467; pSi2<-0.031;
   n<-length(mc)
    m<-numeric(n+2)
     m[1]<-mc[1]*pSi0;
     m[2]<-mc[2]*pSi0+mc[1]*pSi1;
  for(i in 3:n) m[i]<-mc[i]*pSi0+mc[i-1]*pSi1+mc[i-2]*pSi2;
     m[n+1]<-mc[n]*pSi1+mc[n-1]*pSi2;
       m[n+2]<-mc[n]*pSi2;
         m }
      
sulfur<-function(mc){#correction of natural mass distribution accounting for S
  pS0<-0.9504; pS1<-0.0075; pS2<-0.0421;
   n<-length(mc)
    m<-numeric(n+2)
     m[1]<-mc[1]*pS0;
     m[2]<-mc[2]*pS0+mc[1]*pS1;
  for(i in 3:n) m[i]<-mc[i]*pS0+mc[i-1]*pS1+mc[i-2]*pS2;
     m[n+1]<-mc[n]*pS1+mc[n-1]*pS2;
       m[n+2]<-mc[n]*pS2;
         m }
      
tdist<-function(nc,nsi,ns){ #final natural mass distribution
   bcc<-binom(nc);
    m<-carbons(nc, bcc);
 if(nsi>0) for(i in 1:nsi)  m<-silicium(m);
 if(ns>0) for(i in 1:ns)  m<-sulfur(m);
      m }
 
elim<-function(msd,f1){#normalization of GCMS data accounting for the loss of protons
  msd1<-matrix(nrow=nrow(msd),ncol=ncol(msd),0)
    for (i in 2:(ncol(msd)-1)) msd1[,i-1]<-as.numeric(msd[,i]*(1+f1)-msd[,i+1]*f1);
      msd1[,ncol(msd1)-1]<-as.numeric(msd[,ncol(msd)])
     return(msd1)  }
   
mtr<-function(nfrg,nmass,nC,nSi,nS){ #various numbers of labeled 13C with tails of natural distributions
    mt<-tdist(nC,nSi,nS); 
     lem<-length(mt); if(lem<nmass) for(i in (lem+1):nmass) mt[i]<-0.0;
      mt<-mt[1:nmass];
       mt<-mt/sum(mt);
        mm<-numeric(nmass*(nfrg+1));          
         attr(mm,"dim")<-c(nfrg+1,nmass);
          mm[1,]<-mt;
       for(i in 1:nfrg) { mt<-tdist(nC-i,nSi,nS); lem<-length(mt);
         if(lem<(nmass-i)) {for(k in (lem+1):(nmass-i)) mt[k]<-0.0;}
             mt<-mt[1:(nmass-i)];  mt<-mt/sum(mt);
               for(j in 1:length(mt)) mm[i+1,j+i]<-mt[j];}
              mm}
              
mdistr<-function(nreal,msd,mm,nln){ #label incorporation             
          fr<-msd; colon<-ncol(fr); for(i in 1:colon) fr[,i]<-0;
 for(j in 1:(nreal+1)) {fr[,j]<-msd[,j]/mm[j,j];
   for(i in j:colon) msd[,i]<-msd[,i]-fr[,j]*mm[j,i];}
# normalization:
#     for(i in 2:(colon-2)) fr[i]<- fr[i]/fr[colon];
      fr }


basln<-function(vec,pos=length(vec),ofs=0){# baseline
   bas<-mean(vec[vec<(2*min(vec))])
 return(bas*5)}

readcdf<-function(fi) {
 nc <- nc_open(fi, readunlim=FALSE)  #open cdf file
   rett<-ncvar_get( nc, "scan_acquisition_time" )
   tiv<-ncvar_get( nc, "total_intensity" )
   npoint<-ncvar_get( nc, "point_count" )
     mz<-ncvar_get( nc, "mass_values" )
     iv<-ncvar_get( nc, "intensity_values" )
   nc_close( nc )
        return(list(mz,iv,npoint,rett,tiv))    }
        
  savplt<-function(mm,mm0,nma,plname){
  png(paste("../graf/",plname,"png",sep=""))
  par(mfrow=c(2,1))
   plot(mm[,2],xlim=c(nma-50,nma+50))
   plot(mm0[,1],xlim=c(nma-50,nma+50))
   dev.off()
  }

info<-function(mz,iv,npoint){
#  mz,iv,npoint: mz, intensities and number mz points in every scan
      j<-1
     mzpt<-numeric() # number of m/z points in each pattern
     tpos<-numeric() # initial time position for each m/z pattern 
     mzi<-numeric()  # initial value for each m/z pattern presented in the CDF file
     mzind<-numeric()# index in mz array corresponding to mzi
     mzrang<-list()  # list of mz patterns presented in the .CDF
  mzpt[j]<-npoint[1]; tpos[j]<-1; mzi[j]<-mz[1]; imz<-1; mzind[j]<-imz
  mzrang[[1]]<-mz[1:mzpt[1]];
    for(i in 2:length(npoint)) { imz<-imz+npoint[i-1];# mz index
     if(mzi[j]!=mz[imz]){  j<-j+1; tpos[j]<-i;  mzpt[j]<-npoint[i]; mzi[j]<-mz[imz];
      mzind[j]<-imz; mzrang[[j]]<-mz[(mzind[j]):(mzind[j]-1+mzpt[j])] }
    }
  tpos[length(tpos)+1]<- length(npoint) # add the last timepoint
  return(list(mzpt,tpos,mzind,mzrang))
  }

ftitle<-function(){paste("Raw_Data_File", "cells", "tracer_molecule", "labelled_positions","abundance", "injection","Replicate", "Incubation_time", "Metabolite_name", "CHEBI","atomic_positions", "Empirical_formula", "retention(min)", "mz_monitored", "signal_intensity", "isotopologue", "isotologue_abundance")}

wphen<-function(fi,nm,fragg, formul, rtt, pikmz,delta,corr=delta){
tracer<-c("Gluc","12Glc","Glutamina","UGln","[cC][oO][lL][dD]")
inctime<-c(40,'0[Hh]','6[hH]',24,"-")
cells<-c("A549","NCI","BEAS2B","SW620","HUVEC")
replicate<-c("R[0-9]_","G[0-9]_","_[0-9]_","[ _][0-9][ _]")

for(i in 1:length(tracer)) if(grepl(tracer[i],fi)) break;
trac<-switch(i,c("D-[1,2-C13]-Glucose","1,1,0,0,0,0",50),
c("D-[1,2-C13]-Glucose","1,1,0,0,0,0",50),
c("[3-C13]-Glutamine","0,0,1,0,0",100),
c("[U-C13]-Glutamine","1,1,1,1,1",100),
c("Glucose","0,0,0,0,0,0",100)
)
  for(incub in inctime) if(grepl(incub,fi)) break;
  for(cel in cells) if(grepl(cel,fi)) break;
  for(rep in replicate) if(grepl(rep,fi)) break;
      l=regexpr(rep, fi)+1
           inj<-substr(fi,nchar(fi),nchar(fi))
           rep<-substr(fi,l,l)
      nikiso<-paste("m",c(-1:(length(pikmz)-2)),sep="")
  return(paste(shQuote(paste(fi,'.CDF',sep='')),cel,trac[1],trac[2],trac[3],rep,inj,strsplit(incub,'[',fixed=T)[[1]][1],nm,"chebi",fragg, formul, rtt, pikmz,delta,nikiso,corr))
   }
fitG <-function(x,y,mu,sig,scale){
# x,y: x and y values for fitting
# mu,sig,scale: initial values for patameters to fit  
  f = function(p){
    d = p[3]*dnorm(x,mean=p[1],sd=p[2])
    sum((d-y)^2)
  }
  optim(c(mu,sig,scale),f,method="CG")
 # nlm(f, c(mu,sig,scale))
# output: optimized parameters
   }
  
fitdist<-function(ymat,nma,pint=5,cini=2,fsig=1.5,fsc=2.){ # fits distributions
# x: vector of x-values
# ymat: matrix of experimental values where columns are time courses for sequential mz
# nma: point of maximal value
# pint: half interval taken for fitting
# cini: initial column number
  cfin<-ncol(ymat)#cini+nmi-1;
  nmi<-cfin-cini+1 #ncol(ymat)-1;
   fscale<-numeric()
   xe<-c((nma-pint):(nma+pint));    facin<-max(ymat[nma,]);
   yemat<-ymat[(nma-pint):(nma+pint),cini:cfin]/facin
      yfmat<-yemat
          mu<-xe[pint+1]
          sig<-(xe[2*pint]-xe[2])/fsig
   for(i in 1:nmi){
          scale<-yemat[pint,i]*sig/fsc
   fp<-fitG(xe,yemat[,i],mu,sig,scale)
    fscale[i]<-fp$par[3]*facin
    yfmat[,i]<-fp$par[3]*dnorm(xe,mean=fp$par[1],sd=fp$par[2])
#    fscale[i]<-fp$estimate[3]*facin
#    yfmat[,i]<-fp$estimate[3]*dnorm(xe,mean=fp$estimate[1],sd=fp$estimate[2])
#   mu<-fp$par[1];  sig<-fp$par[2];# scale<-fp$par[3]
   }
   list(xe,yemat,yfmat,fscale)
#   xe: x-values used for fit
#   yemat: matrix of experimental intensities
#   yfmat: matrix of fitted intensities
#   fscale: areas of peaks
}
     
plal<-function(fi,x,me,mf){# plots intensities from matrix mm; nma - position of peaks; abs - 0 or 1 depending on mm
# fi: file to plot in
# x: vector of x-values
# me: matrix of experimental values where columns are time courses for sequential mz
# mf: matrix of fittings corresponding to me
    fi<-strsplit(fi,"CDF")[[1]][1]
  png(paste("../graf/",fi,"png",sep=""))
  x_range<-range(x[1],x[length(x)])
  g_range <- range(0,1)
  nkriv<-ncol(me); sleg<-"m0"
  plot(x,me[,1], xlim=x_range, ylim=g_range,col=1)
  lines(x,mf[,1],col=1, lty=1)
   for(i in 2:nkriv){ sleg<-c(sleg,paste("m",i-1))
    points(x,me[,i],pch=i,col=i)
    lines(x,mf[,i],col=i, lty=i)
  }
  legend("topright",sleg,col = 1:length(sleg),lty=1:length(sleg))
   dev.off()
   }

