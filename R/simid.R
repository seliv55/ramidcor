#  setwd(oldi)
# oldi<-getwd()
#  source(infile)
# print(infile)
metan<-function(infile="sw620",cdfdir="SW620/",fiout="out.csv"){ #infile="metdata",  main function; evaluates MID for a set of CDF files specified by pat
#  and metabolite specified by m/z M-1 referred as ms
# call: metan()
# temp <- tempdir()#paste(,"/",sep="")  #"data/ttt/"  #
# lcdf<-unzip(cdfzip,exdir=temp)
#   setwd(cdfdir)
   start.time <- Sys.time()
   pat=".CDF"
   lcdf<-dir(path = cdfdir,pattern=pat)
   outdir="files/"
   intab<-read.table(infile,header=T,sep=" ")
     
title<-paste("Raw_Data_File", "cells", "tracer_molecule", "labelled_positions","abundance", "injection","Replicate", "Incubation_time", "Metabolite_name", "CHEBI","atomic_positions", "Empirical_formula", "retention(min)", "mz_monitored", "signal_intensity", "isotopologue", "isotologue_abundance")

     df0<-data.frame(); # data frame to write Ramid output in PhenoMeNal format
     res<-character(); res1<-character(); res2<-character(); phen<-""

       for(fi in lcdf){ fi<-paste(cdfdir,fi, sep="")
     a <-getdistr(fi,intab);
           res<-c(res,a[[1]]); res1<-c(res1,a[[2]]); res2<-c(res2,a[[3]]); phen<-c(phen,a[[4]])
       }
      
#       outfile<-paste(outdir,"all_info",sep="")
       write(title,fiout)
       write(phen,fiout,append=T)
       a<-strsplit(cdfdir,"/")[[1]];  len<-length(a)
       if(len==1) { if(!(file.exists(outdir))) dir.create(outdir)
       celdir<-paste(outdir,a[length(a)],"/",sep='') }
         else if(len>1){ celdir<-paste("/",a[len-1],"/",outdir,sep='')
           if(!(file.exists(celdir))) dir.create(celdir)
        celdir<-paste(celdir,a[len],"/",sep='')
       } 
        print(paste(celdir," len=",len))
       if(!(file.exists(celdir))) dir.create(celdir)
       
       for(nam in intab$Name) {ofi<-paste(celdir,nam,sep="")
          tmp<-subset(res,(grepl(as.character(nam),res)))
          if(length(tmp)){
            mzrow<-subset(tmp,(grepl("mz:",tmp)))
           mzs<- strsplit(strsplit(mzrow[1],"c: ")[[1]][1],"mz: ")[[1]][2];
           titl<-paste("absolute_values,CDF_file: Max",mzs)
           titl1<-paste("relative_values,CDF_file: Max",mzs)
          tmp<-gsub(as.character(nam)," ",tmp)
          tmp<-gsub("  ","",tmp)
        write(tmp,ofi)
          tmp<-subset(res1,(grepl(paste(as.character(nam)),res1)))
          tmp<-gsub(as.character(nam)," ",tmp)
          tmp<-gsub("  ","",tmp)
          tmp<-c(titl,tmp)
        write("\n",ofi,append=T)
        write(tmp,ofi,append=T)
          tmp<-subset(res2,(grepl(paste(as.character(nam)),res2)))
          tmp<-gsub(as.character(nam)," ",tmp)
          tmp<-gsub("  ","",tmp)
          tmp<-c(titl1,tmp)
        write("\n",ofi,append=T)
        write(tmp,ofi,append=T)
        }}
   Sys.time() - start.time
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
  
wphen<-function(fi,nm,fragg, formul, rtt, pikmz,delta){
      a<-strsplit(fi,'_')[[1]]
      cel<-a[1]; incub<-substr(a[2],1,nchar(a[2]))
      trac<-switch(a[3],
                  "UGln"=c("[U-C13]-Glutamine","1,1,1,1,1",100),
                  "12Glc"=c("[1,2-C13]-Glucose","1,1,0,0,0,0",50),
                  default=c(0,0,0))
      repl<-a[4]; inj<-a[length(a)]
      nikiso<-paste("m",c(-1:(length(pikmz)-2)),sep="")
  return(paste(fi,cel,trac[1],trac[2],trac[3],repl,inj,incub,nm,"chebi",fragg, formul, rtt, pikmz,delta,nikiso))
   }
   
getdistr<-function(fi,intab){
# fi: file name
# intab: parameters of metabolite (mz for m0, retention time)
    a<-readcdf(fi); 
     result<-character(); res1<-character(); res2<-character()
     mz<-round(a[[1]],1); iv<-a[[2]]   # all mz and respecive intensities
     npoint<-a[[3]];      rett<-a[[4]] # number of mz points and respective rt
     totiv<-a[[5]]                     # sum of intensities at each rt
    a<-strsplit(strsplit(fi,".CDF")[[1]][1],"/")[[1]];
    fi<-a[length(a)]; print(fi);  phenom<-""
#    summary: 
 a<-info(mz,iv,npoint); mzpt<-a[[1]] # number of points in all the mz ranges
  tpos<-a[[2]]# time points of beginning of registration corresponding mz ranges
  mzind<-a[[3]]                      # corresponding mz indexes
  mzrang<-a[[4]]                     # all registered mz ranges 
     rts<-intab$RT*60.; mz0<-round(intab$mz0,1); mzcon<-round(intab$control,1)
#  search for specified metabolites
 for(i in 1:nrow(intab)) {nm<-as.character(intab$Name[i]); tlim=50
        ltp<- rts[i]<rett[tpos]         # time interval that includes rts
        ranum<-(c(1:length(tpos))[ltp])[1]-1
   if((mz0[i] %in% mzrang[[ranum]])&(mzcon[i] %in% mzrang[[ranum]])) {
#   check whether mid for a given metabolite is presented in the found time interval
        tpclose<-which.min(abs(rett-rts[i]))
#        tlim<-min(tlim,tpclose-tpos[ranum],tpos[ranum+1]-tpclose)
        tplow<-max(tpclose-tlim,tpos[ranum]); tpup<-min(tpclose+tlim,tpos[ranum+1])     # boundaries that include desired peak
   mzi<-mzind[ranum]+mzpt[ranum]*(tplow-tpos[ranum])#index of initial mz point
   mzfi<-mzi+mzpt[ranum]*(tpup-tplow)   	#index of final mz point
   rtpeak<-rett[tplow:tpup] # retention times within the boundaries
        tpclose<-which.min(abs(rtpeak-rts[i]))
# additional peak
      nmass<-3; rtdev<-15;
    misoc<-c(intab$control[i],intab$control[i]+1,intab$control[i]+2)#desired mz values
    lmisoc<-mzrang[[ranum]] %in% misoc
    intens<-matrix(ncol=nmass,nrow=(tpup-tplow),0)
    intens<-sweep(intens,2,iv[mzi:mzfi][lmisoc],'+')
    pospiks<-apply(intens,2,which.max)
    pikintc<-apply(intens,2,max)
   if(max(abs(diff(pospiks)))>9) goodiso<-which.min(abs(pospiks-tpclose))  else goodiso<-which.max(pikintc)
        pikposc<-pospiks[goodiso]
  if((abs(pikposc-tpclose)<rtdev)&(pikposc>2)&(pikposc<(nrow(intens)-2))) {
        maxpikc<-pikintc[goodiso]
    for(k in 1:nmass) pikintc[k]<-sum(intens[(pikposc-2):(pikposc+2),k])
      basc<-round(apply(intens,2,basln,pos=pikposc,ofs=5))
                deltac<-round(pikintc-basc)
                ratc<-deltac/basc

# main peak
  if(ratc[goodiso]>9){ a<-as.character(intab$Fragment[i])
    nCfrg<-as.numeric(substr(a,4,nchar(a)))-as.numeric(substr(a,2,2))+1
                   a<-as.character(intab$Formula[i])
          Cpos<-regexpr("C",a); Hpos<-regexpr("H",a); Spos<-regexpr("S",a); Sipos<-regexpr("Si",a)
    nCder<-as.numeric(substr(a,Cpos+1,Hpos-1))
    if(Sipos>0) nSi<-as.numeric(substr(a,Sipos+2,Sipos+2)) else nSi<-0
    if((Spos>0)&(Spos!=Sipos)) nS<-as.numeric(substr(a,Spos+1,Spos+1)) else nS<-0
        isomax<-nCfrg+5 # number of isotopomers to present calculated from formula
    ppos<-which((mz0[i]-1)==mzrang[[ranum]]) # position of desired MID in the selected mzrang
    nmass<-min(isomax,length(mzrang[[ranum]])-ppos+1) # number of isotopores to present in the spectrum
    misofin<-array((mz0[i]-1):(mz0[i]+nmass-2)) # isotopores to present in the spectrum
    lmisofin<-mzrang[[ranum]] %in% misofin # do they are present in the given mzrang?
    pikmz<-mzrang[[ranum]][lmisofin] # extrat those that are present
    nmass<-length(pikmz)
    intens<-matrix(ncol=nmass,nrow=(tpup-tplow),0)
    intens<-sweep(intens,2,iv[mzi:mzfi][lmisofin],'+') # create matrix iv(col=mz,row=rt) that includes the peak
     piklim<-min(11,pikposc-1,nrow(intens)-pikposc-1)
    intens<-intens[(pikposc-piklim):(pikposc+piklim),]
    
 
    pikint<-apply(intens,2,max)
    isomax<-which.max(pikint)
    pikpos<-which.max(intens[,isomax])
    maxpik<-intens[pikpos,isomax]; smaxpik<-"max_peak:";
     if(maxpik>8300000) {smaxpik<-"**** !?MAX_PEAK:"; print(paste("** max=",maxpik,"   ",nm,"   **")); break;}
    bas<-apply(intens,2,basln,pos=pikpos,ofs=5)
  if((pikpos>2)&(pikpos<(nrow(intens)-2))){
     for(k in 1:nmass) pikint[k]<-sum(intens[(pikpos-2):(pikpos+2),k])
   }
    delta<-round(pikint-bas); s5tp<-"5_timepoints:"
     if(delta[1]/delta[2] > 0.05) { s5tp<-"*!?* 5_timepoints:";
      print(paste("+++ m-1=",delta[1],"  m0= ",delta[2],"   +++ ",nm)); break; }
                rat<-delta/bas
                rel<-round(delta/max(delta),4)      # normalization
    pikpos<-pikposc-piklim+pikpos-1;

    a<- wphen(fi,nm,intab$Fragment[i], intab$Formula[i], intab$RT[i], pikmz,delta)
    phenom<-c(phenom,a)

   archar<-paste(c(fi,nm),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"peak_index:",pikpos,"c:",pikposc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,smaxpik,maxpik,"c:",maxpikc),collapse=" ") 
         result<-c(result,archar)
   archar<-paste(c(nm,"mz:",pikmz,"c:",misoc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,s5tp,pikint,"c:",pikintc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"base:",round(bas),"c:",round(basc)),collapse=" ")
         result<-c(result,archar)
#   archar<-paste(c(nm,"max-base:",delta),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(fi,nm,maxpik,delta),collapse=" ")
         res1<-c(res1,archar)
#   archar<-paste(c(nm,"relative:",rel),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(fi,nm,maxpik,rel),collapse=" ")
         res2<-c(res2,archar)
                }
              }
           }
           }
 return(list(result,res1,res2,phenom))}
 

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
     

