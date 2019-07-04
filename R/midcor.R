#source("R/lib.R")
#source("R/midcor.R")
#correct("GluC2C4b","var")
case2<-function(tmp,mmlab,corr,ff,fr,nfrg){ nln=nrow(tmp); nmass=ncol(tmp)-2
    for(j in 1:nln) { k=1; for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k]);
               for(k in 2:(nfrg+1)) for(i in 1:(nmass+1-k)) tmp[j,i+1+k]<-tmp[j,i+1+k]-corr[i]*(fr[j,1+k])*ff;#
    }
      fr<-mdistr(nfrg,tmp,mmlab,nln); return(fr)
}

corrm<-function(tmp,mmlab,corr,ff,fr,nfrg){ nln=nrow(tmp); nmass=ncol(tmp)-2;
     k=1; for(i in 1:(nmass+1-k)) tmp[i+1+k]<-tmp[i+1+k]-corr[i]*(fr[1+k]);
   for(k in 2:(nfrg+1)) for(i in 1:(nmass+1-k)) tmp[i+1+k]<-tmp[i+1+k]-corr[i]*(fr[1+k])*ff;#
       return(mdistr(nfrg,tmp,mmlab,1))
}

fitf<-function(tmp,mmlab,corr,ff,fr,nfrg){ nln=nrow(tmp); nmass=ncol(tmp)-2;
  dpred=numeric(nfrg); dpred[2]=1.; xisq=0.; modf=1.01; ff1=ff
    fr<-corrm(tmp,mmlab,corr,ff,fr,nfrg);   xisq0=sum((fr[2:(1+nfrg)]-dpred)**2)
   iw=0;   while(iw<2){ ff1=ff1*modf; print(paste("iw=",iw))
    fr<-corrm(tmp,mmlab,corr,ff,fr,nfrg);   xisq=sum((fr[2:(1+nfrg)]-dpred)**2)
                if(xisq<xisq0) {ff=ff1; }
                 else {modf=1./modf;  iw=iw+1}
      }
      
       return(list(ff,xisq,fr))
}
#midcor(infile="../johanna/MediaHELNormoxia/MediaHELNorLong",dadir="../johanna/files/Media HEL Nor Long/")
#midcor(infile='../johanna/MediaHELHypoxia/Media HEL Hyp long.txt',dadir='../johanna/files/Media HEL Hyp long/')
midcor<-function(infile,dadir){
   a <-read.table(infile, skip=1,nrows=3); info<-as.character(a[,2])
   intab<-read.table(infile, skip=4,header=T); phenom<-""
   for(imet in 1:nrow(intab)){ #  imet<-6
   met <- as.character(intab$Name[imet])
    fn<-paste(dadir,met,sep="")
    if(file.exists(fn)&file.size(fn)>500){
      con<-file(fn,open="r")
      line<-readLines(con)
      dati<-which(grepl("absolute_v",line))
      mzi<-as.numeric(strsplit(line[dati]," ")[[1]][-1:-2])
      datf<-which(grepl("END",line))
      df0<-NULL
      for(ln in line[(dati+1):(datf-1)]){
         a<-strsplit(ln,' ')[[1]]
         rbind(df0,a)->df0
      }
      line=""; close(con)
      phenom<- c(phenom,correct(fn,dfi=df0,mzi,metdat=intab[imet,],info)) 
            }    }
         fiout<-"fenform"               
       write(ftitle(),fiout)
       write(phenom,fiout,append=T)
}

getcor<-function(dfic,dfil,frcold,gcold,glab,name,mmlab,fil,nfrg,mzis) {
       if(!is.matrix(dfic)) {frcold<-t(frcold); gcold<-t(gcold); dfic<-t(dfic); }
       if(!is.matrix(dfil)) {glab<-t(glab); dfil<-t(dfil); }
        dficold<-as.numeric(dfic[,2]); out=''
        dfilab<-as.numeric(dfil[,2]);
      scor<-paste('CF_m',0:nfrg,sep='')
    df<-lapply(dfilab,'-',dficold) #
    adf<-lapply(df,abs)
    madf<-sapply(adf,which.min,simplify=T)
    for(i in 1:length(dficold)) {
         ils<-which(madf==i); if(length(ils)<1) next
    if(frcold[i,1]<0.95) print(paste(name,': m0 in unlabeled=',round(frcold[i,1],4)))
    corr<-mmlab[1,]-gcold[i,1:ncol(mmlab)]
    gci<-glab[ils,]
    if(length(ils)==1) {tmp<-gci[1:ncol(mmlab)]+corr; nln=1;  } else {
    tmp<-apply(gci[,1:ncol(mmlab)],1,'+',corr);
    }
    tmp<-t(tmp)
     fr<-mdistr(nfrg,tmp,mmlab,nln);
     if(length(ils)==1) res<-cbind(t(dfil[ils,]),t(round(fr[,1:(nfrg+1)],4))) else
                        res<-cbind(dfil[ils,],round(fr[,1:(nfrg+1)],4))
     row1<-cbind(t(dfic[i,]),t(round(frcold[i,1:(nfrg+1)],4)))
     res<-rbind(row1,res)
     out<-rbind(out,res)
  write("\n*** Correction factor: **",fil,append=TRUE)
  write(paste(scor,collapse=' '),fil,append=T)
  write.table(format(t(corr[1:(nfrg+1)]),digits=4),fil,quote=F,append=T,col.names=F, row.names = F);
  write(paste(mzis,collapse=' '),fil,append=T);
  write.table(res,fil,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F); 
   }
   return(out)}
   

correct<-function(fn,dfi,mzi,metdat,info){#fname is the name of file with raw data;
# samb,samf,cndb,cndf are the positions of the initial and final characters in the row name
# designating the biological sample and conditions correspondingly
 frag<-as.character(metdat$Fragment)
 formula<-as.character(metdat$Formula)
 pCf<-gregexpr("C",frag)[[1]]
       iCb<- as.numeric(substr(frag,pCf[1]+1,pCf[1]+1))
       iCe<- as.numeric(substr(frag,pCf[2]+1,pCf[2]+1))
       nfrg<- iCe-iCb+1
    Cpos<-regexpr("C",formula); Hpos<-regexpr("H",formula);
    nC<-as.numeric(substr(formula,Cpos+1,Hpos-1))
      nSi<-0
    Sipos<-regexpr("Si",formula); 
    if(Sipos>0) { nnSi<-as.numeric(substr(formula,Sipos+2,Sipos+2))
                  if(is.na(nnSi)) nSi<-1 else nSi<-nnSi   }  
      nS<-0
    Spos<-regexpr("S[1-9]",formula)
    if(Spos>0) {
    nS=as.numeric(substr(formula,Spos+1,Spos+1))
      } else { if((regexpr("S[A-Z]",formula)>0) | (regexpr("S$",formula)>0)) nS=1 }    
    
      coln<-ncol(dfi)-2;  nln<-nrow(dfi); nmass<-coln-1;
      gcmss<-matrix(nrow=nln,ncol=coln)
      
      gcmss[,]<-dfi[,-1:-2]
      dfi<-dfi[,1:2]
      gcms<-NULL
    for(i in 1:coln) cbind(gcms, as.numeric( gcmss[,i]))->gcms
# theoretic distribution:
          mmlab<-mtr(nfrg,nmass,nC,nSi,nS);
# normalization
    if(metdat$mz0>=mzi[2]) {ef<-sum(gcms[,1])/sum(gcms[,2]); gcm<-elim(gcms,ef)
        } else gcm<-as.matrix(cbind(gcms[,-1],0)); 
    gcmsn<-gcm
    gcmsn[1:nln,]<-gcm[1:nln,]/apply(gcm,1,sum)
# mass fractions
  if(nln>1) fr<-mdistr(nreal=nfrg,msd=gcmsn[,1:ncol(mmlab)],mm=mmlab,nln) else {
      fr<-mdistr(nreal=nfrg,msd=t(as.matrix(gcmsn[,1:ncol(mmlab)])),mm=mmlab,nln)
  }
  
 fn1<-paste(fn,".txt",sep="");
 mzis<-c("Sample_file Max_intensity",  paste('m',0:nfrg, sep=''))
 write(paste('###\tMIDCOR version 1.0:',tools::md5sum('R/midcor.R'),Sys.time(),'\t***\n'),fn1)
 write("*** MID for each injection, corrected for natural enrichment using empirical formula ***",fn1,append=T)
 write(paste(mzis,collapse=' '),fn1,append=T)
 if(nln>1) write.table(cbind(dfi,round(fr[,1:(nfrg+1)],4)),fn1,quote=F,append=T,col.names=F, row.names = F) else {
  write.table(cbind(t(dfi),round(t(fr[,1:(nfrg+1)]),4)),fn1,quote=F,append=T,col.names=F, row.names = F)
 }
 if(nln==1) {print(paste(metdat$Name,': Only one sample provided')); return()}
# correction
     cold<- grep(info[2],dfi[,1])
 if(length(cold)<1) {print(paste(metdat$Name,': No unlabeled provided')); return()}
     incold<-which(as.numeric(dfi[cold,2])<as.numeric(info[1]))
    if(length(incold)<1) {print(paste(metdat$Name,': Max in unlabeled >',info[1])); return()}
  write("\n*** Samples additionally corrected for possible matrix effects **",fn1,append=TRUE)
     cold<-cold[incold]
     frcold<-fr[cold,]; gcold<-gcmsn[cold,]; dfic<-dfi[cold,]
#     if(!is.matrix(dfic)){frcold<-t(frcold); gcold<-t(gcold); dfic<-t(dfic)}
     frlab<-fr[-cold,]; glab<-gcmsn[-cold,]; dfil<-dfi[-cold,]
#     if(!is.matrix(dfil)){frlab<-t(frlab); glab<-t(glab); dfil<-t(dfil)}
     thr<-as.numeric(info[3])
     hicold<-integer(); locold<-integer(); hilab<-integer(); lolab<-integer()
     if(length(cold)>1){
     hicold<-which(as.numeric(dfic[,2])>thr)
     hilab<-which(as.numeric(dfil[,2])>thr)
     locold<-which(as.numeric(dfic[,2])<=thr)
     lolab<-which(as.numeric(dfil[,2])<=thr)
     }
  if((length(hicold)==0)|(length(locold)==0)|(length(cold)==1)) 
       getcor(dfic=dfic,dfil=dfil,frcold=frcold,gcold=gcold,glab=glab,name=metdat$Name,mmlab=mmlab,fil=fn1,nfrg,mzis) else {
  if(length(hilab)==0) getcor(dfic[locold,], dfil, frcold[locold,], gcold[locold,], glab, metdat$Name, mmlab, fil=fn1,nfrg,mzis)
    else {if(length(lolab)==0)
     getcor(dfic[hicold,],dfil, frcold[hicold,], gcold[hicold,], glab, metdat$Name, mmlab, fil=fn1,nfrg,mzis) else {
     getcor(dfic[hicold,],dfil[hicold,], frcold[hicold,], gcold[hicold,], glab[hicold,], metdat$Name, mmlab, fil=fn1,nfrg,mzis)
     getcor(dfic[locold,],dfil[locold,], frcold[locold,], gcold[locold,], glab[locold,], metdat$Name, mmlab, fil=fn1,nfrg,mzis)
     }
     }
     }
}

 isoform<-function(isofi='../filesimid/files/toIsodyn',dadir='../cdf2mid/files/cdfcase3/',marca=2){
        a<-readLines(isofi)
# basic data:  
  fnam<-strsplit(a[1],' ')[[1]] # metabolite(file) name
  tinc<-strsplit(a[2],' ')[[1]] # incubation times
  trac<-strsplit(a[3],' ')[[1]] # tracer used
   trr<-trac[marca]; metm<-paste(c(tinc),collapse=" ")
    tmp4<-paste(c('tracer',a[marca+2]),collapse=" ")
   
 for(met in fnam[2:length(fnam)]){
   a<-readLines(paste(dadir,strsplit(met,',')[[1]][1],'.txt',sep=''))
   beg<-grep(' corrected',a)
   suba<-a[beg[length(beg)]:length(a)]
     tr<-vybor(suba,trr) #select tracer
  if(length(tr)>0){
  metm<-c(metm,paste("name:",met));
   for(tii in tinc[2:length(tinc)]){
      tm<-paste(tii,'h',sep="")
      ttr<-vybor(tr,tm) #select incubation time
  if(length(ttr)>0){
  a<-msdlist(ttr)
        tmp2<-paste(c("t=",tii,"mean:",a[1][[1]],sum(a[1][[1]])),collapse=" ")
        tmp3<-paste(c("sd:",a[2][[1]]),collapse=" ")
        tmp3<-gsub("NA","0",tmp3)
        metm<-c(metm,tmp2,tmp3)
  }   }  } }
        metm<-c(metm,tmp4)
  write(noquote(metm),paste(dadir,"mark",marca,sep=""))
   return (metm)}
