#source("lib.R")
#source("midcor.R")
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

rumidcor<-function(infile="sw620",dadir="files/SW620/"){
   intab<-read.table(infile,header=T,sep=" "); phenom<-""
   for(imet in 1:nrow(intab)){
   met <- as.character(intab$Name[imet])
    fn<-paste(dadir,met,sep="")
    if((file.exists(fn))&(file.size(fn)>1000)){
      con<-file(fn,open="r")
      line<-readLines(con)
      dati<-which(grepl("absolute_v",line))
      mzi<-as.numeric(strsplit(line[dati]," ")[[1]][-1:-2])
      datf<-which(grepl("relative_v",line[dati:length(line)]))
      line=""; close(con)
      con<-file(fn,open="rt")
      df0<-read.table(con,skip=(dati-1),header=T,nrows=(datf-4))
      close(con)
      phenom<- c(phenom,correct(fn,df0,mzi,intab[imet,])) 
            }    }
         fiout<-"fenform"               
       write(ftitle(),fiout)
       write(phenom,fiout,append=T)
}

correct<-function(fn,dfi,mzi,metdat){#fname is the name of file with raw data;
# samb,samf,cndb,cndf are the positions of the initial and final characters in the row name
# designating the biological sample and conditions correspondingly
 frag<-as.character(metdat$Fragment)
 formula<-as.character(metdat$Formula)
 nfrg<-as.numeric(substr(frag,4,nchar(frag)))-as.numeric(substr(frag,2,2))+1
   Cpos<-regexpr("C",formula); Hpos<-regexpr("H",formula);
   Spos<-regexpr("S",formula); Sipos<-regexpr("Si",formula)
    nC<-as.numeric(substr(formula,Cpos+1,Hpos-1))
    if(Sipos>0) nSi<-as.numeric(substr(formula,Sipos+2,Sipos+2)) else nSi<-0
    if((Spos>0)&(Spos!=Sipos)) nS<-as.numeric(substr(formula,Spos+1,Spos+1)) else nS<-0
      gcms<-dfi[-1:-2]
            coln<-ncol(gcms);  nln<-nrow(gcms); nmass<-coln-1;
# theoretic distribution:
          mmlab<-mtr(nmass,coln,nC,nSi,nS);
# normalization
    if(metdat$mz0>=mzi[2]) {ef<-sum(gcms[1])/sum(gcms[2]); gcm<-elim(gcms,ef)
        } else gcm<-as.matrix(cbind(gcm[,-1],0))
    gcmsn<-gcm[1:nln,]/apply(gcm,1,sum)
# mass fractions
   fr<-mdistr(nmass,gcmsn,mmlab,nln);# write mass fractions without correction:
 fn1<-paste(fn,"_c",sep="");
 write("*** MID for each injection, corrected only for natural 13C, 29,30Si, 33,34S ***",fn1)
  write.table(format(cbind(dfi[1],fr),digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
# correction
           corr<-numeric(nmass+1); corr1=numeric(nmass+1); icomm=0;
    lcon<-grep('cold',dfi[,1])
if(length(lcon)>0){
 if(length(lcon)>1){
    ncon<-which.min(abs(1-fr[lcon,1]))
    corr<-mmlab[1,]-gcmsn[lcon,][ncon,] 
 } else if(length(lcon)==1){
  corr<-mmlab[1,]-gcmsn[lcon,]
 }
     tmp<-t(apply(gcmsn,1,'+',corr)); 
     fr<-mdistr(nmass,tmp,mmlab,nln);
     res<-cbind(dfi[1],fr)
 write("\n*** Samples fully corrected **",fn1,append=TRUE)
    write.table(format(res,digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F); 
 write("*** Correction factor: **",fn1,append=TRUE)
     write.table(format(t(corr),digits=4),fn1,quote=FALSE,append=TRUE,col.names=FALSE, row.names = F);
}
  phen<-""; fr<-cbind(0,fr[,-ncol(fr)]); gc<-cbind(0,gcms[,-ncol(gcms)])
  for(i in 1:nrow(dfi)){
    a<- wphen(as.character(dfi[i,1]),metdat$Name,frag, formula, metdat$RT, mzi,round(gcms[i,]),fr[i,])
    phen<-c(phen,a)
  }
return(phen)}
