# Program to combine various outputs of Midcor to a single input for Isodyn
msdlist<-function(trati){nln<-length(trati)
        ntrati<-lapply(strsplit(trati, " "),as.numeric)
        lnumv<-!lapply(ntrati,is.na)[[1]]
        mdis<-matrix(ncol=length(ntrati[[1]][lnumv]),nrow=nln,0)
        mdis[1,]<-ntrati[[1]][lnumv]
      if(nln>1)for(i in 2:nln){
        lnumv<-!lapply(ntrati,is.na)[[i]]
        mdis[i,]<-ntrati[[i]][lnumv]
      } 
        mval<-round(apply(mdis,2,mean),4)
        sdv<-round(apply(mdis,2,sd))
 return(list(mval,sdv))}
 
 vybor<-function(spis,obr){
 	lspis<-grepl(obr,spis)
 return(spis[lspis]) }
 
 isoform<-function(isofi='toIsodyn',dor="SW620/",marca=3){
        a<-readLines(isofi)
# basic data:  
  fnam<-strsplit(a[1],' ')[[1]] # metabolite(file) name
  tinc<-strsplit(a[2],' ')[[1]] # incubation times
  trac<-strsplit(a[3],' ')[[1]] # tracer used
   trr<-trac[marca]; metm<-paste(c(tinc),collapse=" ")
    tmp4<-paste(c('tracer',a[marca+2]),collapse=" ")
   
 for(met in fnam[2:length(fnam)]){
   a<-readLines(paste(dor,strsplit(met,',')[[1]][1],'_c',sep=''))
   beg<-grepl(' corrected',a)
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
        metm<-c(metm,tmp2,tmp3)
  }   }  } }
        metm<-c(metm,tmp4)
  write(metm,paste("mark",marca,sep=""))
   return (metm)}
