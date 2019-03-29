#  setwd(oldi)
# oldi<-getwd()
#  source(infile)
# print(infile)
#metan(infile="tere1",cdfile='teresa/86stds_mix.cdf',fiout="out.csv")
    library(ncdf4)
icms<-function(infile="tere1",cdfile='teresa/Refsample_12C_Glucose_ICMS.cdf',fiout="out"){
   start.time <- Sys.time()
   outdir="files/"
   intab<-read.table(infile,header=T)
   nc <- nc_open(cdfile, readunlim=FALSE)  #open cdf file
   ret<-ncvar_get( nc, "scan_acquisition_time" )
   twide<-50
   rett<-list();
   for(i in 1:nrow(intab)) rett[[i]]<-getRtInt(retal=ret,rt=intab$RT[i],tlim=twide)
   remove(ret)
   mzz<-ncvar_get( nc, "mass_values" )
   mzbeg<-ncvar_get( nc, "scan_index" )+1 # index of beginning of mz scan interval for each timepoint
   mzpos<-list(); nCfrg<-numeric()
   for(i in 1:nrow(intab)){
     frag<-as.character(intab$Fragment[i]) # string of fragnent name
     frpos<-gregexpr("C[0-9]",frag)[[1]]+1 # C-positions in the fragment
     c1=as.numeric(substr(frag,frpos[1],frpos[1]));  c2=as.numeric(substr(frag,frpos[2],nchar(frag)))
     nCfrg[i]<-c2-c1+2 # number of mass isotopomers in fragment
     isopos<-list();
     for(j in 1:nCfrg[i]) # get valid mz for the selected scans and mass isotopomers
       isopos[[j]]<-getMzInt(mzal=mzz,mzalbeg=mzbeg,tl=rett[[i]][[2]],numscans=twide*2,lmz=intab[i,3], hmz=intab[i,4],iso=(j-1))
     mzpos[[i]]<-isopos
   }
   remove(mzz)
   ival<-ncvar_get( nc, "intensity_values" )
   peaks<-list(); ivabs<-list();  ivrel<-list(); alinfo<-character()
   for(i in 1:nrow(intab)) {
      peak<-list(); disabs<-numeric(); disrel<-numeric();
      for(j in 1:nCfrg[i]) {
      peak[[j]]<-round(geIVsum(iv=ival,mzgood=mzpos[[i]][[j]]))
      ima<-which.max(peak[[j]]);  ma<-peak[[j]][ima]
      if((ima<7)|(ima>twide-7)){disabs[j]<-0; next} 
          disabs[j]<-basecor(peak[[j]],ima,ma) # correct for baseline
          alinfo<-paste(c(alinfo,j-1,'#',ima,':',ma,'\n'),collapse=" ")
      } # sum of intensities for good mz 
      disabs[is.na(disabs)]<-0; ivabs[[i]]<-c(as.character(intab[i,1]),round(disabs))
      sabs<-sum(disabs); disrel<-disabs/sabs; ivrel[[i]]<-c(as.character(intab[i,1]),round(disrel,4))
      
      peaks[[i]]<-peak; alinfo<-paste(c(alinfo,ivabs[[i]],'\n',ivrel[[i]],'\n\n'),collapse=" ")
   }
   remove(ival)
        write(alinfo,fiout)
  print(Sys.time() - start.time)
  return(peaks)
  }
#             rm(list=ls(all=TRUE))

 dc<-1.003355
 dh<-1.006277
 dn<-0.997035
 do1<-1.004217
 do2<-2.004245
 ds1<-0.9993878
 ds2<-1.995796
 chekiso<-function(delta,l0,h0,nc=3,niso=1){
 ((l0+delta)<(h0+niso*dc))&((h0+delta)>(l0+niso*dc))
 }
 
