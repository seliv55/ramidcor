#  setwd(oldi)
# oldi<-getwd()
#  source(infile)
# print(infile)
#metan(infile="../filescamid/sw620",cdfdir="../filescamid/SW620/",fiout="out.csv",md='scan')
# metan(infile="../INES/ScanList.csv",cdfdir="../INES/PIM/KO_Hypoxia/KO_Hypoxia_SCANLAC.AIA/",fiout="out.csv",md='scan')  
# metan(infile="../INES/SimList.csv",cdfdir="../INES/PIM/Parental_Hypoxia/Parental_Hypoxia_SIMLAC.AIA/",fiout="out.csv",md='sim')
 library(ncdf4)
metan<-function(infile="../filesimid/sw620",cdfdir="../filescamid/SW620/",fiout="out.csv",md='scan'){
   if(md=='uhr') {a<-icms(); return(a)}
   start.time <- Sys.time()
   pat=".CDF"
   lcdf<-dir(path = cdfdir,pattern=pat)
   outdir="files/"
   intab<-read.table(infile,header=T)

title<-ftitle()

     df0<-data.frame(); # data frame to write Ramid output in PhenoMeNal format
     res<-character(); res1<-character(); res2<-character(); phen<-""
     if(md=='scan') for(fi in lcdf){
                                    fi<-paste(cdfdir,fi, sep="")
                                    a <-discan(fi,intab) 
                                    res<-c(res,a[[1]])
                                    res1<-c(res1,a[[2]])
                                    res2<-c(res2,a[[3]])
                                    phen<-c(phen,a[[4]])
       } else for(fi in lcdf){
                              fi<-paste(cdfdir,fi, sep="")
                              a <-getdistr(fi,intab)
                              res<-c(res,a[[1]])
                              res1<-c(res1,a[[2]])
                              res2<-c(res2,a[[3]])
                              phen<-c(phen,a[[4]])
       }
       write(title,fiout)
       write(phen,fiout,append=T)
       a<-strsplit(cdfdir,"/")[[1]];  len<-length(a)
       if(len==1) { if(!(file.exists(outdir))) dir.create(outdir)
           celdir<-paste(outdir,a[length(a)],"/",sep='') } else if(len>1){ celdir<-paste("../",a[2],"/",outdir,sep='')
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
       
getdistr<-function(fi,intab, tlim=100){
# fi: file name
# intab: parameters of metabolite (mz for m0, retention time)
    a<-readcdf(fi); 
     result<-character(); res1<-character(); res2<-character()
     mz<-round(a[[1]],1); iv<-a[[2]]   # all mz and respecive intensities
     npoint<-a[[3]];      rett<-a[[4]] # number of mz points and respective rt
#     totiv<-a[[5]]                     # sum of intensities at each rt
    a<-strsplit(strsplit(fi,".CDF")[[1]][1],"/")[[1]];
    fi<-a[length(a)]; print(fi);  phenom<-""
     rts<-intab$RT*60.; mz0<-round(intab$mz0,1); mzcon<-round(intab$control,1)
#  search for specified metabolites
 for(imet in 1:nrow(intab)) {nm<-as.character(intab$Name[imet])
      nmass<-3; rtdev<-15; 
#  cutting the desired peak inside tlim
        tpclose<-which.min(abs(rett-rts[imet]))   # index (scan number) of retention time for the peak
        tplow<-tpclose-tlim; tpup<-tpclose+tlim-1 # time boundaries that include desired peak
   mzi<-sum(npoint[1:tplow])+1              # index of m/z corresponding to left time boundary of peak
   mzpik<-sum(npoint[(tplow+1):tpup])              # number of mz values inside the time boundaries
   mzfi<-mzi+mzpik-1                        # index of m/z corresponding to right time boundary of peak
   misoc<-c(intab$control[imet],intab$control[imet]+1,intab$control[imet]+2)  #desired mz values
    lmisoc<-mz[mzi:mzfi] %in% misoc               # logical array indicating the positions of desired mz 
    pic<-which(lmisoc>0)                          # array of positions of desired mz
    if(length(pic)<3) next
#    checking if the left boundary of actual desired mz registration < tlim
    ibmz<-npoint[tplow+1]; ibmz0<-0; i<-1;
    beg<-pic[1];
    while(ibmz<beg) { i<-i+1; ibmz0<-ibmz; ibmz<-ibmz+npoint[tplow+i]}
    tplow<-tplow+i-1; mzi<-mzi+ibmz0; tpclose<-tlim-i+1
#    checking if the right boundary of actual desired mz registration < tlim
    ifmz0<-length(mz[mzi:mzfi]); ifmz<-ifmz0-npoint[tpup-1]
    fin<-pic[length(pic)]; i<-1
    while(ifmz>fin) { i<-i+1; ifmz0<-ifmz; ifmz<-ifmz-npoint[tpup-i]} 
    tpup<-tpup-i+1; mzfi<-mzi+ifmz0-1
#    putting the intensities of desired mz into a matrix (scan number, mz)    
    lmisoc<-mz[mzi:mzfi] %in% misoc               # logical array indicating the positions of desired mz 
    pic<-which(lmisoc>0)                          # array of positions of desired mz
    intens<-matrix(ncol=nmass,nrow=(tpup-tplow),0)
    intens<-sweep(intens,2,iv[mzi:mzfi][lmisoc],'+')
#    determine positions of peaks for each mz, their intensities and the most correct position
    pospiks<-apply(intens,2,which.max)
    pikintc<-apply(intens,2,max)
    if(max(abs(diff(pospiks)))>9) goodiso<-which.min(abs(pospiks-tpclose))  else goodiso<-which.max(pikintc)
        pikposc<-pospiks[goodiso]
        
  if(abs(rett[pikposc]-rett[tpclose])<10) { # actual rt deviates from the expected in < 10 sec.
    maxpikc<-pikintc[goodiso]
    for(k in 1:nmass) pikintc[k]<-sum(intens[(pikposc-2):(pikposc+2),k])
      basc<-apply(intens,2,basln,pos=pikposc,ofs=11)
                deltac<-round(pikintc-basc)
                ratc<-deltac/basc
# main peak
   if(ratc[goodiso]>3){ frag<-as.character(intab$Fragment[imet])
     frpos<-gregexpr("C[0-9]",frag)[[1]]+1
     c1=as.numeric(substr(frag,frpos[1],frpos[1]));  c2=as.numeric(substr(frag,frpos[2],nchar(frag)))
     nCfrg<-c2-c1+1
     nmass<-nCfrg+5       # number of desired isotopomers calculated from formula
    misodes<-array((mz0[imet]-1):(mz0[imet]+nmass-2)) # desired isotopomers 
    lmiso<-misodes %in% mz[mzi:(mzi+npoint[tplow+1])] # do they actually present?
    miso<-misodes[lmiso]                              # actually presented isotopomers
    lalmiso<-mz[mzi:mzfi] %in% miso # right positions of actually presented isotopomers
    
    nmass<-length(miso) 
    intens<-matrix(ncol=nmass,nrow=(tpup-tplow),0)
    intens<-sweep(intens,2,iv[mzi:mzfi][lalmiso],'+') # create matrix iv(col=mz,row=rt) that includes the peak

    pospiks<-apply(intens,2,which.max)
    difpos<-abs(pikposc - pospiks)
    if(min(difpos)<3) {
     goodpos<-pospiks[which.min(difpos)]
     pikint<-intens[goodpos,]
     pikpos<-which.max(pikint)
     maxpik<-intens[goodpos,pikpos]; smaxpik<-"max_peak:";
     if(maxpik>8300000) {smaxpik<-"**** !?MAX_PEAK:"; print(paste("** max=",maxpik,"   ",nm,"   **")); next;}
     bas<-apply(intens,2,basln,pos=goodpos,ofs=11)
     if((goodpos>2)&(goodpos<(nrow(intens)-2))){
       for(k in 1:nmass) pikint[k]<-sum(intens[(goodpos-2):(goodpos+2),k]) }
     delta<-round(pikint-bas); s5tp<-"5_timepoints:"
    if((delta[1]/delta[2] > 0.075)) { s5tp<-"*!?* 5_timepoints:";
      print(paste("+++ m-1=",delta[1],"  m0= ",delta[2],"   +++ ",nm)); next }
          
                rat<-delta/bas
                rel<-round(delta/max(delta),4)      # normalization

    a<- wphen(fi,nm,intab$Fragment[imet], intab$Formula[imet], intab$RT[imet], miso,delta)
    phenom<-c(phenom,a)

   archar<-paste(c(gsub(' ','_',fi),nm),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"peak_index:",goodpos,"c:",goodpos),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,smaxpik,maxpik,"c:",maxpikc),collapse=" ") 
         result<-c(result,archar)
   archar<-paste(c(nm,"mz:",miso,"c:",misoc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,s5tp,pikint,"c:",pikintc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"base:",round(bas),"c:",round(basc)),collapse=" ")
         result<-c(result,archar)
#   archar<-paste(c(nm,"max-base:",delta),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(gsub(' ','_',fi),nm,maxpik,delta),collapse=" ")
         res1<-c(res1,archar)
#   archar<-paste(c(nm,"relative:",rel),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(gsub(' ','_',fi),nm,maxpik,rel),collapse=" ")
         res2<-c(res2,archar)
    }
    
#     { }      else miso<-c(0,miso)
                }
              }
           }
 return(list(result,res1,res2,phenom))}     
 
discan<-function(fi,intab, tlim=50){
# fi: file name
# intab: parameters of metabolite (mz for m0, retention time)
    a<-readcdf(fi); 
     result<-character(); res1<-character(); res2<-character()
     mz<-round(a[[1]],1); iv<-a[[2]]   # all mz and respecive intensities
     npoint<-a[[3]];      rett<-a[[4]] # number of mz points and respective rt
#     totiv<-a[[5]]                     # sum of intensities at each rt
    a<-strsplit(strsplit(fi,".CDF")[[1]][1],"/")[[1]];
    fi<-a[length(a)]; print(fi);  phenom<-""
     rts<-intab$RT*60.; mz0<-round(intab$mz0,1); mzcon<-round(intab$control,1)
#     totiv<-a[[5]]                     # sum of intensities at each rt
      dmz=0.49  # index of beginning of mz scan interval for each timepoint
#  search for specified metabolites
 for(imet in 1:nrow(intab)) if(max(rett)>rts[imet]){nm<-as.character(intab$Name[imet])
   tpclose<-which.min(abs(rett-rts[imet]))  # index of timepoint closest to theoretical retention time
   tplow<-tpclose-tlim; tpup<-tpclose+tlim # indexes of peak boundaries
   rtpeak<-rett[tplow:tpup] # retention times within the boundaries

   mzi<-sum(npoint[1:tplow])+1              # index of m/z corresponding to left time boundary of peak
   mzpik<-sum(npoint[(tplow+1):tpup])              # number of mz values inside the time boundaries
   mzfi<-mzi+mzpik-1                        # index of m/z corresponding to right time boundary of peak
# additional peak
   mzpeak<-mz[mzi:mzfi] # all mz between the timepoints limiting the peak
   ivpeak<-iv[mzi:mzfi] # all intensity between the timepoints limiting the peak
      nmass<-3; rtdev<-15; intens<-matrix(); selmz<-matrix()
   a<-psimat(nr=(length(rtpeak)+2), nmass, mzpeak, ivpeak, mzz0=mzcon[imet], dmzz=dmz, lefb=1, rigb=length(rtpeak)-1, ofs=1)
    intens<-a[[2]]; selmz<-a[[3]]; intens[is.na(intens)]<-0; selmz[is.na(selmz)]<-0;
  pikmzc<-numeric(); 
    pikintc<-apply(intens,2,max)
    pospiks<-apply(intens,2,which.max)
   if(max(abs(diff(pospiks)))>9) goodiso<-which.min(abs(pospiks-tlim))  else goodiso<-which.max(pikintc)
        pikposc<-which.max(intens[,goodiso])
  if(abs(pikposc-tlim)<rtdev) {
       maxpikc<-intens[pikposc,goodiso]
   for(k in 1:nmass) {
   pikmzc[k]<-selmz[pikposc,k] # peak mz
   pikintc[k]<-sum(intens[(pikposc-2):(pikposc+2),k])}
      basc<-apply(intens,2,basln,pos=pikposc,ofs=5)
                ratc<-round(pikintc-basc)/basc
# main peak
  if(ratc[goodiso]>9){  piklim<-9
        a<-as.character(intab$Fragment[imet])
    nCfrg<-as.numeric(substr(a,4,nchar(a)))-as.numeric(substr(a,2,2))+1
    nmass <-nCfrg+5 # number of isotopores to present calculated from formula
#        a<-as.character(intab$Formula[imet])
#     Cpos<-regexpr("C",a);  Hpos<-regexpr("H",a)
#     Spos<-regexpr("S",a); Sipos<-regexpr("Si",a)
#    nCder<-as.numeric(substr(a,Cpos+1,Hpos-1))
#    if(Sipos>0) nSi<-as.numeric(substr(a,Sipos+2,Sipos+2)) else nSi<-0
#    if((Spos>0)&(Spos!=Sipos)) nS<-as.numeric(substr(a,Spos+1,Spos+1)) else nS<-0
#       tit1<-paste("m",c(1:nmass),sep="")
    a<-psimat(nr=(2*piklim+1),nmass,mzpeak,ivpeak,mzz0=mz0[imet],dmzz=dmz,lefb=(pikposc-piklim),rigb=(pikposc+piklim),ofs=2)
    intens<-a[[2]]; selmz<-a[[3]];
#	lapply(intens,length)
    pikint<-apply(intens,2,max)
    isomax<-which.max(pikint)
    pikpos<-which.max(intens[,isomax])
    maxpik<-intens[pikpos,isomax];  smaxpik<-"max_peak:";
     if(maxpik>8000000) {smaxpik<-"**** !?MAX_PEAK:"; print(paste("** max=",maxpik,"   ",nm,"   **"))}
    bas<-apply(intens,2,basln,pos=pikpos,ofs=5)
  if((pikpos>2)&(pikpos<(nrow(intens)-2))){
      pikmz<-numeric(); 
     for(k in 1:nmass) {
     pikmz[k]<-round(selmz[pikpos,k]) # peak mz
     pikint[k]<-sum(intens[(pikpos-2):(pikpos+2),k])
   }
    delta<-round(pikint-bas); s5tp<-"5_timepoints:"
     if(delta[1]/delta[2] > 0.05) { s5tp<-"*!?* 5_timepoints:"; print(paste("+++ m-1=",delta[1],"  m0= ",delta[2],"   +++ ",nm))}
    a<- wphen(fi,nm,intab$Fragment[imet], intab$Formula[imet], intab$RT[imet], pikmz,delta)
    phenom<-c(phenom,a)
                rat<-delta/bas
                rel<-round(delta/max(delta),4)      # normalization
    pikpos<-pikposc-piklim+pikpos-1;
   archar<-paste(c(gsub(' ','_',fi),nm),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"peak_index:",pikpos,"c:",pikposc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,smaxpik,maxpik,"c:",maxpikc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"mz:",round(pikmz,1),"c:",round(pikmzc,1)),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,s5tp,pikint,"c:",pikintc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"base:",round(bas),"c:",round(basc)),collapse=" ")
         result<-c(result,archar)
#   archar<-paste(c(nm,"max-base:",delta),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(gsub(' ','_',fi),nm,maxpik,delta),collapse=" ")
         res1<-c(res1,archar)
#   archar<-paste(c(nm,"relative:",rel),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(gsub(' ','_',fi),nm,maxpik,rel),collapse=" ")
         res2<-c(res2,archar)
            }
         }
       } # if additional peak exists
     } # change of metabolite (imet)
 return(list(result,res1,res2,phenom))}
  
icms<-function(infile="tere1",cdfile='teresa/Refsample_12C_Glucose_ICMS.cdf',fiout="out"){
   start.time <- Sys.time()
   outdir="files/"
   intab<-read.table(infile,header=T)
   nc <- nc_open(cdfile, readunlim=FALSE)  #open cdf file
   ret<-ncvar_get( nc, "scan_acquisition_time" )
   twide<-50
   rett<-list();
   for(imet in 1:nrow(intab)) rett[[imet]]<-getRtInt(retal=ret,rt=intab$RT[imet],tlim=twide)
   remove(ret)
   mzz<-ncvar_get( nc, "mass_values" )
   mzbeg<-which(diff(mzz)<0)+1 # index of beginning of mz scan interval for each timepoint
   mzpos<-list(); nCfrg<-numeric()
   for(imet in 1:nrow(intab)){
     frag<-as.character(intab$Fragment[imet]) # string of fragnent name
     frpos<-gregexpr("C[0-9]",frag)[[1]]+1 # C-positions in the fragment
     c1=as.numeric(substr(frag,frpos[1],frpos[1]));  c2=as.numeric(substr(frag,frpos[2],nchar(frag)))
     nCfrg[imet]<-c2-c1+2 # number of carbons in fragment
     isopos<-list();
     for(j in 1:nCfrg[imet]) # get valid mz for the selected scans and mass isotopomers
       isopos[[j]]<-getMzInt(mzal=mzz,mzalbeg=mzbeg,tl=rett[[imet]][[2]],numscans=twide*2,lmz=intab[imet,3], hmz=intab[imet,4],iso=(j-1))
     mzpos[[imet]]<-isopos
   }
   remove(mzz)
   ival<-ncvar_get( nc, "intensity_values" )
   peaks<-list(); ivabs<-list();  ivrel<-list(); alinfo<-character()
   for(imet in 1:nrow(intab)) {
      peak<-list(); disabs<-numeric(); disrel<-numeric();
      for(j in 1:nCfrg[imet]) {
      peak[[j]]<-round(geIVsum(iv=ival,mzgood=mzpos[[imet]][[j]]))
      ima<-which.max(peak[[j]]);  ma<-peak[[j]][ima]
      if((ima<7)|(ima>twide-7)){disabs[j]<-0; next} 
          disabs[j]<-basecor(peak[[j]],ima,ma) # correct for baseline
          alinfo<-paste(c(alinfo,j-1,'#',ima,':',ma,'\n'),collapse=" ")
      } # sum of intensities for good mz 
      disabs[is.na(disabs)]<-0; ivabs[[imet]]<-c(as.character(intab[imet,1]),round(disabs))
      sabs<-sum(disabs); disrel<-disabs/sabs; ivrel[[imet]]<-c(as.character(intab[imet,1]),round(disrel,4))
      peaks[[imet]]<-peak; alinfo<-paste(c(alinfo,ivabs[[imet]],'\n',ivrel[[imet]],'\n\n'),collapse=" ")
   }
   remove(ival)
        write(alinfo,fiout)
  a<-Sys.time() - start.time
  return(list(peaks,a))
  }

