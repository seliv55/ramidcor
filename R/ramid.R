#  setwd(oldi)
# oldi<-getwd()
#  source('R/lib.R')
# print(infile)
#ramid(infile="../filescamid/sw620",cdfdir="../filescamid/SW620/",fiout="out.csv",md='scan')
# ramid(infile="../INES/ScanList.csv",cdfdir="../INES/PIM/KO_Hypoxia/KO_Hypoxia_SCANLAC.AIA/",fiout="out.csv",md='scan')  
# ramid(infile="../INES/SimList.csv",cdfdir="../INES/PIM/Parental_Hypoxia/Parental_Hypoxia_SIMLAC.AIA/",fiout="out.csv",md='sim')
#ramid(infile='../johanna/MediaHELNormoxia/MediaHELNorLong',cdfdir='../MediaHELNormoxia/Media HEL Nor Long/',fiout="out.csv",md='scan')
#ramid(infile='../johanna/MediaHELHypoxia/Media HEL Hyp long.txt',cdfdir='../johanna/MediaHELHypoxia/Media HEL Hyp long/',fiout="out.csv",md='scan')

 library(ncdf4)
ramid<-function(infile="../filesimid/sw620",cdfdir="../filescamid/SW620/",fiout="out.csv",md='scan'){
   if(md=='uhr') {a<-icms(); return(a)}
   start.time <- Sys.time()
   pat=".CDF"
   lcdf<-dir(path = cdfdir,pattern=pat)
   outdir="files/"
   intab<-read.table(infile,header=T)

title<-ftitle()

     df0<-data.frame(); # data frame to write Ramid output in PhenoMeNal format
     res<-character(); res1<-character(); res2<-character(); phen<-""
     if(md=='scan') for(fi in lcdf){ # fi <- lcdf[1]
            fi<-paste(cdfdir,fi, sep="");     fi1<-fi
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
        write(paste('###\t\t\tRAMID version 1.0',Sys.time(),'\t***\n'),ofi)
          tmp<-subset(res,(grepl(as.character(nam),res)))
          if(length(tmp)){
            mzrow<-subset(tmp,(grepl("mz:",tmp)))
           mzs<- strsplit(strsplit(mzrow[1],"c: ")[[1]][1],"mz: ")[[1]][2];
           titl<-paste("absolute_values,CDF_file: Max",mzs)
           titl1<-paste("relative_values,CDF_file: Max",mzs)
          tmp<-gsub(as.character(nam)," ",tmp)
          tmp<-gsub("  ","",tmp)
        write(tmp,ofi,append=T)
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
                deltac<-round(pikintc-basc*5)
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
     delta<-round(pikint-bas*5); s5tp<-"5_timepoints:"
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

rmbadpik<-function(mat,vekc,piklim,ls){
    lpikpos<-which.max(mat[,2]); rpikpos<-lpikpos; len<-length(mat[,2]) 
if(rpikpos>(4*len/5)) {
      mat<-mat[-(rpikpos-piklim):-len,]
      lpikpos<-which.max(mat[,2]); rpikpos<-lpikpos
      } else { if(lpikpos<(len/5)) {mat<-mat[-1:-(rpikpos+piklim),] 
       vekc<-vekc[-1:-(rpikpos+piklim)]
            lpikpos<-which.max(mat[,2]); rpikpos<-lpikpos
      }      }
    while(mat[rpikpos,2]>ls) if(rpikpos<len) {rpikpos=rpikpos+1} else {break}
return(list(lpikpos,rpikpos,mat,vekc))}

 
discan<-function(fi,intab, tlim=20){
# fi: file name
# intab: parameters of metabolite (mz for m0, retention time)
    limsens<-8000000

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
 for(imet in 1:nrow(intab)) if(max(rett)>rts[imet]){#imet<-23#729fcf
   nm<-as.character(intab$Name[imet])
   itpeak<-which(abs(rett-rts[imet])<tlim)  # index of timepoint closest to theoretical retention time
   ltpeak<-length(itpeak)
   tpeak<-rett[itpeak]
   
   imzi<-numeric()
   imzi[1]<-sum(npoint[1:(itpeak[1]-1)])+1      # index of m/z corresponding to left time boundary of peak
  for(i in 2:(ltpeak+1))  imzi[i]<-imzi[i-1]+npoint[itpeak[i-1]]
   imzfi<-imzi[ltpeak+1] # index of m/z corresponding to right time boundary of peak
   mzpeak<-mz[imzi[1]:imzfi] # all mz between the timepoints limiting the peak
   ivpeak<-iv[imzi[1]:imzfi] # all intensity between the timepoints limiting the peak 
  for(i in 2:(ltpeak+1))  imzi[i]<-imzi[i]-imzi[1]+1;  imzi[1]<-1

# main peak
        a<-as.character(intab$Fragment[imet])
    nCfrg<-as.numeric(substr(a,4,nchar(a)))-as.numeric(substr(a,2,2))+1
    nmassm <-nCfrg+4 # number of isotopores to present calculated from formula
   a<-psimat(nr=ltpeak, nmass=nmassm, imzi, mzpeak, ivpeak, mzz0=mz0[imet], dmzz=dmz, ofs=2)
    intensm<-a[[2]]; selmzm<-a[[3]]; intensm[is.na(intensm)]<-0; selmzm[is.na(selmzm)]<-0;
      
    piklim<-7;  nmassc<-1
   a<-psimat(nr=ltpeak, nmass=nmassc, imzi, mzpeak, ivpeak, mzz0=mzcon[imet], dmzz=dmz, ofs=1)
    intensc<-a[[2]]; selmzc<-a[[3]]; intensc[is.na(intensc)]<-0; selmzc[is.na(selmzc)]<-0;

      a<-rmbadpik(mat=intensm,vekc=intensc,piklim, ls=limsens)
       pikposm<-a[[1]]; pikposr<-a[[2]]; intensm<-a[[3]]; intensc<-a[[4]]
       ltpeak<-length(intensm[,2])
      
      if((pikposm<piklim)|(pikposm>(length(intensm[,2])-piklim))) next
       bs<-basln(intensm[,2])
      if(intensm[pikposm,2]< 5*bs) next
# control peak
         pikposc<-which.max(intensc[(pikposm-piklim):(pikposm+piklim)])
        
  if(abs(pikposc-piklim)>5) next
       maxpikc<-intensc[pikposc]
       maxpikm<-max(intensm[pikposm,])
       isomax<-which.max(intensm[pikposm,])
     pikmzm<-numeric(); pikintm<-numeric();  basm<-numeric(); ilim=0
  for(k in 1:nmassm) {
      pikmzm[k]<-selmzm[pikposm,k] # peak mz
      if(intensm[pikposm,k]>limsens) { ilim=ilim+1;
         while(intensm[pikposm+ilim,k]>limsens) ilim=ilim+1
         pikposr<-pikposm+ilim
         pikintm[k]<-((intensm[(pikposm-1),k]+intensm[pikposr,k]))/2
      }  else pikintm[k]<- sum(intensm[(pikposm-2):(pikposm+2),k])/5
   }
      basm<-apply(intensm,2,basln,pos=pikposm,ofs=9)
  delta<-round(pikintm-basm)
       smaxpik<-"max_peak:";
  if(maxpikm>limsens) {smaxpik<-"**** !?MAX_PEAK:"; print(paste("** max=",maxpikm,"   ",nm,"   **"))}
        s5tp<-"5_timepoints:"
  if(delta[1]/delta[2] > 0.05) { s5tp<-"*!?* 5_timepoints:"; 
      print(paste("+++ m-1=",delta[1],"  m0= ",delta[2],"   +++ ",nm))}
    a<- wphen(fi,nm,intab$Fragment[imet], intab$Formula[imet], intab$RT[imet], pikmzm,delta)
    phenom<-c(phenom,a)
                rat<-delta/basm
                rel<-round(delta/max(delta),4)      # normalization
    pikposc<-pikposm-piklim+pikposc-1;
   archar<-paste(c(gsub(' ','_',fi),nm),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"RT(min):",round(tpeak[pikposm]/60,3),"c:",round(tpeak[pikposc]/60,3)),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,smaxpik,maxpikm,"c:",maxpikc),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"mz:",round(pikmzm,1)),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,s5tp,pikintm),collapse=" ")
         result<-c(result,archar)
   archar<-paste(c(nm,"base:",round(basm)),collapse=" ")
         result<-c(result,archar)
#   archar<-paste(c(nm,"max-base:",delta),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(gsub(' ','_',fi),nm,maxpikm,delta),collapse=" ")
         res1<-c(res1,archar)
#   archar<-paste(c(nm,"relative:",rel),collapse=" ")
#         result<-c(result,archar)
   archar<-paste(c(gsub(' ','_',fi),nm,maxpikm,rel),collapse=" ")
         res2<-c(res2,archar)
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

