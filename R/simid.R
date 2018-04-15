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

title<-ftitle()

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
  if(ratc[goodiso]>3){ a<-as.character(intab$Fragment[i])
    nCfrg<-as.numeric(substr(a,4,nchar(a)))-as.numeric(substr(a,2,2))+1
#                   a<-as.character(intab$Formula[i])
#          Cpos<-regexpr("C",a); Hpos<-regexpr("H",a); Spos<-regexpr("S",a); Sipos<-regexpr("Si",a)
#    nCder<-as.numeric(substr(a,Cpos+1,Hpos-1))
#    if(Sipos>0) nSi<-as.numeric(substr(a,Sipos+2,Sipos+2)) else nSi<-0
#    if((Spos>0)&(Spos!=Sipos)) nS<-as.numeric(substr(a,Spos+1,Spos+1)) else nS<-0
        nmass<-nCfrg+5 # number of isotopomers to present calculated from formula
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
#     { }      else pikmz<-c(0,pikmz)
    if((misofin[1]==pikmz[1])&(delta[1]/delta[2] > 0.075)) { s5tp<-"*!?* 5_timepoints:";
      print(paste("+++ m-1=",delta[1],"  m0= ",delta[2],"   +++ ",nm)); break }
          
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
        if(misofin[1]!=pikmz[1]) {delta<-c(0,delta); rel<-c(0,rel)}
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

