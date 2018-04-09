
basln<-function(vec,pos=length(vec),ofs=0){# baseline
   basl<--1; basr<--1;bas<-0
  if(pos>ofs) basl<-mean(vec[1:(pos-ofs)])
  if(pos<(length(vec)-ofs)) basr<-mean(vec[(pos+ofs):length(vec)])
  if((basl>0)&(basr>0)) bas<-min(basl,basr)
  else if(basl<0) bas<-basr
  else if(basr<0) bas<-basl
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
    
