- 1) Create a library of functions:
   
```
 sudo R
 library(devtools)
 build() 
 install() 
 library(sim2mid) 
 library(ncdf4)
```
- 2) read directly the necessary functions:
  
```
 R 
 source("R/simid.R") 
 source("R/libsim.R") 
 library(ncdf4)
 
  metan(infile="sw620",cdfdir="SW620/",outdir="files/")
```
 source("R/midcor.R") 
 source("R/lib.R") 
 
 correct("filename",samb=1,samf=3,cndb=4,cndf=9)

 
  
```

The file containing the results provided by cdf2mid (here "cdf2midout.csv") can be used by RaMID, or directly proceed for further correction by MIDcor.

