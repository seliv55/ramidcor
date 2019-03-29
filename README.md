# simid
To extract mass isotopomer distribution from CDF files saved by mass spectrometers in SIM mode

![Logo](figs/logo.png)

# SIMID
Version: 1.0

## Short description
R-program to read CDF files containing time course of mass spectra of 13C-labeled metabolites, and write the extracted spectra in a format appropriate for further analysis.

## Description
SIMID is a computer program designed to read the machine-generated files saved in netCDF format containing registered time course of m/z chromatograms. It evaluates the mass isotopomer distribution (MID) at the moment when peaks are reached, and saves the obtained information in a table, making it ready for further correction for natural isotope occurrence.
SIMID is written in “R”, uses library “ncdf4” (it should be installed before the first use of siMID)  and contains several functions, located in the files “simid.R” and "lib.R", designed to read cdf-files, and analyze and visualize the spectra that they contain. The functionality of SIMID is similar to that of RaMID, with the difference that it does not need the previously prepared table with a list of CDF files and additional information, but uses only a short description of conditions provided by data owners.

## Key features
- primary processing of 13C mass isotopomer data obtained with GCMS

## Functionality
- Preprocessing of raw data
- initiation of workflows of the data analysis

## Approaches
- Isotopic Labeling Analysis / 13C
    
## Instrument Data Types
- MS

## Data Analysis
siMID reads the CDF files presented in the working directory, and then
- separates the time courses for selected m/z peaks corresponding to specific mass isotopomers;
- corrects baseline for each selected mz;
- choses the time points where the distribution of peaks is less contaminated by other compounds and thus is the most representative of the real analyzed distribution of mass isotopomers;
- evaluates this distribution, and saves it in files readable by MIDcor, a program, which performs the next step of analysis, i.e. correction of the RaMID spectra for natural isotope occurrence, which is necessary to perform a fluxomic analysis.

## Screenshots
- screenshot of input data (format Metabolights), output is the same format with one more column added: corrected mass spectrum

![screenshot]()

## Tool Authors
- Vitaly Selivanov (Universitat de Barcelona)

## Container Contributors
- [Pablo Moreno](EBI)

## Website
- N/A

## Git Repository
- https://github.com/seliv55/simid

## Installation

- As independent program, SIMID itself does not require installation.  There are two ways of using it: either creating a library "simid", or reading source files containing the implemented functions. Standing in the sim2mid directory:

- 1) Create a library of functions:
   
```
 sudo R

 library(devtools)
 
 build() 
 
 install() 
 
 library(simid) 
 
 library(ncdf4)
```

- 2) read directly the necessary functions:
  
```
 R 
 
 source("R/simid.R") 
 
 source("R/midcor.R") 
 
 source("R/lib.R") 
 
 source("R/stadyn.R") 
 
 library(ncdf4)
```

- a directory should contain the .cdf files that are to be analyzed.

## Usage Instructions

- The programs support three subsequent steps of preparation of raw MS data for the fluxomic analysis. 
1. Extraction of raw mass spectra from NetCDF files saved by a mass spectrometer, containing recordings in selected intervals of m/z. It is performed using SIMID executing the  command (default arguments are shown):

```
  metan(iinfile="sw620",cdfdir="SW620/",fiout="out.csv")
```
- here the parameters are (i) the path to an input file with a short description of input data explained in detail below; (ii) the path to a directory containing the .CDF files with raw mass spectrometer data, i.e. registration of the injections into the mass spectrometer performed in the course of the given analyzed experiment; (iii) a path to the directory for the output results (extracted relative intensities for all m/z constituting the mass spectra (or mass isotopomer distributions (MID)) of the metabolites of interests).

2. The file containing the results provided by SIMID (here "out.csv") can be used to fill the database Metabolights. The other outputs prepared by SIMID and saved in "files/SW620/" are made convenient for visual checking and further correction for natural isotope abundance by MIDcor.
 
```
 rumidcor(infile="sw620",dadir="files/SW620/")
```
3. Then the corrected data can be checked and wrong data edited/eliminated manually. The tool Stadyn reads the checked files performs simple statistical analysis for repeated measurements and prepares the data for simulation.
```
 isoform(isofi='toIsodyn',dor="SW620/",marca=3)
```
## The format of input data description

- The input data description file (here "sw620") contains the additional information prepared by the data provider that is necessary for the analysis and for the output table to write in the format accepted as exchangeable with the Metabolights database. It contains the following columns: 
(i) names of metabolites of interest, which spectra should be extracted from the provided CDF files; 
(ii) retention time (RT); 
(iii) m/z value of the lightest isotopomer (mz0) corresponding to the resolved derivatized fragment of metabolite of interest; 
(iv) position of the resolved carbon fragment in the parent molecule; 
(v) chemical formula of the derivatized compound containing the given fragment; 
(vi) m/z value of the lightest isotopomer corresponding to another fragment of the same metabolite (control).

## The names of CDF files provided

- The input file (here "sw620") provides the general informarion used by the program. Moreover, the names of CDF files should contain a specific information referred to each separate measurement. Here is an example of filename: "SW620_6h_12Glc_R1_PIM_SIM_01.CDF". SW620 is the type of analyzed cells 6h is the time of incubation 12Glc indicates the artificially labeled substrate applied. R1 is a number of biological replicate 01 is a number of ingection to MS machine from the same biological replicate

Based on this information and that extracted from the CDF files presented in the working directory siMID evaluates the mass spectra of the metabolites listed in "metdata", and saves it in tables accepted as exchangeable with Metabolights database.

## An example provided

- Run the provided example using the command:

```
  metan(infile="sw620",cdfdir="SW620/",fiout="out.csv")
```


