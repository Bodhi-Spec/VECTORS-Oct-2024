#This is the main file that calls each individual metrics
#This file will source functions of each metric from their respective folders
#This file loads in packages
#This file loads in the input csv data
#This file outputs the results of each called metric

#____________________Read in Packages_________________
library(ggplot2)
library(ggpmisc)
library(splus2R)

#_________ Define code and data paths_____________________________________
codedir<-"" #define path to the folder that contain the code files of the other functions
datadir<-"" #define path to the folder that contians the data (so the folder with the COORD files and the converted CSV files)


#_______Source functions from other files______________
codename_list<-list.files(codedir,pattern="\\.R") #list of all names of R scripts in code folder that contian the indivudal emtrics
codename_list<-codename_list[codename_list!="GetData.R"] #remove this current script from the list if it is in it
for (codename in codename_list){ #source each function
  source(paste0(codedir,codename))
  print(codename)
}

#__________Get the list of filenames of the data___________
filename_list<-c(list.files(datadir,pattern="\\.COORD"),list.files(datadir,pattern="\\.coord"))#first step to find list of filename which are the csv files we are reading ins
filename_list<-gsub("\\.COORD","",filename_list)
filename_list<-gsub("\\.coord","",filename_list)#list of filename's



#__________________MAIN FUNCTION___________________
#Call function providing specific filename and purpose
#Function will read the csv files coresponding to the filename
#Then calls the function that that preforms the metric corresponding to the purpose
#This function will output the results
#For specific format of the output check the code of the function

getdata<-function(filename,purpose){ #purpose is the purpose of this function
  
  #Read in csv files and format the data for the given filename
  spectra_data<<-read.csv(paste0(datadir,filename,"_coord.csv"))#read spectra data which is the spectral datapoints
  spectra_data$Residual<<-spectra_data$RawData-spectra_data$ProcessedData #Calculate a residual
  maxamp <- max(subset(spectra_data, ppm >= 1.5 & ppm <= 4)$ProcessedData) #To calibrate the relative amplitude scale, use the highest peak from 1.5ppm to 4.0 ppm
  spectra_data$RelativeAmp<<-(spectra_data$ProcessedData)/(maxamp) #Make relative amplitude (highest point in graph will have RelativeAmp=1)
  spectra_data$RelativeBase<<-(spectra_data$Background)/(maxamp) #Make relative baseline measure
  
  met_data<<-read.csv(paste0(datadir,filename,"_Metabolites_coord.csv"))#read Metabolites Data (containing CRLBs, SD, etc)
  misc_data<<-read.csv(paste0(datadir,filename,"_misc_coord.csv"))#read Miscelanous data (containing SNR, FWHM, etc)
  #note R will display the csv files' headers with weird notation for special characters. Make sure to check with csv file if unsure of header in data frames.
  
  
  #Calls other functions to preform metric given purpose
  if (purpose=="Plot"){
    return(plotgraph(spectra_data,met_data,misc_data,filename))
  }else if (purpose=="DuplicatePeaks"){
    return(DuplicatePeaks(spectra_data,met_data,misc_data,filename))
  }else if (purpose=="Glx_mI_Peaks"){
    return(append(Glxpeaks(spectra_data),mIpeaks(spectra_data)))
  }else if (purpose =="VerticalShifts"){
    return(append(anyNegative(spectra_data,met_data,misc_data,filename),belowBaseline(spectra_data,met_data,misc_data,filename)))
  }
}







#OPTIONAL EXAMPLE TASKS________________________________________________________:
#list of purposes (parameter in the getdata() function) and the metrics reported:
#purpose="VerticalShifts" --> reports anyNegative and belowBaseline metric
#purpose="DuplicatePeaks" --> reports existDuplicate metric and graphs the spectrum with duplicate peaks outlined
#purpose = "Glx_mI_Peaks" --> reports Glx Merge, Glx Distinct, mI Merge, and mI Distinct
#purpose = "Plot" --> plots the spectrum for visualization



#To find the doublet function (detecting movement artifact) for the first file:
getdata(filename_list[1],"DuplicatePeaks")

#To plot all the data in the folder and compile the spectra in a pdf:
outputpath<- "" #define output path
plotcompilation<-function(filename_list,outputpath){
  pdf(outputpath,width=12,height=8) #create pdf, can specify dimensions
  for (filename in filename_list){ #iterates the function for every spectra in the big folder
    print(getdata(filename,"Plot"))
  }
  dev.off()
}
plotcompilation(filename_list, outputpath)


