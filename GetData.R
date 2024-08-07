#ABout this file:
#This is the main file. All operations are called from this file
#This file will source functions from other files.
#This file loads in pacakges, loads in the input data, and outputs the result
#this file calls operations for one input file, or can compile the results of the operation for many input files in a folder

#____________________Read in Packages___________
library(ggplot2)
library(ggpmisc)
library(nortest)
library(car)
library(splus2R)
library(progress)
library(stringr)
library(e1071)
library(openxlsx2)
library(dplyr)
library(ggrepel)

#_________ FInd list of files from Folder_____________________________________
codedir<-"/Users/bodhiberoukhim/Documents/SummerInternship2023/QAProject/Code/" #directory with the code files

datadir<-"/Users/bodhiberoukhim/Documents/SummerInternship2023/QAProject/Data/COORD/" #directory with the coord data files SKYLERS DATA
#datadir<-"/Volumes/SPECTRO_PROC$/CCS_Proc/Study_Bipolar_7T/SVS-7T/LCModel_Bodhi_Subset/" #Subset Jessica Data (best pipeline)
#datadir<-"/Volumes/SPECTRO_PROC$/CCS_Proc/Study_Bipolar_7T/SVS-7T/LCModel/" #All Jessica Data
#datadir<-"/Users/bodhiberoukhim/Documents/SummerInternship2023/QAProject/Data/NEUROFIT_SVS_COORD/" #Neurofit SVS data
#datadir<-"/Users/bodhiberoukhim/Documents/SummerInternship2023/QAProject/Data/NEUROPRECISE_COORD/StudyNEUROPRECISE30ms/" #Neuroprecise 30ms
#datadir<-"/Users/bodhiberoukhim/Documents/SummerInternship2023/QAProject/Data/NEUROPRECISE_COORD/StudyNEUROPRECISE40ms/" #Neuroprecise 40ms
#datadir<- "/Users/bodhiberoukhim/Documents/SummerInternship2023/QAProject/Data/U01_CSI_COORD/" #U01 CSI

#Technical Paper Subset data dir
datadir<- "/Users/bodhiberoukhim/Documents/SummerInternship2023/TechnicalPaper/Data/NEUROFIT_SVS/" #NEUROFIT SVS Data
#datadir<- "/Users/bodhiberoukhim/Documents/SummerInternship2023/TechnicalPaper/Data/NEUROPRECISE30ms_SVS/" #NEUROPRECISE 30ms
#datadir<- "/Users/bodhiberoukhim/Documents/SummerInternship2023/TechnicalPaper/Data/NEUROPRECISE40ms_SVS/" #Neuroprecise 40ms

#Test Sample dataset:
datadir<- "/Users/bodhiberoukhim/Desktop/TestSubset/COORD/"
datadir<- "/Users/bodhiberoukhim/Desktop/ZhouCOORD/"

filename_list<-c(list.files(datadir,pattern="\\.COORD"),list.files(datadir,pattern="\\.coord"))#first step to find list of filename which are the csv files we are reading ins
filename_list<-gsub("\\.COORD","",filename_list)
filename_list<-gsub("\\.coord","",filename_list)#list of filename's

#set one filename_list to Skyler_filename_list and other to Jessica_filename_list for quick organization


#_______Source functions form other files
codename_list<-list.files(codedir,pattern="\\.R")#list all names of R scripts in COde folder bc they contain functions
codename_list<-codename_list[codename_list!="GetData.R"]#remove this current script from the lsit
codename_list<-codename_list[codename_list!="ACFplot.R"]
for (codename in codename_list){#source each R script
  source(paste0(codedir,codename))
  print(codename)
}



#_______Given filename+purpose, reads csv file, then calls other functions, sending the data_________
getdata<-function(filename,purpose){ #purpose is the purpose of this function
  spectra_data<<-read.csv(paste0(datadir,filename,"_coord.csv"))#read spectra data
  spectra_data$Residual<<-spectra_data$RawData-spectra_data$ProcessedData
  #To calibrate the relative amplitude scale, use the highest peak from 1.5ppm to 4.0 ppm. (tentative)
  maxamp <- max(subset(spectra_data, ppm >= 1.5 & ppm <= 4)$ProcessedData)
  spectra_data$RelativeAmp<<-(spectra_data$ProcessedData)/(maxamp) #Make relative amplitude unit with highest amplitude being 1
  spectra_data$RelativeBase<<-(spectra_data$Background)/(maxamp) #Make Relative Bbase line measure
  
  met_data<<-read.csv(paste0(datadir,filename,"_Metabolites_coord.csv"))#read Metabolites Data
  misc_data<<-read.csv(paste0(datadir,filename,"_misc_coord.csv"))#read Miscelanous data
  
  
  
  #note R will display the csv files' headers with weird notation for special characters. Make sure to check with csv file if unsure of header in data frames.
  
  #Calls other functions given purpose
  if (purpose=="Plot"){
    return(plotgraph(spectra_data,met_data,misc_data,filename))
  }
  else if (purpose=="Means"){
    wlx_p<-as.numeric(wilcoxon(spectra_data,met_data,misc_data,filename))
    tt_p<-as.numeric(ttest(spectra_data,met_data,misc_data,filename))
    return(list(wlx_p=wlx_p,tt_p=tt_p))
  }else if (purpose=="Normality"){
    shapiro_p<-shapiro(spectra_data,met_data,misc_data,filename)
    andersonD_p<-andersonD(spectra_data,met_data,misc_data,filename)
    plotnormality(spectra_data,met_data,misc_data,filename,shapiro_p,andersonD_p)
    return(list(shapiro_p=shapiro_p,andersonD_p=andersonD_p))
  }else if (purpose=="Autocorrelation"){
    DB_p<-(Durbinwatson(spectra_data,met_data,misc_data,filename,lag=20))
    return(DB_p)
  }else if (purpose=="Peaks"){
    return(met_peaks(spectra_data,met_data,misc_data,filename))
  }else if (purpose=="Double"){
    return(artifact_double(spectra_data,met_data,misc_data,filename))
  }else if (purpose=="ACF"){
    return(ACFfunc(spectra_data,met_data,misc_data,filename,lag=100))
  }else if (purpose=="Ghost"){
    return(ArtifactGhosting(spectra_data,filename))
  }else if (purpose=="Glx"){
    return(Glxpeaks(spectra_data))
  }else if (purpose=="Ins"){
    return(Inspeaks(spectra_data))
  }else if (purpose == "PlotGlxIns"){
    Glxdata<-(getdata(filename,"Glx"))# to test one file, input file name and purpose
    Glxscore=Glxdata$score
    Glxextrascore1=Glxdata$extrascore1
    gammascore=Glxdata$gammascore
    Insdata<-(getdata(filename,"Ins"))# to test one file, input file name and purpose
    Insrow=Insdata$row
    Insscore=Insdata$score
    Insextrascore1=Insdata$extrascore1
    return(plotgraphGlxIns(spectra_data,met_data,misc_data,filename,Glxscore,Glxextrascore1,gammascore,Insscore,Insextrascore1))
  }else if (purpose == "Peaksppmtest"){
    PeaksppmTest(spectra_data,met_data,misc_data,filename)
  }else if (purpose == "PlotPeaksOffset"){
    dataframe<-PeaksppmTest(spectra_data,met_data,misc_data,filename)
    PlotPeaksOffset(spectra_data,filename,dataframe)
  }else if (purpose == "Range"){
    rangetest(spectra_data,met_data,misc_data,filename)
  }else if (purpose == "BigLeftRight"){
    BigLeftRight(spectra_data,filename)
  }else if (purpose =="Negative"){
    return(append(artifact_negative_graph(spectra_data,met_data,misc_data,filename),artifact_negative_baseline(spectra_data,met_data,misc_data,filename)))
  }else if (purpose=="SNR_FWHM"){
    return(misc_data[c(1,2)])
  }else if (purpose=="PlotNegative"){
    return(PlotNegative(spectra_data,met_data,misc_data,filename))
  }else if (purpose=="CRLB"){
    Glu<-met_data$X.SD[which(met_data$Metabolite=="Glu")]
    Gln<-met_data$X.SD[which(met_data$Metabolite=="Gln")]
    Ins<-met_data$X.SD[which(met_data$Metabolite=="Ins")]
    Glx<-met_data$X.SD[which(met_data$Metabolite=="Glu+Gln")]
    return(list(Glu=Glu,Gln=Gln,Ins=Ins,Glx=Glx))
  }else if (purpose=="Gauss"){
    return(gaussiantest(filename,spectra_data,"Ins"))
  }
}

getdata(filename_list[40],"Plot")
getdata("S14-NEUROFIT_104-xx-5-yy-7","Plot")

filename_list[449:476]

set.seed(123)
pdf("/Users/bodhiberoukhim/Desktop/GaussIns.pdf",width=12,height=7)
pb<-progress_bar$new(total=length(filename_list[sample(1:7180,size=1000,replace=FALSE)]))
Gauss_df<-data.frame()
for (filename in filename_list[sample(1:7180,size=1000,replace=FALSE)]){
  row<-getdata(filename,"Gauss")
  Gauss_df<-rbind(Gauss_df,row)
  rownames(Gauss_df)[nrow(Gauss_df)]<-filename
  pb$tick()
}
dev.off()


Gauss_df[which(is.na(Gauss_df$AcheivedConvergenceTolerance_double)),1:7]<-NA
Gauss_df[which(is.na(Gauss_df$AcheivedConvergenceTolerance_single)),10:13]<-NA

write.csv(Gauss_df,"/Users/bodhiberoukhim/Desktop/GaussIns_df.csv")



pb<-progress_bar$new(total=length(filename_list))
crlb_df<-data.frame()
for (filename in filename_list){
  row<-getdata(filename,"CRLB")
  crlb_df<-rbind(crlb_df,row)
  pb$tick()
}

#PLot BIGLEFTRIGHT + make data frame
pdf("/Users/bodhiberoukhim/Desktop/BigLeftRight_plot.pdf",width=12,height=8)
bigleftright_df<-data.frame()
pb<-progress_bar$new(total=length(filename_list))
for (filename in filename_list){
  data<-getdata(filename,"BigLeftRight")
  print(data$plot)
  bigleftright_df<-rbind(bigleftright_df, data[2:5])
  rownames(bigleftright_df)[nrow(bigleftright_df)]<-filename
  pb$tick()
}
dev.off()
write.csv(bigleftright_df,"/Users/bodhiberoukhim/Desktop/BigLeftRight_df.csv")

#_____pekas Ppm test --> to see if a faze/doublet artifact shifts peak location:
peaksppmdf<-data.frame()
pb<-progress_bar$new(total=1444)
for (filename in filename_list){
  row<-getdata(filename,"Peaksppmtest")
  peaksppmdf<-rbind(peaksppmdf, row)
  rownames(peaksppmdf)[nrow(peaksppmdf)]<-filename
  pb$tick()
}
write.csv(peaksppmdf,"/Users/bodhiberoukhim/Desktop/PeaksPPMTest.csv")

#Plot the bad peaks:
pdf("/Users/bodhiberoukhim/Desktop/OffsetPeaks.pdf",width=12,height=8)
offsetpeaks<-filename_list[c(which(peaksppmdf$NAAppm<=1.998839),which(peaksppmdf$NAAppm>=2.022606))]
pb<-progress_bar$new(total=39)
for (filename in offsetpeaks){
  dataframe<-getdata(filename,"Peaksppmtest")
  print(getdata(filename,"PlotPeaksOffset"))
  pb$tick()
}
dev.off()

#Runs the Glx and Ins tests, makes a fodler with the plot, and two csv files with Glx darta and Ins data
pdf("/Users/bodhiberoukhim/Desktop/GlxInsData/GlxInsScoresplot.pdf",width=12,height=8)
Glxdf<-data.frame()
Insdf<-data.frame()
pb<-progress_bar$new(total=length(filename_list))
#baddata<-c("S68-NEUROFIT_57-xx-8-yy-4","S35-NEUROFIT_23-xx-6-yy-4","S21-NEUROFIT_111-xx-11-yy-6","S104-NEUROFIT_94-xx-10-yy-11", "S98-NEUROFIT_88-xx-4-yy-8","S39-NEUROFIT_27-xx-11-yy-4")  
#badglx<-c("S103-NEUROFIT_93-xx-11-yy-9" ,"S56-NEUROFIT_44-xx-11-yy-9" ,"S104-NEUROFIT_94-xx-9-yy-5" ,"S39-NEUROFIT_27-xx-11-yy-11" ,"S18-NEUROFIT_108-xx-5-yy-4"  , "S75-NEUROFIT_64-xx-10-yy-4"  ,"S104-NEUROFIT_94-xx-9-yy-4" ,  "S6-NEUROFIT_07-xx-11-yy-11"  ,"S39-NEUROFIT_27-xx-11-yy-4" ,  "S68-NEUROFIT_57-xx-8-yy-4"   ,"S35-NEUROFIT_23-xx-6-yy-4"  ,  "S21-NEUROFIT_111-xx-11-yy-6" ,"S104-NEUROFIT_94-xx-10-yy-11", "S98-NEUROFIT_88-xx-4-yy-8")
for (filename in filename_list){
  Glxdata<-(getdata(filename,"Glx"))# to test one file, input file name and purpose
  Glxrow=Glxdata$row
  Glxscore=Glxdata$score
  Glxextrascore1=Glxdata$extrascore1
  Glxextrascore2=Glxdata$extrascore2
  gammascore=Glxdata$gammascore
  Glxdf<-rbind(Glxdf,c(Glxscore,Glxextrascore1,Glxextrascore2,gammascore,Glxrow))
  rownames(Glxdf)[nrow(Glxdf)]<-filename
  colnames(Glxdf)<-c("score","extrascore1", "extrascore2","gammascore",seq(from = 3, to = 37, by = 2))
  
  Insdata<-(getdata(filename,"Ins"))# to test one file, input file name and purpose
  Insrow=Insdata$row
  Insscore=Insdata$score
  Insextrascore1=Insdata$extrascore1
  Insextrascore2=Insdata$extrascore2
  Insdf<-rbind(Insdf,c(Insscore,Insextrascore1,Insextrascore2,Insrow))
  rownames(Insdf)[nrow(Insdf)]<-filename
  colnames(Insdf)<-c("score","extrascore1", "extrascore2", seq(from = 3, to = 57, by = 2))
  pb$tick()
  print(getdata(filename,"PlotGlxIns"))
}
dev.off()
write.csv(Glxdf,"/Users/bodhiberoukhim/Desktop/GlxInsData/GlxScoresdf.csv")
write.csv(Insdf,"/Users/bodhiberoukhim/Desktop/GlxInsData/InsScoresdf.csv")



pdf("/Users/bodhiberoukhim/Desktop/RightSidePlot.pdf",width=12,height=8)
for (filename in listrightside){
  print(getdata(filename,"Plot"))# to test one file, input file name and purpose
}
dev.off()

pdf("/Users/bodhiberoukhim/Desktop/NEUROFIT_SVS_fazedublet.pdf",width=12,height=8)
pb<-progress_bar$new(total=length(filename_list))
for (filename in filename_list){
  print(getdata(filename,"Double")$plot)# to test one file, input file name and purpose
  pb$tick()
}
dev.off()


pdf("/Users/bodhiberoukhim/Desktop/DoublePlot.pdf",width=12,height=8)
for (filename in filename_list[which(DoubleFazedf$existduplicate!=summary$existDublet)]){
  plot<-getdata(filename,"Double")$plot
  print(plot)
}
dev.off()


#__________Function that compiles Durbin watson p-value, D-W test stat, and AutoC test stat for in one plot (per lag)____
pb<-progress_bar$new(total=400)
AC_p<-data.frame()
AC_dw<-data.frame()
AC_ac<-data.frame()
ACplotcomp<-function(filename_list){
  for (filename in filename_list[sample(seq(1,length(filename_list)),100)]){#here I'm only choosing 400 random samples
    for (i in c(1,2,3,4)){
      AutoC_triad<-getdata(filename,"Autocorrelation",i)
      n<-length(AutoC_triad)
      AutoC_row<-AutoC_triad[1:(n/3)]
      AC_p<-rbind(AC_p, AutoC_row)
      rownames(AC_p)[nrow(AC_p)]<-paste0(filename,"_",i)
      
      AutoC_row<-AutoC_triad[((n/3)+1):(2*n/3)]
      AC_dw<-rbind(AC_dw, AutoC_row)
      rownames(AC_dw)[nrow(AC_dw)]<-paste0(filename,"_",i)
      
      AutoC_row<-AutoC_triad[((2*n/3)+1):n]
      AC_ac<-rbind(AC_ac, AutoC_row)
      rownames(AC_ac)[nrow(AC_ac)]<-paste0(filename,"_",i)
      pb$tick()
    }
  }
  
  #data frames are finsihed
  colnames(AC_p)<-names(AutoC_row)
  colnames(AC_dw)<-names(AutoC_row)
  colnames(AC_ac)<-names(AutoC_row)
  #now make pdf with the plots for each lag.
  pdf("/Users/bodhiberoukhim/Desktop/ACplot_lag20_ppm200/Plot.pdf",width=12,height=8)
  ACplot(AC_p,AC_dw,AC_ac)
  dev.off()
  write.csv(AC_p,"/Users/bodhiberoukhim/Desktop/ACplot_lag20_ppm200/AC_p.csv")
  write.csv(AC_dw,"/Users/bodhiberoukhim/Desktop/ACplot_lag20_ppm200/AC_dw.csv")
  write.csv(AC_ac,"/Users/bodhiberoukhim/Desktop/ACplot_lag20_ppm200/AC_ac.csv")
 
}
ACplotcomp(filename_list)




#__________BELOW we have functions that iterate getdata() over many files/filelist. Fucntions copile data in pdfs or lists______

#_______Function that compiles all graphs into one pdf______________

plotcompilation<-function(filename_list){
  pb<-progress_bar$new(total=length(filename_list))
  pdf("/Users/bodhiberoukhim/Desktop/Plot.pdf",width=12,height=8)
  for (filename in filename_list){#iterates the function for every spectra in the big folder
    print(getdata(filename,"Plot"))
    pb$tick()
  }
  dev.off()
}
plotcompilation(filename_list)



#___________Function that compiles all p-values for (all spectra) for the t test and wilcoxon test into one data frame_________
mean_list<-data.frame(Wilcoxon_Pvalue=numeric(0),Ttest_Pvalue=numeric(0))
meanscompilation<-function(filename_list){
  for (filename in filename_list){#iterates the function for every spectra in the big folder
    mean_row<-as.data.frame(getdata(filename,"Means"))  #finds means and changes them to data frame
    mean_list<-rbind(mean_list,mean_row)
    rownames(mean_list)[nrow(mean_list)]<-filename#make row name filename
    print(paste("Calculating p-values of means of",filename))
  }
  return(mean_list)
}
#mean_list<-meanscompilation(filename_list)


#__________Function that compiles normality pvalue fore all spectra in one data frame. Also prints out histogram____
normality_list<-data.frame(Column1=numeric(0),Column2=numeric(0))
normalitycompilation<-function(filename_list){
  for (filename in filename_list){
    normality_row<-as.data.frame(getdata(filename,"Normality"))
    normality_list<-rbind(normality_list,normality_row)
    rownames(normality_list)[nrow(normality_list)]<-filename
    print(paste("Calculating normality of", filename))
  }
  return(normality_list)
}
#normality_list<-normalitycompilation(filename_list)
  


#__________Function that compiles Autocorrelation pvalue fore all spectra in one data frame____
pb<-progress_bar$new(total=400)
AutoC_list<-data.frame()
AutoCcompilation<-function(filename_list){
  for (filename in filename_list[sample(seq(1,length(filename_list)),400)]){
    AutoC_triad<-getdata(filename,"Autocorrelation")
    n<-length(AutoC_triad)
    AutoC_row<-AutoC_triad[1:(n/3)]
    AutoC_list<-rbind(AutoC_list, AutoC_row)
    rownames(AutoC_list)[nrow(AutoC_list)]<-paste(filename,"_P")
    
    AutoC_row<-AutoC_triad[((n/3)+1):(2*n/3)]
    AutoC_list<-rbind(AutoC_list, AutoC_row)
    rownames(AutoC_list)[nrow(AutoC_list)]<-paste(filename,"_DW")
    
    AutoC_row<-AutoC_triad[((2*n/3)+1):n]
    AutoC_list<-rbind(AutoC_list, AutoC_row)
    rownames(AutoC_list)[nrow(AutoC_list)]<-paste(filename,"_AC")
    pb$tick()
  }
  
  colnames(AutoC_list)<-names(AutoC_row)
  return(AutoC_list)
}
AutoC_list<-AutoCcompilation(filename_list)


#__________Function that iterates the ACF plot over all spectra and makes pdf_____
ACFcompilation<-function(filename_list){
  acf_df<-data.frame()
  pb<-progress_bar$new(total=length(filename_list))
  pdf("/Users/bodhiberoukhim/Desktop/ACFplot.pdf",width=12,height=8)
  for (filename in filename_list){
    acfrow<-getdata(filename,"ACF")
    acf_df<-rbind(acf_df,acfrow)
    rownames(acf_df)[nrow(acf_df)]<-filename
    pb$tick()
  }
  dev.off()
  return(acf_df)
}
acf_df<-ACFcompilation(filename_list)



peakcompilation<-function(filename_list){
  pb<-progress_bar$new(total=20)
  pdf("/Users/bodhiberoukhim/Desktop/Peaksartifacts.pdf",width=12,height=8)
  for (filename in filename_list[sample(1:7180, 20)]){
    print(getdata(filename,"Peaks")$plot)
    pb$tick()
  }
  dev.off()
}
peakcompilation(filename_list)


pdf("/Users/bodhiberoukhim/Desktop/Plot_Rating4.pdf",width=12,height=8)
for (filename in summary4$name[sample(1:nrow(summary4), 30)]){
  print(getdata(filename,"Plot"))
}
dev.off()

#________Ghosting compilation fucntion. This function interates ArtifactGhosting() over each spectra
#Does two things: (1)makes pdf of all spectra with ghosting artifacts outlined
#(2) returns a data frame with ghoting info for all spectra

ghostcompilationOLD<-function(filename_list){
  pb<-progress_bar$new(total=length(filename_list))
  pdf("/Users/bodhiberoukhim/Desktop/Ghost_JessicaPipeline_*05_1*5_Height,Velocity,Energy.pdf",width=12,height=8)
  ghostcomp_df<-data.frame()
  for (filename in filename_list){
    result<-getdata(filename,"Ghost")
    row<-result$ghost_df
    plotHeight<-result$plotHeight
    plotVelocity<-result$plotVelocity
    plotEnergy<-result$plotEnergy
    
    ghostcomp_df<-rbind(ghostcomp_df,row)
    print(plotHeight)
    print(plotVelocity)
    print(plotEnergy)
    pb$tick()
  }
  dev.off()
  return(ghostcomp_df)
}
ghostcomp_df<-ghostcompilationOLD(filename_list)


#new ghost compilation msystem
pdf("/Users/bodhiberoukhim/Desktop/GhostNew.pdf",width=12,height=8)
ghostcomp_detailed_df<-data.frame()
ghostcomp_summary_df<-data.frame()
for (filename in filename_list){
  result<-getdata(filename,"Ghost")
  row<-result$ghost_df
  plotVelocity<-result$plotVelocity
  rowsummary<-result$rowsummary
  
  
  ghostcomp_detailed_df<-rbind(ghostcomp_detailed_df,row)
  ghostcomp_summary_df<-rbind(ghostcomp_summary_df,rowsummary)
  rownames(ghostcomp_summary_df)[nrow(ghostcomp_summary_df)]<-filename
  print(plotVelocity)
}
dev.off()

write.csv(ghostcomp_summary_df,"/Users/bodhiberoukhim/Desktop/GhostSummary.csv")
write.csv(ghostcomp_detailed_df,"/Users/bodhiberoukhim/Desktop/GhostDetailed.csv")

