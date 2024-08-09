#This file supports the Glx and mI Peaks tests --> the functions (Glxpeaks() and mIpeaks()) that detects the Merge and Distinct metrics for both Glx and mI
#Both Glxpeaks() and mI() are called through the GetData() function which provides it with the data in the parameters. To call these functions, use the command getdata(filename,"Glx_mI_Peaks")
#Outputs a list of all four scores 

Glxpeaks<-function(spectra_data){
  bigdata<-subset(spectra_data,ppm>2.1 & ppm<2.5) #this is the data frame where I'm going to search for the peaks 
  row<-c() #will be the row containing the number of peaks for all the spans
  Merge=1  #Merge metric: the final score should be the max span where there are still two (or more) peaks. Large implies better DQC
  Distinct=1 #Distinct metric: final score should maximum span that #peaks>2 exist. Lower implies better DQC
  for (span in 2*1:18+1){
    vector_peak<-peaks(bigdata$RelativeAmp,span)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak. Note peaks() function within Splus2R package
    index_peaks<-which(vector_peak==TRUE) #index of peaks
    peaks_df<-data.frame(index= index_peaks, ppm=bigdata$ppm[index_peaks], amp=bigdata$RelativeAmp[index_peaks],base=bigdata$RelativeBase[index_peaks])#data frame with amplitude and ppm of allpeaks and the row index of these data points in spectra_data
    Glxraw<-subset(peaks_df,(ppm>2.245  & ppm<2.385))#subset of the peaks within the Glx tentative range
    numpeaks<-nrow(Glxraw)
    row<-c(row,numpeaks)#add the number of peaks for the span in this given iteration
    if (numpeaks>1){
      Merge=span 
    }
    if (numpeaks>2){
      Distinct=span
    }
  }
  
  return(list(GlxMerge=Merge,GlxDistinct=Distinct)) #Returns the Merge and Distinct scores
}

mIpeaks<-function(spectra_data){
  bigdata<-subset(spectra_data,ppm>3.4 & ppm<3.8) #this is the data frame where I'm going to search for the peaks 
  row<-c() #will be the row containing the number of peaks for all the spans
  Merge=1  #Merge metric: the final score should be the max span where there are still two (or more) peaks. Large implies better DQC
  Distinct=1 #Distinct metric: final score should maximum span that #peaks>2 exist. Lower implies better DQC
  
  for (span in 2*1:28+1){
    vector_peak<-peaks(bigdata$RelativeAmp,span)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak. Note peaks() function within Splus2R package
    index_peaks<-which(vector_peak==TRUE) #index of peaks
    peaks_df<-data.frame(index= index_peaks, ppm=bigdata$ppm[index_peaks], amp=bigdata$RelativeAmp[index_peaks],base=bigdata$RelativeBase[index_peaks])#data frame with amplitude and ppm of allpeaks and the row index of these data points in spectra_data
    mIraw<-subset(peaks_df,(ppm>3.52  & ppm<3.68))#subset of the peaks within the mI tentative range
    numpeaks<-nrow(mIraw)
    row<-c(row,numpeaks)#add the number of peaks for the span in this given iteration
    if (numpeaks>1){
      Merge=span 
    }
    if (numpeaks>2){
      Distinct=span
    }
  }

  return(list(mIMerge=Merge,mIDistinct=Distinct)) #Returns the Merge and Distinct scores
}
