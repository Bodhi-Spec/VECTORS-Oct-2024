# Contians 3 functions.
# Glxpeaks() calculates the quality of Glx and its gamma peak
# mIpeaks() calculates the quality of the mI peak
# plotgraphGlxmI() plots the spectrum and its 5 different scores



#________________Function that plots the graph given a single filename________________
plotgraphGlxmI<-function(spectra_data,met_data,misc_data,filename,GlxMerge,GlxDistinct,mIMerge,mIDistinct){
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=ProcessedData),size=.5,color="red")+geom_line(aes(y=RawData),size=.05,color="black")+geom_line(aes(y=Background),color="blue")#graphs RELATIVE data (procesed+raw) + base
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Amplitude")#Add labels
  plot<-plot+geom_vline(xintercept=c(2.245,2.385,3.52,3.68), linetype = "dashed", color = "black")#plot ranges for GLx/mI
  
  plot<-plot+ geom_label_npc(label=paste0("Glx Merge = ",GlxMerge),npcx=.93,npcy=.95)
  plot<-plot+ geom_label_npc(label=paste0("Glx Distinct = ",GlxDistinct),npcx=.93,npcy=.87)
  plot<-plot+ geom_label_npc(label=paste0("mI Merge = ",mIMerge),npcx=.93,npcy=.79)
  plot<-plot+ geom_label_npc(label=paste0("mI Distinct = ",mIDistinct),npcx=.93,npcy=.71)
  return(plot)
}

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
  
  return(list(Merge=Merge,Distinct=Distinct)) #Returns the Merge and Distinct scores
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

  return(list(Merge=Merge,Distinct=Distinct)) #Returns the Merge and Distinct scores
}
