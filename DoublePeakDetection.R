# Contians 3 functions.
# Glxpeaks() calculates the quality of Glx and its gamma peak
# Inspeaks() calculates the quality of the Ins peak
# plotgraphGlxIns() plots the spectrum and its 5 different scores



#________________Function that plots the graph given a single filename________________
plotgraphGlxIns<-function(spectra_data,met_data,misc_data,filename,Glxscore,Glxextrascore1,Glxgammascore,Insscore,Insextrascore1){
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=ProcessedData),size=.5,color="red")+geom_line(aes(y=RawData),size=.05,color="black")+geom_line(aes(y=Background),color="blue")#graphs RELATIVE data (procesed+raw) + base
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Amplitude")#Add labels
  plot<-plot+geom_vline(xintercept=c(2.245,2.385,3.52,3.68), linetype = "dashed", color = "black")#plot ranges for GLx/Ins
  
  plot<-plot+ geom_label_npc(label=paste0("Glx Merge = ",Glxscore),npcx=.93,npcy=.95)
  plot<-plot+ geom_label_npc(label=paste0("Glx Distinct = ",Glxextrascore1),npcx=.93,npcy=.87)
  plot<-plot+ geom_label_npc(label=paste0("Gamma = ",Glxgammascore),npcx=.93,npcy=.79)
  plot<-plot+ geom_label_npc(label=paste0("mI Merge = ",Insscore),npcx=.93,npcy=.71)
  plot<-plot+ geom_label_npc(label=paste0("mI Distinct = ",Insextrascore1),npcx=.93,npcy=.63)
  return(plot)
}

Glxpeaks<-function(spectra_data){
  bigdata<-subset(spectra_data,ppm>2.1 & ppm<2.5) #this is the data frame where I'm going to search for the peaks of 3<span<.1. These numbers will change
  #filename1
  row<-c() #will be the row containing the number of peaks for all the spans
  score=1#the final score should be the max span where tehre are still two (or more) peaks
  extrascore1=1#final score should maximum span that #peaks>2 exist. so 1 is optimal
  for (span in 2*1:18+1){
    vector_peak<-peaks(bigdata$RelativeAmp,span)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak
    index_peaks<-which(vector_peak==TRUE) #index of peaks
    peaks_df<-data.frame(index= index_peaks, ppm=bigdata$ppm[index_peaks], amp=bigdata$RelativeAmp[index_peaks],base=bigdata$RelativeBase[index_peaks])#data frame with amplitude and ppm of allpeaks and the row index of these data points in spectra_data
    Glxraw<-subset(peaks_df,(ppm>2.245  & ppm<2.385))#subset of the peaks within the GLX tentative range
    numpeaks<-nrow(Glxraw)
    row<-c(row,numpeaks)#add the number of peaks for the span in this given iteration
    if (numpeaks>1){
      score=span 
    }
    if (numpeaks>2){
      extrascore1=span
    }
  }
  
  compositescore1<-score/extrascore1
  
  if (extrascore1>1){
    compositescore2=1
  }else{
    compositescore2=score
  }
  
  gamma_data<-subset(spectra_data,ppm>3.65 & ppm<3.85)
  gammascore=1
  #now check for gamma peak:
  for (span in 2*1:18+1){
    vector_peak<-peaks(gamma_data$RelativeAmp,span)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak
    index_peaks<-which(vector_peak==TRUE) #index of peaks
    peaks_df<-data.frame(index= index_peaks, ppm=gamma_data$ppm[index_peaks], amp=gamma_data$RelativeAmp[index_peaks],base=gamma_data$RelativeBase[index_peaks])#data frame with amplitude and ppm of allpeaks and the row index of these data points in spectra_data
    Gammaraw<-subset(peaks_df,(ppm>3.7  & ppm<3.8))#subset of the peaks within the GLX tentative range
    numpeaks<-nrow(Gammaraw)
    if (numpeaks==1){
      gammascore=span 
    }
  }
  return(list(score=score,extrascore1=extrascore1,gammascore=gammascore,compositescore1=compositescore1,compositescore2=compositescore2))#retunrs the entire list of spans and #pekas and also the score
}

Inspeaks<-function(spectra_data){
  bigdata<-subset(spectra_data,ppm>3.4 & ppm<3.8) #this is the data frame where I'm going to search for the peaks of 3<span<.1. These numbers will change
  #filename1
  row<-c() #will be the row containing the number of peaks for all the spans
  score=1#the final score should be the max span where tehre are still two (or more) peaks
  extrascore1=1#final score should maximum span that #peaks>2 exist. so 0 is optimal
  for (span in 2*1:28+1){
    vector_peak<-peaks(bigdata$RelativeAmp,span)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak
    index_peaks<-which(vector_peak==TRUE) #index of peaks
    peaks_df<-data.frame(index= index_peaks, ppm=bigdata$ppm[index_peaks], amp=bigdata$RelativeAmp[index_peaks],base=bigdata$RelativeBase[index_peaks])#data frame with amplitude and ppm of allpeaks and the row index of these data points in spectra_data
    Insraw<-subset(peaks_df,(ppm>3.52  & ppm<3.68))#subset of the peaks within the Ins tentative range
    numpeaks<-nrow(Insraw)
    row<-c(row,numpeaks)#add the number of peaks for the span in this given iteration
    if (numpeaks>1){
      score=span 
    }
    if (numpeaks>2){
      extrascore1=span
    }
  }
  compositescore1<-score/extrascore1
  
  if (extrascore1>1){
    compositescore2=1
  }else{
    compositescore2=score
  }
  return(list(score=score,extrascore1=extrascore1,compositescore1=compositescore1))#retunrs the entire list of spans and #pekas and also the score
}
