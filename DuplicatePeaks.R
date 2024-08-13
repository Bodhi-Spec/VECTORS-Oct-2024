### Author: Bodhi Beroukhim, bberoukhim@bwh.harvard.edu
### Contact: Alexander P. Lin, aplin@bwh.harvard.edu
### Date last updated: August 12, 2024

### License: GNU GPL2, Copyright (C) 2024 Brigham and Women's Hospital



#This file contains the fucntion (called DuplicatePeaks()) that detects movement artifacts that present as duplicate peaks
#DuplicatePeaks() is called through the GetData() function which provides it with the data in the parameters. To call these functions, use the command getdata(filename,"DuplicatePeaks")
#Outputs the spectra outlining any duplicate peaks and a TRUE/FALSE statement of if a movement artifact was detected

DuplicatePeaks<-function(spectra_data,met_data,misc_data,filename){
  #Prepare the plot of the spectrum
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=ProcessedData),size=.5,color="red")+geom_line(aes(y=RawData),size=.05,color="black")+geom_line(aes(y=Background),color="blue") #graphs processed, raw, and baseline data
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal() #reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Amplitude")#Add labels
  
  #Criteria for what constitutes a doublet peak/movement artifact (only checking doublet peaks of NAA, CHo Cr). The peak has to:
    #(1) have a span > 3. 
    #(2) be in the ppm range of metabolite (exact range elucidated below)

  #Criteria A (orange vertical lines). 
    #(3) If a criteria (1+2) is found in all three metabolites all, then doublet peaks must: Difference in ppm from a doublet peak to the main met peak (+/- ppm diffrence of two adjacent data points) is found in all three metabolitees
  
  #Criteria B (green vertical lines): Note only apply this criteria to Cho and Cr because withina NAA often falsely detects NAAG as a dublet peak
    #(4) For Cho, Cr, peaks need to be above 1/3(Amp-Baseline)+Baseline.  
    #(5) The minimium in the valley between adjacent peaks must be .01 less in amplitude than the smaller of the two peaks for Cho and Cr 
     
  

  #In other words, criteria (1+2) always have to be met, but either (3 -> Green) OR (4+5 -> Orange) have to be met to be a double


   
  
  #_____Finds all peaks that fulfill Criteria 1+2 __________
  #Note met range is (metabolite exact location +/- (FWHM+.05). But between Cho and Cr FWHM has max of (Cho-Cr)/2 But on upper bound of Cho and lower bound of Cr met FWHM is FWHM!=FWHM+.05) 
  FWHM<-misc_data$FWHM..ppm.
  FWHM1<-FWHM+.05
  if (FWHM1>=(3.210943-3.0270-.05)){
    FWHM1<-3.210943-3.0270-.05 #this .05 is so that Cho doesnt incljude an offset (bY .05ppm) Cr peal
  }

  #note the peaks() function is from the splus2R package
  vector_peak<-peaks(spectra_data$RelativeAmp,span=3)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak
  index_peaks<-which(vector_peak==TRUE) #index of peaks
  peaks_df<-data.frame(index= index_peaks, ppm=spectra_data$ppm[index_peaks], amp=spectra_data$RelativeAmp[index_peaks],base=spectra_data$RelativeBase[index_peaks])#data frame with amplitude and ppm of allpeaks and the row index of these data points in spectra_data
  range<-c(2.0080-FWHM1,2.0080+FWHM1,3.0270-FWHM1,3.0270+FWHM1,3.210944-FWHM1,3.210944+FWHM1)#DELETE/MOVETODOUBLET
  NAAdfraw<-subset(peaks_df,(ppm>(2.0080-FWHM1)  & ppm<(2.0080+FWHM1)))#dataframe of unfiltere dublet NAA peaks.
  Crdfraw<-subset(peaks_df,(ppm>=(3.0270-FWHM1) & ppm<=(3.0270+FWHM1)))#dataframe of unfilterd dublet Creatine peaks. Note use FWHM1 for boundary between Cho and Cr. Note use higher lower bound
  Chodfraw<-subset(peaks_df,(ppm>=(3.210944-FWHM1) & ppm<=(3.210944+FWHM1)))#dataframe of unfiltered dublet CHoline peaks. Note use FWHM1 for boundary between Cho and Cr. Note use lower high bound.
  
  
  #ERROR MEASURES ---> will report an error if no peaks are found in a metabolite range, implying that the spectrum is super extraneous
  ploterror<-plot+ geom_label_npc(label=paste0("exist faze= error"),npcx=.93,npcy=.95)
  ploterror<-ploterror+ geom_label_npc(label=paste0("exist duplicate= error"),npcx=.93,npcy=1)
  if (nrow(NAAdfraw)==0||nrow(Crdfraw)==0||nrow(Chodfraw)==0){
    error=TRUE
    return(list(plot=ploterror, existfaze="error",existduplicate="error")) #returns empty plot
  }

  
  NAAdf<-NAAdfraw
  Crdf<-Crdfraw
  Chodf<-Chodfraw

  #___Now find the specific metabolite mainpeak
  NAApeakdf<-subset(NAAdf,ppm>(2.0080-.05)  & ppm<(2.0080+.05))#Smaller range to find NAApeak
  Crpeakdf<-subset(Crdf,ppm>(3.0270-.05)  & ppm<(3.0270+.05))#smaller range to find Crpeak
  Chopeakdf<-subset(Chodf,ppm>(3.210944-.05)  & ppm<(3.210944+.05))#Smaller range to find Chopeak 
  
  #ERROR MEASURES pt2 ---> if no main metabolite peak is detected
  if (nrow(NAApeakdf)==0||nrow(Crpeakdf)==0||nrow(Chopeakdf)==0){
    error=TRUE
    return(list(plot=ploterror, existfaze="error",existduplicate="error")) #returns empty plot
  }
  
  
  #_____________Criteria 3_______
  NAAmaxindex<-NAApeakdf$index[which.max(NAApeakdf$amp)]
  NAAmax<-which(NAAdf$index==NAAmaxindex)
  NAAdifference<-c() #difference is defined by index[max]-index[i]. This list holds the difference in ideces from a given NAA dublet peak to the NAA peak
  for (i in 1:NAAmax-1){
    NAAdifference<-c(NAAdifference, NAAmaxindex-NAAdf$index[i])
  }
  for (i in (NAAmax+1):nrow(NAAdf)){
    NAAdifference<-c(NAAdifference, NAAmaxindex-NAAdf$index[i])
  }
  
  Crmaxindex<-Crpeakdf$index[which.max(Crpeakdf$amp)]
  Crmax<-which(Crdf$index==Crmaxindex)
  Crdifference<-c() #difference is defined by index[max]-index[i]
  for (i in 1:Crmax-1){
    Crdifference<-c(Crdifference, Crmaxindex-Crdf$index[i])
  }
  for (i in (Crmax+1):nrow(Crdf)){
    Crdifference<-c(Crdifference, Crmaxindex-Crdf$index[i])
  }
  
  Chomaxindex<-Chopeakdf$index[which.max(Chopeakdf$amp)]
  Chomax<-which(Chodf$index==Chomaxindex)
  Chodifference<-c() #difference is defined by index[max]-index[i]
  for (i in 1:Chomax-1){
    Chodifference<-c(Chodifference, Chomaxindex-Chodf$index[i])
  }
  for (i in (Chomax+1):nrow(Chodf)){
    Chodifference<-c(Chodifference, Chomaxindex-Chodf$index[i])
  }
  
  #Get rid of main peak itself
  NAAdifference<-NAAdifference[NAAdifference!=0]
  Crdifference<-Crdifference[Crdifference!=0]
  Chodifference<-Chodifference[Chodifference!=0]
  
  
  #now need to find if there is a value +/-2 that is found in all three lists. 
  rangeNAA<-c()#All possible values that will indicate that NAAdifference element is overlap
  for (element in NAAdifference){
    rangeNAA<-c(rangeNAA,element-1,element,element+1) #+/- 1 because that menas maximum of +/-2 difference (ie. if NAA is +1 and Cho is -1 and they mathc)
  }
  rangeCho<-c()
  for (element in Chodifference){
    rangeCho<-c(rangeCho,element-1,element,element+1)
  }
  rangeCr<-c()
  for (element in Crdifference){
    rangeCr<-c(rangeCr,element-1,element,element+1)
  }
  
  common<- intersect(intersect(rangeNAA,rangeCho),rangeCr) #see what differences/indexes are contained in all three lists
  NAAselectdiff<-c(intersect(NAAdifference,common),intersect(NAAdifference-1,common)+1,intersect(NAAdifference+1,common)-1)#NAA selected elements, which are the differences which are in all metbaolites
  Choselectdiff<-c(intersect(Chodifference,common),intersect(Chodifference-1,common)+1,intersect(Chodifference+1,common)-1)
  Crselectdiff<-c(intersect(Crdifference,common),intersect(Crdifference-1,common)+1,intersect(Crdifference+1,common)-1)
  NAAselectdiff<-unique(NAAselectdiff)
  Choselectdiff<-unique(Choselectdiff)
  Crselectdiff<-unique(Crselectdiff)
  #sometimes we assign the three main peaks to be duplicates under this critera 5. So remove these points, which have a element=0
  NAAselectdiff<-NAAselectdiff[NAAselectdiff!=0]
  Choselectdiff<-Choselectdiff[Choselectdiff!=0]
  Crselectdiff<-Crselectdiff[Crselectdiff!=0]
  
  #now, we have the differneces for each metabolite that are replicated across all metbaolites. To get indexes as following. since maxpeak-index=difference, maxpeak-difference=index
  NAAselectindex<-NAAdf$index[NAAmax]-NAAselectdiff
  Choselectindex<-Chodf$index[Chomax]-Choselectdiff
  Crselectindex<-Crdf$index[Crmax]-Crselectdiff
  
  #now we need to make new data frame for the verticle lines that are duplicates in all metbaolites. we have the indexes from NAAdf, Chodf, and Crdf
  NAAcriteria3df<-subset(NAAdf,index %in% NAAselectindex)
  Chocriteria3df<-subset(Chodf,index %in% Choselectindex)
  Crcriteria3df<-subset(Crdf,index %in% Crselectindex)
  
  



  
  #___________________Criteria 4_-___
  #NAAabove<- max(NAAdf$amp)/3+2*(NAAdf$base[which.max(NAAdf$amp)])/3
  Crabove<- max(Crdf$amp)/3+2*(Crdf$base[which.max(Crdf$amp)])/3
  Choabove<- max(Chodf$amp)/3+2*(Chodf$base[which.max(Chodf$amp)])/3
  
  #NAAdf<-subset(NAAdf,amp> NAAabove)
  Crdf<-subset(Crdf,amp> Crabove)
  Chodf<-subset(Chodf,amp> Choabove)
  
  
  

  
  #this for loop fulfilles crieteria (5). it finds the valley between two adjacent peaks and makes sure its significantly less than both peaks
  #Note NAA is not checkde for this critera
  #For Cr
  offset<-0
  if (nrow(Crdf)>1){
    for (i in 1:(nrow(Crdf)-1)){
      start<-Crdf$index[i-offset]
      end<-Crdf$index[i+1-offset]
      valley<-min(spectra_data$RelativeAmp[start:end])
      if ((valley+.01) > min(Crdf$amp[i-offset],Crdf$amp[i+1-offset])){ #if valley is within .01 of lowest peak, dont count this peak
        Crdf<-Crdf[-which(Crdf$amp==min(Crdf$amp[i-offset],Crdf$amp[i+1-offset])),]
        offset<-offset+1
      }
    }
  }else if (nrow(Crdf)==0){
    print("ERROR, no Cr detected")
  }
  #For Cho
  offset<-0
  if (nrow(Chodf)>1){
    for (i in 1:(nrow(Chodf)-1)){
      start<-Chodf$index[i-offset]
      end<-Chodf$index[i+1-offset]
      valley<-min(spectra_data$RelativeAmp[start:end])
      if ((valley+.01) > min(Chodf$amp[i-offset],Chodf$amp[i+1-offset])){ #if valley is within .01 of lowest peak, dont count this peak. The spectra will have some movement artifact, just not substantial
        Chodf<-Chodf[-which(Chodf$amp==min(Chodf$amp[i-offset],Chodf$amp[i+1-offset])),]
        offset<-offset+1
      }
    }
  }else if (nrow(Chodf)==0){
    print("ERROR, no Cho detected")
  }

  #find maxes (for later). 
  NAAdf_max<-subset(NAApeakdf,NAApeakdf$amp==max(NAApeakdf$amp))
  Chodf_max<-subset(Chopeakdf,Chopeakdf$amp==max(Chopeakdf$amp))
  Crdf_max<-subset(Crpeakdf,Crpeakdf$amp==max(Crpeakdf$amp))
  
  #Now need to remove the actual met peak from NAAdf, Chodf, Crdf
  #NAAdf<-subset(NAAdf,NAAdf$index!=NAAdf_max$index)
  Chodf<-subset(Chodf,Chodf$index!=Chodf_max$index)
  Crdf<-subset(Crdf,Crdf$index!=Crdf_max$index)

  criteria3peaks<-c(NAAcriteria3df$ppm,Crcriteria3df$ppm,Chocriteria3df$ppm) #ppm values of the peaks that show up in all metabolites (criteria 1,2,5)
  filtered_x_peaks<-c(Crdf$ppm,Chodf$ppm) #ppm of pekas that survive crtieria 1,2,4,5
  
  
  #add verticle lines identifying peaks on the spectrum
  plot<-plot+geom_vline(xintercept =filtered_x_peaks, linetype = "dashed", color = "green") #Green lines are for criteria 1,2,4,5
  plot<-plot+geom_vline(xintercept =criteria3peaks, linetype = "dashed", color = "orange") #orange lines are for criteria 1,2,3
  plot<-plot+geom_vline(xintercept =range, linetype = "dashed", color = "black") #black vertical lines denote the range in which we are scanning for doublet peaks
  if ((length(filtered_x_peaks)>0)||(length(criteria3peaks)>0)){ #if any peaks are identifiied
    existduplicate<-TRUE
  }else{
    existduplicate<-FALSE
  }
  plot<-plot+ geom_label_npc(label=paste0("exist duplicate= ",existduplicate),npcx=.93,npcy=.9) #add label of if movment artifact on the plot
  return(list(plot=plot,existDuplicate=existduplicate))
}

