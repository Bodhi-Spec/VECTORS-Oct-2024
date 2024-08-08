#__________checks for faze artifact, where eachof the three main peaks (NAA,Cho,Cr) have negative feet.____
#Note this function uses the doublet peaks from artifact_dublet() so this funciton must be called within artifact_doublet()
artifact_negative_faze<-function(spectra_data,met_data,misc_data,filename,NAAdouble,Chodouble,Crdouble,NAAdf_max,Chodf_max,Crdf_max,plot){
  #criteria for negative faze:
  #(1) valley only needs span=3
  #(2) valley within ppm range (same as for doublet peak)
  #(3) valley needs to be outside of any doublet peak (because if else would detect the valley between met peak and doublet peak as a negative faze.)
  #(4) Vallye has to occur in all three metabilites, in other words be same distance from the met peak (similarly to criteria 5 of doublet)

  #Criteria 1
  vector_valley<-peaks(-spectra_data$RelativeAmp,span=3)#returns vector of length processed data, with only TRUE or FALSE elements. TRUE being that element is a peak
  index_valley<-which(vector_valley==TRUE) #index of valleys
  valley_df<-data.frame(index= index_valley, ppm=spectra_data$ppm[index_valley], amp=spectra_data$RelativeAmp[index_valley],base=spectra_data$RelativeBase[index_valley])#data frame with amplitude and ppm of all valleys and the row index of these data points in spectra_data
  
  # #Criteria 2 (with the relative met peaks)
  # FWHM<-misc_data$FWHM..ppm.
  # FWHM<-FWHM+.05#often times the FWHm falls short by like .002 of the cuttoff so make it slightly bigger
  # FWHM1<-FWHM#create FWHM1, which makes sure FWHM does not exceed .1 in the cuttoff between Cho and Cr
  # if (FWHM>((Chodf_max$ppm-Crdf_max$ppm)/2)){#to prevent overlap betwee met ranges of Cho and Cr
  #   FWHM1=((Chodf_max$ppm-Crdf_max$ppm)/2)
  # }
  # range<-c(NAAdf_max$ppm-FWHM,NAAdf_max$ppm+FWHM,(Crdf_max$ppm-(FWHM-.02)),(Crdf_max$ppm+FWHM1),(Chodf_max$ppm-FWHM1),(Chodf_max$ppm+FWHM-.05))#DELETE/MOVETODOUBLET
  # NAAdfvalley<-subset(valley_df,(ppm>(NAAdf_max$ppm-FWHM)  & ppm<(NAAdf_max$ppm+FWHM)))#dataframe of unfiltere dublet NAA valleys
  # Crdfvalley<-subset(valley_df,(ppm>=(Crdf_max$ppm-(FWHM-.02)) & ppm<=(Crdf_max$ppm+FWHM1)))#dataframe of unfilterd dublet Creatine valleys Note use FWHM1 for boundary between Cho and Cr. Note use higher lower bound
  # Chodfvalley<-subset(valley_df,(ppm>=(Chodf_max$ppm-FWHM1) & ppm<=(Chodf_max$ppm+FWHM-.05)))#dataframe of unfiltered dublet CHoline valleys Note use FWHM1 for boundary between Cho and Cr. Note use lower high bound.
  # 
  
  #Criteria 2 (with literature met peaks)
  FWHM<-misc_data$FWHM..ppm.
  FWHM1<-FWHM+.05
  if (FWHM1>=(3.210943-3.0270-.05)){
    FWHM1<-3.210943-3.0270-.05 #this .05 is so that Cho doesnt incljude an offset (bY .05ppm) Cr peal
  }
  
  NAAdfvalley<-subset(valley_df,(ppm>(2.0080-FWHM1)  & ppm<(2.0080+FWHM1)))#dataframe of unfiltere dublet NAA valleys
  Crdfvalley<-subset(valley_df,(ppm>=(3.0270-FWHM1) & ppm<=(3.0270+FWHM1)))#dataframe of unfilterd dublet Creatine valleys
  Chodfvalley<-subset(valley_df,(ppm>=(3.210944-FWHM1) & ppm<=(3.210944+FWHM1)))#dataframe of unfiltered dublet CHoline valleys

  
  #Criteria 3
  NAAppm<-NAAdf_max$ppm
  Choppm<-Chodf_max$ppm
  Crppm<-Crdf_max$ppm
  NAAdouble$difference<-NAAdouble$ppm-NAAppm#difference in ppm between each double and the met peak
  Chodouble$difference<-Chodouble$ppm-Choppm
  Crdouble$difference<-Crdouble$ppm-Crppm
  #we only want to consider the doublets that are fartherst from met peak, so have max(abs(difference))
  NAAdouble_pos<-subset(NAAdouble,NAAdouble$difference>0)#df with all dublet peaks to the right of the metpeak
  NAAdouble_neg<-subset(NAAdouble,NAAdouble$difference<0)#df with all dublet peaks to the left of the metpeak
  Chodouble_pos<-subset(Chodouble,Chodouble$difference>0)
  Chodouble_neg<-subset(Chodouble,Chodouble$difference<0)
  Crdouble_pos<-subset(Crdouble,Crdouble$difference>0)
  Crdouble_neg<-subset(Crdouble,Crdouble$difference<0)

  

  
  #Now, for each metabolite, there are four sceneriaos 
  #(1) when no doublet peaks length(double_pos,double_neg=0) 
  #(2) when dublet peak(s) are greater than met length(doublet_neg=0, doublet_pos>=0) 
  #(3) when doublet peak(s) are less than metpeak (length(doublet_neg>0, doublet_pos=0)) 
  #(4) when doublet peak(s) exist both greater and less than metpeak (length(doublet_neg,doublet_pos>0))
  #However we can simplify these cases, and just say that if length(dublet_pos,doublet_neg=0), then (doublet_neg,doublet_pos)$difference=0
  #Adn then we simply remove the range between the doublet min or max from the ppm of vallyes we are dtecting
  if (nrow(NAAdouble_pos)==0){
    NAAdouble_max<-data.frame(index=NA,ppm=NA,amp=NA,base=NA,difference=0)
  }else if (nrow(NAAdouble_pos>0)){
    NAAdouble_max<-subset(NAAdouble_pos,NAAdouble_pos$difference==max(NAAdouble_pos$difference))#this is the row containing the dublet peak farthest from the metpeak (closest to the outskirts) in the right side (greater than metpeak)
  }
  if (nrow(NAAdouble_neg)==0){
    NAAdouble_min<-data.frame(index=NA,ppm=NA,amp=NA,base=NA,difference=0)
  }else if (nrow(NAAdouble_neg)>0){
    NAAdouble_min<-subset(NAAdouble_neg,NAAdouble_neg$difference==min(NAAdouble_neg$difference))#this is the row contianing the dublet peak farthest from the metpeak on the left side (less that metpeak)
  }
  
  if (nrow(Chodouble_pos)==0){
    Chodouble_max<-data.frame(index=NA,ppm=NA,amp=NA,base=NA,difference=0)
  }else if (nrow(Chodouble_pos)>0){
    Chodouble_max<-subset(Chodouble_pos,Chodouble_pos$difference==max(Chodouble_pos$difference))
  }
  if (nrow(Chodouble_neg)==0){
    Chodouble_min<-data.frame(index=NA,ppm=NA,amp=NA,base=NA,difference=0)
  }else if (nrow(Chodouble_neg)>0){
    Chodouble_min<-subset(Chodouble_neg,Chodouble_neg$difference==min(Chodouble_neg$difference))
  }
  
  if (nrow(Crdouble_pos)==0){
    Crdouble_max<-data.frame(index=NA,ppm=NA,amp=NA,base=NA,difference=0)
  }else if (nrow(Crdouble_pos)>0){
    Crdouble_max<-subset(Crdouble_pos,Crdouble_pos$difference==max(Crdouble_pos$difference))
  }
  if (nrow(Crdouble_neg)==0){
    Crdouble_min<-data.frame(index=NA,ppm=NA,amp=NA,base=NA,difference=0)
  }else if (nrow(Crdouble_neg)>0){
    Crdouble_min<-subset(Crdouble_neg,Crdouble_neg$difference==min(Crdouble_neg$difference))
  }
  
  #Remove from dfvalley the range within dublets
  NAAdfvalley<-subset(NAAdfvalley,!(ppm<(NAAdouble_max$difference+NAAppm)  & ppm>(NAAdouble_min$difference+NAAppm)))#dataframe of unfiltere dublet NAA valleys
  Chodfvalley<-subset(Chodfvalley,!(ppm<(Chodouble_max$difference+Choppm)  & ppm>(Chodouble_min$difference+Choppm)))#dataframe of unfiltere dublet Cho valleys
  Crdfvalley<-subset(Crdfvalley,!(ppm<(Crdouble_max$difference+Crppm)  & ppm>(Crdouble_min$difference+Crppm)))#dataframe of unfiltere dublet Cr valleys
  
  #Criteria 4
  NAAindex<-NAAdf_max$index
  Choindex<-Chodf_max$index
  Crindex<-Crdf_max$index
  NAAdifference<-NAAdfvalley$index-NAAindex #difference is defined by i$index-max$index. This list holds the difference in ideces from a given NAA negative valley to the NAA peak
  Chodifference<-Chodfvalley$index-Choindex
  Crdifference<-Crdfvalley$index-Crindex
  
  
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
  
  #now, we have the differneces for each metabolite that are replicated across all metbaolites. To get indexes as following. since maxpeak-index=difference, maxpeak-difference=index
  NAAselectindex<-NAAindex+NAAselectdiff
  Choselectindex<-Choindex+Choselectdiff
  Crselectindex<-Crindex+Crselectdiff
  
  #now we need to make new data frame for the verticle lines that are duplicates in all metbaolites. we have the indexes from NAAdf, Chodf, and Crdf
  NAAvalley<-subset(NAAdfvalley,index %in% NAAselectindex)
  Chovalley<-subset(Chodfvalley,index %in% Choselectindex)
  Crvalley<-subset(Crdfvalley,index %in% Crselectindex)
  
  
  
  criteria123<-c(NAAdfvalley$ppm,Crdfvalley$ppm,Chodfvalley$ppm)
  faze_artifact<-c(NAAvalley$ppm,Crvalley$ppm,Chovalley$ppm)#ppm of faze artifact peaks
 
  #plot<-plot+geom_vline(xintercept =criteria123, linetype = "dashed", color ="brown")#brown is criteria 1+2+3
  plot<-plot+geom_vline(xintercept =faze_artifact, linetype = "dashed", color ="blue")#blue is all criteria
  if (length(faze_artifact)>0){
    existfaze<-TRUE
  }else{
    existfaze<-FALSE
  }
  plot<-plot+ geom_label_npc(label=paste0("exist faze= ",existfaze),npcx=.93,npcy=.8)
  return(list(plot=plot,existfaze=existfaze))
}


#_____Chekcs if any part oif the graph is negative______
#Further, specifies if this negative is oversupression artifact OR something else
artifact_negative_graph<-function(spectra_data,met_data,misc_data,filename){
  negative_ppm<-spectra_data$ppm[spectra_data$ProcessedData<=0]
  negative_ppm<-negative_ppm[negative_ppm>.35] #somtime the very left side of the spectra dips below 0 slightly
    #For skyler's QA she only cares about negatives to the left of NAAish
    #negative_ppm<-negative_ppm[negative_ppm>1.8]#this
  
  oversupressionnegative<-negative_ppm[negative_ppm>3.4]#alter the cuttogg
  if (length(negative_ppm)>0){
    anynegative<-TRUE
  }else{
    anynegative<-FALSE
  }
  if (length(oversupressionnegative)>0){
    oversupressionnegative<-TRUE
  }else{
    oversupressionnegative<-FALSE
  }
  return(list(anynegative=anynegative,oversupressionnegative=oversupressionnegative))
}


#____Checls if any pafrt of the the graph is below the baseline
artifact_negative_baseline<-function(spectra_data,met_data,misc_data,filename){
  spectra_data<-subset(spectra_data,ppm>.35) #same as above remove very left side of graph
  belowbaseline<-any(spectra_data$ProcessedData<=spectra_data$Background)
  return(list(belowbaseline=belowbaseline))
}

PlotNegative<-function(spectra_data,met_data,misc_data,filename){
  belowBaseline<-artifact_negative_baseline(spectra_data,met_data,misc_data,filename)$belowbaseline
  anyNegative<-artifact_negative_graph(spectra_data,met_data,misc_data,filename)$anynegative
  
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=ProcessedData),size=.5,color="red")+geom_line(aes(y=RawData),size=.02,color="black")+geom_line(aes(y=Background),color="blue")#make line plot of the raw data (black),proccessed data (red),background(blue)
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+geom_hline(yintercept = 0, color = "black", linetype = "solid")#add line in xaxis
  plot<-plot+ggtitle(filename)+ylab("Amplitude")+geom_label_npc(label=paste0("anyNegative: ", anyNegative),npcx=.93,npcy=.95)+geom_label_npc(label=paste0("belowBaseline: ", belowBaseline),npcx=.93,npcy=.87)#Add labels
  return(plot)
}
