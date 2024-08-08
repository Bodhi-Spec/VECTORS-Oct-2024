#This file supports the Vertical Shift tests --> the functions (anyNegative() and belowBaseline()) that detects if the processed data is negative and below the baseline respectively
#Both anyNegative() and belowBaseline() is called through the GetData() function which provides it with the data in the parameters

#
anyNegative<-function(spectra_data,met_data,misc_data,filename){
  negative_ppm<-spectra_data$ppm[spectra_data$ProcessedData<=0] #subset all the points which are negative
  negative_ppm<-negative_ppm[negative_ppm>1.8] #subsets all the points which are in the metabolite ppm range
  if (length(negative_ppm)>0){
    anynegative<-TRUE
  }else{
    anynegative<-FALSE
  }
  return(list(anynegative=anynegative))
}


#____Checls if any pafrt of the the graph is below the baseline
belowBaseline<-function(spectra_data,met_data,misc_data,filename){
  spectra_data<-subset(spectra_data,ppm>1.8) #only detects points in the metbaolite ppm range
  belowbaseline<-any(spectra_data$ProcessedData<=spectra_data$Background) #checks if any of the processed data dips below the baseline
  return(list(belowbaseline=belowbaseline))
}
