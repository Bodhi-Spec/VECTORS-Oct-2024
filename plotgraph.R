#________________Function that plots the spectrum given a single filename________________
#Includes the SNR and FWHM on the graph.
#This function is called and given the parameters from the getdata() function
#Note this function is not part of the QA, but is provided to view the spectra


plotgraph<-function(spectra_data,met_data,misc_data,filename){  
  label_data<-data.frame(xcord=c(1,1),ycord=c(.95,1),text=c(paste("SNR: ",misc_data$S.N),paste("FWHM:",misc_data$FWHM..ppm.)))
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=ProcessedData),size=.7,color="red")+geom_line(aes(y=RawData),size=.2,color="black")+geom_line(aes(y=Background),color="blue")+geom_line(aes(y=Residual),size=.3,color="purple")#make line plot of the raw data (black),proccessed data (red),background(blue),residuals(purple)
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Amplitude")+geom_text_npc(data = label_data, aes(npcx = xcord, npcy = ycord, label = text),size=4)#Add labels
  return(plot)
}
