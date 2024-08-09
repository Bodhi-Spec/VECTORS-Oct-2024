#________________Function that plots the graph given a single filename________________
plotgraph<-function(spectra_data,met_data,misc_data,filename){
  model <- lm(Residual~ ppm, data = spectra_data)#linear regression model of residuals
  Slope<-signif(coef(model)["ppm"],digits =3)#Takes slope to 3 digits
  yintercept<-signif(coef(model)["(Intercept)"],digits=3)#takes interceot to 3 digits
  
  label_data<-data.frame(xcord=c(1,1,1),ycord=c(.9,.95,1),text=c(paste("SNR: ",misc_data$S.N),paste("FWHM:",misc_data$FWHM..ppm.),paste("y=",Slope,"x+",yintercept)))
  
  
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=ProcessedData),size=.7,color="red")+geom_line(aes(y=RawData),size=.2,color="black")+geom_line(aes(y=Background),color="blue")+geom_line(aes(y=Residual),size=.3,color="purple")#make line plot of the raw data (black),proccessed data (red),background(blue),residuals(purple)
  plot<-plot+geom_smooth(aes(y=Residual),method = "lm", color="black",se=FALSE)#add regression line of residuals
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Amplitude")+geom_text_npc(data = label_data, aes(npcx = xcord, npcy = ycord, label = text),size=4)#Add labels
  return(plot)
}

#________________Function that plots the metabolites given a single filename________________
plotmet<-function(spectra_data,met_data,misc_data,filename){
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=NAA+NAAG+Cr+PCr+GPC+PCh),size=.7,color="red")+geom_line(aes(y=RawData),size=.2,color="black")+geom_line(aes(y=Background),color="blue")+geom_line(aes(y=Residual),size=.3,color="purple")#make line plot of the raw data (black),proccessed data (red),background(blue),residuals(purple)
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Amplitude")
  return(plot)
}

rangetest<-function(spectra_data,met_data,misc_data,filename){
  FWHM<-misc_data$FWHM..ppm.
  FWHM<-FWHM+.05
  range<-c(2.0080-FWHM,2.0080+FWHM,(3.0270-FWHM),(3.0270+FWHM),(3.210944-FWHM),(3.210944+FWHM))#DELETE/MOVETODOUBLET
  plot<-ggplot(spectra_data,aes(x=ppm))+geom_line(aes(y=RelativeAmp),size=.7,color="red")+geom_line(aes(y=RelativeBase),color="blue")#graphs RELATIVE data + base
  plot<-plot+geom_vline(xintercept =range, linetype = "dashed", color = "black")#range in which we are scanning both for faze and dublet
  plot<-plot+scale_x_reverse(breaks = seq(0, 4, .2),minor_breaks = seq(0, 4, 0.1))+theme_minimal()#reverse x axis (like in LC model)+axis tick mars
  plot<-plot+ggtitle(filename)+ylab("Relative Amplitude")#Add labels
  
  return(plot)
}