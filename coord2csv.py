#This function converts COORD files into csv files
#Run the code. It will prompt in the command line to provide the path to the folder containing the COORD files
#Once provided, the code will create three .csv files that contain 1) the data points 2) the metbolite information (CRLB, SD, etc) and 3) Miscelleneous stats such as SNR and FWHM
#These .csv files will be outputted in the original folder

#______Import packages________________
import sys
from glob import glob
import numpy as np
from numpy import r_, c_
import pandas as pd


#__________Function reads COORD file for miscelleneous data______________
def read_misc(filetxt, infoline, n_info):
    infoline=infoline.strip()#remove whitespace
    FWHM_start=infoline.find("FWHM = ")+len("FWHM = ")#index of where the numerical value of FWHM starts
    FWHM_end=infoline.find(" ppm")#index of where numerical value of FWHM ends
    SNR_start=infoline.find("S/N = ")+len("S/N = ")#index of where the numerical value of S/N starts
    
    #Do the same as above, just for the data shift line
    datashiftline=filetxt[n_info+1]
    datashiftline=datashiftline.strip()
    Datashift_start=datashiftline.find("Data shift = ")+len("Data shift = ")
    Datashift_end=datashiftline.find(" ppm")

    #Do the same for Ph. However will divide Ph into Ph1 and Ph2. Ph1 is in deg, Ph2 is in deg/ppm
    Phline=filetxt[n_info+2]
    Phline=Phline.strip()
    Phlist=Phline.split()#will break up line into each compoenent.the 2nd element is Ph1 value and the 4th is Ph2 value
    Phlist=[item for sublist in Phlist for item in sublist.split(":")]#sometimes the spaces in COORD file are wonky so records "Ph:X" as one item. Need to furhter split by ":"
    Phlist = [item.strip() for item in Phlist if item.strip()] #remove blank eleemnts



    FWHM_value=float(infoline[FWHM_start:FWHM_end])
    SNR_value=float(infoline[SNR_start:])
    Datashift_value=float(datashiftline[Datashift_start:Datashift_end])
    
    Ph1_value=float(Phlist[1])
    Ph2_value=float(Phlist[3])


    return pd.DataFrame({"FWHM (ppm)":[FWHM_value],"S/N":[SNR_value],"Data Shift (ppm)":[Datashift_value], "Ph1 (deg)":[Ph1_value], "Ph2 (deg/ppm)":[Ph2_value]})#returns data frame with FWHM, SNR, Datashift, Ph1, Ph2


#___________Function that reads the COORD file for metabolite data
def read_Metabolites(filetxt):
    for n, linetxt in enumerate(filetxt):#find data abt headers of Metabolite array
        if (("Conc." in linetxt.strip().split()) and ("%SD" in linetxt.strip().split()) and ("Metabolite"in linetxt.strip().split())):
            header_index = n #header_index is the index of the line in which the headers (Conc.,  %SD, /Cr+PCr, Metabolite)for the metabolite array is found. 
            header_line = linetxt #header_line is content of the line^, so in other words the headrs
            break

    for n, linetxt in enumerate(filetxt):#find ending of Metbolite array
        if 'lines in following misc. output table' in linetxt.strip():
            end_index = n #end_inddex is the index of the last line of the metbolite data.
            break

    headers=[] #empty list to contain the headers
    headers+=header_line.strip().split() #appends to the list all the header names

    metabolite_data=[]#empty list that will contain all metabolite data
    for aline in filetxt[header_index+1:end_index]:#appends data to metabolite_data from the lines that contain metabolite data
        metabolite_data+=aline.strip().split()#note metabolite_data is one list, which contains conc, SD%, /Cr+PCD, Metabolite. SO need to reorganize

    #Reorganizing here. Making list of each collumn of the array
    Conc=metabolite_data[0::4]
    SD_perc=metabolite_data[1::4]
    Cr_PCD_ratio=metabolite_data[2::4]
    metabolite_name=metabolite_data[3::4]

    #Create data frame with metabolite data.
    return pd.DataFrame(list(zip(Conc,SD_perc,Cr_PCD_ratio,metabolite_name)),columns=headers), metabolite_name #returns data frame of metabolite data. Also returns all the names of metabolites to be used as column names for the coord data



    
#___________Function that reads the COORD file for the data points
def read_coord(coordfileloc): 
#this function only returns the coordinates of the spectra points
    print('Converting '+coordfileloc+'...')

    with open(coordfileloc) as coordfile:
        filetxt = coordfile.readlines()#list of all the lines

    for n, linetxt in enumerate(filetxt):#n is the line location, linetxt is content of that line
        if linetxt.strip().split(' ')[0] == 'FWHM':
            n_info = n #n_info is the line in which FWHM is found. Note SNR is also this line
            infoline = linetxt
            break

    first_header=n_info+6 #this is the index of the first header which reads "800 points of ppm axis"


    
    misc_df=read_misc(filetxt, infoline, n_info)# data.frame with the misc values. The data frame is made from read_misc() fucniton
    Metabolites_df, metabolite_names =read_Metabolites(filetxt)#data frame with Metabolite infos+list of metabolite names. Called from funciton

    metabolite_headers=["ppm","RawData","ProcessedData", "Background"]
    metabolite_headers+=[item for item in metabolite_names if not '+' in item]#remove any metabolites which are added together (include 'x')
    #this lsit will contian the headers for the csv file

    plotsize = int(filetxt[first_header].strip().split(' ')[0])#Differes, in orginal data was 800
    linestoread = int(np.ceil(plotsize/10))#calculates lines to read

    for n, linetxt in enumerate(filetxt):#find ending of the data using "lines in following diagnostic table"
        if 'lines in following diagnostic table:' in linetxt.strip():
            last_index = n #index of the last line of the data
            break

    for n, linetxt in reversed(list(enumerate(filetxt))):#find ending of the data looking at the last line that has 10 different elements that are seperated by spaces. (because the data is 10 different numbers seperated by spaces). For further security also check if the lien contains more than 100 characters. Going in the reverse way becasue want to start from end of scipt
        if (len(linetxt.strip().split())==10 and len(linetxt.strip())>=60):
            if plotsize%10!=0:#depending on #of data points. If plot size not divisivle by 10, will not take the last line but the second to alst data point
                last_index_check = n+2#index of the last line of the data
                break
            elif plotsize%10==0:
                last_index_check=n+1
                break
            else:
                print("ERROR in line indexing")


    if last_index_check==last_index:
        print("Endpoint of data is confirmed in row "+ str(last_index_check))
    else:
        last_index=int(input("Data endpoint is not confirmed. From detecting 'lines in following diagnostic table' endpoint is "+(last_index) + "From detecting 10 elements in a line endpoint is "+ (last_index_check)+ ". Which row is correct? "))






    
    coorddata=pd.DataFrame(columns=metabolite_headers)#Create empty data frame with metabolite_headers as the headers.
    for column_name in metabolite_headers:
        zeros_list = np.zeros(plotsize, dtype=float)#create a list with plotsize 0's in float form. It's going to be the column for the metabolite that has conc=0 and is not recorded in the data
        coorddata[column_name]=zeros_list
    #so now have the base data frame with correct headers but filled with 0's. Next, We are going to change the 0's for the metbaolites/columns that include data
 

    ppmaxis = []#empty list to store the data points of the ppm (first column in csv file)
    for aline in filetxt[first_header+1:first_header+1+linestoread]:#n_info+7 is the start of the data. n_info+7+linestoread is the end of the data. #aline is the contents of each individual line of each 80 lines of 800 data poitns
        ppmaxis += aline.strip().split()
    ppmaxis = np.array(ppmaxis)
    ppmaxis = ppmaxis.astype(float)#note this array (and the 3 folowing) have 1 collumn and 800 rows. Unlike the metfit arrays
    coorddata["ppm"]=ppmaxis#Change the first collumn

    data = []#does same as above, just stores the next 80 lines containing 800 data poitns.  (second column in csv file)
    for aline in filetxt[first_header+2+linestoread:first_header+2+2*linestoread]:
        data += aline.strip().split()
    data = np.array(data)
    data = data.astype(float)
    coorddata["RawData"]=data#changes raw Data column (right now nothing) to contian data


    fit = [] #Same as above*2. (third column is csv file)
    for aline in filetxt[first_header+3+2*linestoread:first_header+3+3*linestoread]: 
        fit += aline.strip().split()
    fit = np.array(fit)
    fit = fit.astype(float) 
    coorddata["ProcessedData"]=fit #changes Processed Data columng (right now nothing) to contian data

    if "background" in filetxt[first_header+3+3*linestoread]: #some data do not have background value
        bkg = [] #same as above*3. (fourth column in csv file)
        for aline in filetxt[first_header+4+3*linestoread:first_header+4+4*linestoread]: 
            bkg += aline.strip().split()
        bkg = np.array(bkg)
        bkg = bkg.astype(float)
        coorddata["Background"]=bkg #Edits the background column (currently empty) to contain the data points
        bkg_offset=0
    else:
        zeros_list = np.zeros(plotsize, dtype=float)#create a list with plotsize 0's in float form. It's going to be the column for the background 
        coorddata["Background"]=zeros_list
        bkg_offset=1



    visible_met=((last_index-first_header)/(linestoread+1))-4+bkg_offset#Find how many metabolites are listed for data. (or how many metabolites (exluding ones that contian"+") do not have Conc=0 or SD=999%). -4 because first 4 sets of 800 data points are the ppmaxis, raw data, etc
    print("There are {} visible metabolites".format(visible_met))


    met_start=first_header+(4-bkg_offset)*(linestoread+1) #the index of wherethe first metabolite header starts. THe metabolizte data starts on met_start+1
    for n in range(int(visible_met)):#similar to above fucntions, just a for loop that finds the data for each metabolite and appends it to the respective column in csv file
        met = []#empty list for each csv collumn. Will be created for all visible metabolite, but ofc not all metabolite_headers (bc some metbaolite_headers don't have data if Conc=0 and SD=999%)
        #note met_start+(linestoread+1)*n is the nth/ header, NOT the nth start to the data. nth start to data is metstart+(linestoread+1)*n+1
        header_name=filetxt[met_start+(linestoread+1)*n].strip().split()[0]#name of the visible_met that we are abt to read the data for (this should match the header of csv column)
        for aline in filetxt[met_start+(linestoread+1)*n+1:met_start+(linestoread+1)*(n+1)]:#from one after header (start of data) to end of the data for a the nth listed metabolite
            met += aline.strip().split()#create list of all the plotsize data poitns for the given metabolite
        met = np.array(met)
        met = met.astype(float)#this list contains the plotsize data poitns of this nth metabolite listed. Need to chjange into float outside for loop to avoid error from differnt elements in array  
        coorddata[str(header_name)]=met #first mathch the name of the metabolite to the correct column in csv, then append the data poitns to that column
         

    return coorddata,misc_df, Metabolites_df #returns a tuple with all three data frames (each one becomes an individual csv)




#____________Main function that takes the directory of the COORD files, calls the specific functions, and then creates all three csv files.
def main(coorddir): 
    if coorddir[-1] != '/':
        coorddir += '/'
    coords = glob(coorddir+'*.COORD') + glob(coorddir+'*.coord')#for some reason COORD folder has to be in desktop, it cant be within any other folder.
    print(str(len(coords)) + ' .COORD files detected')

    for coord in coords:#creats (using read_coord() function) csv files. Saves csv files in the same folder as original .coord files but with .csv 
        coorddata, miscdata, Metabolitesdata = read_coord(coord)

        #np.savetxt(coord.rsplit('.')[0]+'_coord.csv', coorddata, delimiter=',') #not using this only for arrays
        coorddata.to_csv(coord.rsplit('.')[0]+'_coord.csv', sep=',', header=True, index=False)#create csv file with headers but no
        print(coord.rsplit('.')[0]+'_coord.csv' + ' written')
        
        miscdata.to_csv(coord.rsplit('.')[0]+'_misc_coord.csv', sep=',', header=True, index=False)#create csv file for the FWHM and SNR
        print(coord.rsplit('.')[0]+'_misc_coord.csv' + ' written')   

        Metabolitesdata.to_csv(coord.rsplit('.')[0]+'_Metabolites_coord.csv', sep=',', header=True, index=False)#create csv file for the Metabolites
        print(coord.rsplit('.')[0]+'_Metabolites_coord.csv' + ' written')   


#Here is where we specify the director of the folder with the COORD files and call the main function
coorddir=str(input("What is the directiory of the folder with the coord files? "))#function that prompts user for location of coord files
#Note that for some reason the folder should be in desktop (like not nested in anotehr folder) + the name of the folder should not contain any special charafcters or spaces
print(coorddir)
main(coorddir)#calls function
