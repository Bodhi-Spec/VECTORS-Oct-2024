Described below is the process to utilize the automated DQC pipleine.

To get the data in the correct format:
1) you must create a folder with all the COORD files of the data that you want to QA.
2) run the coord2csv python script, which will convert the COORD files into csv files (and output them to the same folder)

Provided are five R scripts:
- Three of which preform the metrics associated with each heading of the Methods section (Duplicate Peaks, Vertical Shifts, and Glx_mI_peaks). 
    - VerticalShifts reports the metrics anyNegative and belowBaseline
    - Duplicate Peaks reports the metric existDuplicate
    - Glx_mI_peaks report the metrics Glx Distinct, Glx Merge, mI Distinct, and mI Merge.
- A fourth script, plotgraph, can plot the data into a spectrum and is helpful for visualization, but is not necessary for the DQC pipeline. 
- Finally, the last script, GetData is the parent function that 1) loads the packages 2) reads in the csv files and 2) calls the other four R script functions. 
    - Example ways to call the other R functions are provided at the bottom of the script.









