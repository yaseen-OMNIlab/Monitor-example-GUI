/*****************************************************************************
* | File      	:   Readme.txt
* | Author      :   Yuntao Li
* | Function    :   Help with use
* | Info        :   Email: li.yunt@northeastern.edu
*----------------
* |	This version:   V1.0
* | Date        :   2024-11-27
* | Info        :   For uploading to GitHub
******************************************************************************
Correlation analysis_MATLAB:
corr_data.mat: the dataset including 20 OISI and pupil diameter samples for correlation analysis
OISI_Correlation.m: The code for calculating cross correlation between delta_R/R and Pupil diameter

Data acquisition_Python:
CameraGUI.py: open one camera, record pupil, and review the frames after data acquisition
ADS1115.py: ADS1115 analog-to-digital converter configuration code.
CSI_Camera.py: Raspberry Pi Camera module configuration code.
GUI layout.png: A screenshot showing the function of CameraGUI.py.
Simple data: a 1-second frame collected from CameraGUI.py

Data pre-processing_MATLAB:
OISIandPupil.m: screen the motion and preview the results of OISI and pupil diameter.
Data_frational_change_all.mat and intensityMatrix.mat: example OISI results.
mag.csv: example acceleration and pupil diameter results.
******************************************************************************
