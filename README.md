# NE-204
Group project materials for Chris Lamb, Curtis Berger, Jisu Park and Megan Schiferl 

## Contents
 - Purpose
 - Relavent Files
 - Quick Procedure
 - Detailed Notes

## Purpose
To properly calibrate a gamma ray spectrum from a high purity germanium detector, the attached four notebooks address the four tasks listed in the lab1 manual. To recreate the results and to make sure you don't run into any bugs, it is important to read through the markdown cells carefully. For example, creating the trapezoids is time consuming and so they were save as dataframe.mat files upon completion. If you import the already built trapezoids, there is no need to run the cell that creates them in the first place. 

## Relavent Files
Lab1V2.ipnyb

  Pre processes the raw signals and transforms them into trapezoids. Peak finder algorithm finds the appropriate peaks and Performs First order linear        calibration and full energy calibration. 

Lab1_Task2.ipnyb

  Creates the functional relation between peak energy and FWHM resolution
  
Lab1_Task3.ipnyb

  Investigates Ballistic deficit
  
Lab1_Create_Calibrated_Spectrum.ipnyb

  Creates the calibrated by loading in the calibration file.
