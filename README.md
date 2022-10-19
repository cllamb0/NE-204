# NE204 Lab 1: Pulse Processing  
Author: Megan Schiferl (mschiferl@berkeley.edu)

Credit to: Chris Lamb, Curtis Berger, and Jisu Park who helped me pull all this code together

Date: 10/19/2022

## Contents
 - Purpose
 - Relavent Files
 - Quick Procedure
 - Detailed Notes

## Purpose
To properly calibrate a gamma ray spectrum from a high purity germanium detector, these seven notebooks explore optimization of a trapezoidal filter and building a calibration curve. If your goal is only to produce a calibrated spectrum, you'll only need to use the 7th notebook. If your goal is to follow my learning and optimization process, follow along with this readme file. 

## Relavent Files
1-FinalNotes_Lab1_MS_NE204

2-m&kOptimizationRoutine_Lab1_MS_NE204.ipynb

3-m&kComparisons_Lab1_MS_NE204.ipynb
 - __You'll need:__ TrapHeights2, 25 csv files (contains k/m optimization heights)

4-HeightOptRoutine_Lab1_MS_NE204.ipynb

5-HeightComparisons_Lab1_MS_NE204.ipynb
 - __You'll need:__ TrapHeightsFitOpt, 3 csv files (contains height finding method heights)

6-CalibrationCurve_Lab1_MS_NE204.ipynb
 - __You'll need:__ FinalTrapHeights, 4 csv files (heights for calibration sources at optimized k/m and height finding method)

7-FinalRaw2Calib_Lab1_MS_NE204.ipynb
 - __You'll need:__ CalibrationParameters.csv

_TrapHeights2, TrapHeightsFitOpt, and FinalTrapHeights are all in the google drive._

## Quick Procedure

### Step 1: 
'1-FinalNotes_Lab1_MS_NE204' has notes on almost all steps going forward. Come back to this notebook for more details/visuals as needed.

### Step 2: 
'2-m&kOptimizationRoutine_Lab1_MS_NE204' runs an optimization procedure over peaking and gap time using Cs data. Utilizes the trapezoid fitting height finding method. Returns 25 csv files, one for each k-m combination. 

### Step 3:
'3-m&kComparisons_Lab1_MS_NE204' uses the 25 csv files from step 2 to determine the optimal peaking and gap time. 

### Step 4:
'4-HeightOptRoutine_Lab1_MS_NE204' runs an optimization procedure over the three height finding methods: max value, trapezoid fitting, and gradient. Returns 3 csv files, one for each method. Also, the max value method is run for each of the calibration sources (Cs-137, Co-60, Co-57, and Eu-152) and returns 4 csv files with the resultant heights.

### Step 5:
'5-HeightComparisons_Lab1_MS_NE204' uses the 3 csv files from step 4 to determine the optimal height finding method.

### Step 6:
'6-CalibrationCurve_Lab1_MS_NE204' uses the 4 csv files from step 4 to create a calibration curve and produce calibrated spectra for the calibration sources.

### Step 7:
'7-FinalRaw2Calib_Lab1_MS_NE204' takes in one raw data file and prints a calibrated spectrum based on the calibration curve from step 6.

## Detailed Notes (also found at the beginning of each notebook)

### Step 1:

__Purpose & Contents:__ This notebook contains notes on a lot of the functions that were tried as well as more visuals for each step of the process from raw pulses to finding the trapezoid height. All methods are explored, only some are used going forward.

#### Explores:
 - DAC Echo Correction
 - Saving Raw Data to Array
 - Removing Saturation
 - Background Subtraction
 - Applying a Savinsky-Golay Filter
 - Noise Height Investigation
 - Finding Start of Rise (i) and Rise Time (k) Versions 1 - 3
 - Finding Tau
 - Defining the Trapezoidal Filter
 - Applying the Trapezoidal Filter
 - Fitting Trapezoids
 - Finding the Trapezoid Height Max Value, Trapezoid Fit, and Gradient Methods

__To run this notebook:__ you'll need to change the directory and file name, but otherwise you can run all.


### Step 2:

__Purpose and contents:__ This notebook takes in raw pulses and returns the heights of signal trapezoids for 25 different combinations of gap and peaking time. This data will be used in Optimization Routine pt 2 for finding the ideal combination of gap and peaking time. This optimization routine finds trapezoid heights by fitting the signal trapezoid with a trapezoid.

 - Section 1: Imports and Directory
 - Section 2: Managable Chunks of Data
 - Section 3: Functions
 - Section 4: Optimization Routine

__Notes on Running:__ First of all - don't! I've uploaded the optimization routine results from this notebook to google drive, so just move on to the next step to use them. If you still want to try this out, as always make sure the directory and data file are correct, otherwise you can run it right through and you'll end up with 25 csv files containing trapezoid heights. It takes ~90 mins on my computer, but maybe yours is faster.


### Step 3:

__Purpose and contents:__ This notebook takes in the heights of signal trapezoids generated in the step 2 optimization procedure and outputs (non-calibrated) spectra and plots on optimization of peaking and gap time.

 - Section 1: Imports and Directory
 - Section 2: Functions and Parameters
 - Section 3: Visual Comparison of Spectra
 - Section 4: Analytical Comparison of Spectra (plus some notes on how the outputs of the gaussain fit work)

__Notes on Running:__ Ensure that the directory is correct for your system and that you are using the height data provided (otherwise you must manually find the peak for each k/m system).


### Step 4: 

__Purpose and Contents:__ This notebook takes in raw pulses and returns 3 csv files containing trapezoid height. The first file finds the heights using the maximum value method, the second uses a trapezoid fitting method, and the third uses a gradient method. This notebook is also used to save the heights of trapezoids using only the maximum value method for multiple sources (simply change the data file and don't run the last two cells)

  - Section 1: Imports and Directory
  - Section 2: Functions
  - Section 3: Calculating Trapezoids
  - Section 4: Height Finding Methods

__Notes on Running:__ Probably refrain from running this one too. I've uploaded the resulting 3 (height finding methods) + 4 (calibration heights) csv files to google drive, so you can just take those and move on to step 5. If you want to run it howeer, as usual, you'll need to check your directory and the imported data name, but otherwise you should be able to run all. This one takes about 10 mins to run, so it's not too bad.


### Step 5:

__Purpose and Contents:__ This notebook takes in the height finding data from step 4 and compares the methods for finding the height of the trapezoids. It is determined that the methods do not result in statistically significant photopeak channel locations or resolution, so the maximum height method is chosen going forward for it's speed.

 - Section 1: Imports
 - Section 2: Functions
 - Section 3: Plotting Height Finding Methods

__Notes on Running:__ You'll need to change to your directory and load in the 3 files from step 4, but otherwise, run all!

### Step 6:

__Purpose and Contents:__ This notebook takes in the heights for Cs137, Co-60, Co-57, and Eu-152 using max value method from step 4. This data is used to construct a calibration curve and return calibrated spectra for each of these isotopes.

 - Section 1: Imports
 - Section 2: Functions
 - Section 3: Find Photopeaks
 - Section 4: Create Calibration Curves
 - Section 5: Produce Calibrated Spectra

__Notes on Running:__ Change to your directory that contains the csv files with the heights and you should be able to run all if you use the csv files I uploaded. If you use other isotopes or other csv files, you may have to change the window around where the photopeaks are to get an accurate calibration curve.


### Step 7:

__Purpose and Content:__ To take in raw pulses from a high purity germanium detector and return a calibrated spectrum, using prevoiusly optimized filtering parameters and calibration curve.

 - Section 1: Imports
 - Section 2: Splitting the Data
 - Section 3: Defining Functions
 - Section 4: Calculating Trapezoid Heights from Raw Data (Using the Max Value Method)
 - Section 5: Calibrated Spectrum

__Note on Running:__ To run this notebook successfully, you'll need to check/change the following within the notebook:

 - In Section 1: Imports - the second cell is where you'll import the raw pulses from the detector and the calibration paramters. You'll need to change the directory and file names.

 - In Section 2: Splitting the data - depending on the number of pulses in the file, you'll need to change the number of pulsesX arrays to hold ~5000 pulses. This step was necessary for my (Megan's) computer to handle the large file sizes. Thus, the code is built around this structure. Section 1 has an output that will tell you the number of pulseX arrays to use.

 - Once you have these parts checked and changed, you can run the rest of the notebook through to get your spectrum!

----

Have fun, and good luck! If you have any questions, feel free to email me.
