# MPFI-Neuroimaging-Workshop
Thanks for joining our 2018 MPFI Neuroimaging Workshop this year! This code repository contains:
<br />**1.)** A guide to setup ImageJ, and allow this program to interface with MATLAB. I have included a guide for setting up ImageJ for Windows and OSX users. 
<br />**2.)** MATLAB code to allow you to fully process raw two-photon calcium imaging data from start to finish. The code will allow you to:
<br /> - read in multi-page TIFF stacks into MATLAB.
<br /> - use image registration to eliminate motion artifacts in an entire imaging stack. The included image registration algorithms are based on either Manuel Guizar's efficient subpixel image registration by cross-correlation (https://www.mathworks.com/matlabcentral/fileexchange/18401-efficient-subpixel-image-registration-by-cross-correlation) or Theo Walker's faster image registration with cross-correlation (https://github.com/naroom/downsampleReg).
<br /> - segment out cellular ROIs in a semi-automated way with the Cell Magic Wand tool written by Theo Walker
<br /> - extract out fluorescence timecourses for both our cells and the surrounding neuropil.
<br /> - convert raw fluorescence timecourses to ΔF/Fo. Fo is accurately computed by via a moving percentile-filter.
<br /> - basic analysis code that will allow you to inspect the individual response traces of each cell, the impact of neuropil contamination on cellular responses, the degree of tuning in each cell for our visual stimulus (here the the stimulus is a drifting grating of different orientations), and the fall-off of response correlations between neurons (including signal and noise correlations).

If you have any questions, please feel free to email me at david.whitney@mpfi.org.
