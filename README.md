# GP snap - Position-dependent flexible beam simulation
Several files for performing FRF, ILCBF, GP snap feedforward identifcation and vaidation. This is done for a flexible beam setup, as in the DCT lab at TU/e.

## Requirements
For both these toolboxes, please refer to https://github.com/MJvanHaren/Additonal-MATLAB-Functions.git 
 1. GPML Toolbox from Rasmussen
 2.  'General Functions' - some generic functions for MATLAB

## file structure
Several files are found in the folder. The files ILCBFSimscape.m and FRFSimscape.m are the files which should be run (directly from matlab).
Requires Matlab 2021a for the simulink files!

1. GPSnap.m - MAIN FILE executes full framework and validates (ILCBF, GP, evaluation)
1. ModelFlexibleBeamFirstPrinciple.m - model of first principles free-free beam
2. FBcontroller_firstPrinciplesBeam.mat  - shapeit file for first principles free-free beam feedback controller
3. flexibleBeamILCBF.slx - simulink file for simulating either first principles or simscape beam (only used for ILCBF Right now)
4. (NOT USED RIGHT NOW) ILCBFSimscape.m - runs ILCBF on flexible beam 
	- flexibleBeamILCBF.slx - simulink file for ILCBF on flexible beam
5. (NOT USED RIGHT NOW)FRFSimscape.m - performs FRF on flexible beam, with output either center point or POI
	- flexibleBeamFRF.slx - simulink file for FRF on flexible beam
	- flexibleBeamFRF.mat - data of FRF (performed on 09-06-2021)
6. (NOT USED RIGHT NOW) FBcontroller_09062021.mat - FB controller designed on 09-06-2021 for simscape beam
7. (NOT USED RIGHT NOW) ModelFlexibleBeamPOIx.slx - referenced simulink model for simscape beam

## folders
In the folder 'functions' the following additonal used functions are found
1. FFbeam.m - calculates mode shapes and eigen frequencies of free-free beam

