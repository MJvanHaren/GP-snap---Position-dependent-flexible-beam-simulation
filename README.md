# GP snap - Position-dependent flexible beam simulation
Several files for performing FRF, ILCBF, GP snap feedforward identifcation and vaidation. This is done for a flexible beam setup, as in the DCT lab at TU/e.

## file structure
Several files are found in the folder. The files ILCBFSimscape.m and FRFSimscape.m are the files which should be run (directly from matlab).
Requires Matlab 2021a for the simulink files!

1. ModelFlexibleBeamPOIL.slx - model of flexible beam when POI is left of center
2. ILCBFSimscape.m - runs ILCBF on flexible beam
	- flexibleBeamILCBF.slx - simulink file for ILCBF on flexible beam
3. FRFSimscape.m - performs FRF on flexible beam, with output either center point or POI
	- flexibleBeamFRF.slx - simulink file for FRF on flexible beam
	- flexibleBeamFRF.mat - data of FRF (performed on 09-06-2021)
4. FBcontroller_09062021.mat - FB controller designed on 09-06-2021

## folders
In the folder 'functions' several files are found
1. make4.m - used for 4th order motion profile
2. profile4.m - used for 4th order motion profile
3. PlotTrialData.m - used for plotting trial data for ILC(BF)
4. SetPlotLatexStyle.m - set plots to latex font/painters renderer etc.