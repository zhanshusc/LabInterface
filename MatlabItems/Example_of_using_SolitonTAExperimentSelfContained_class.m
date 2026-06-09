%% Example_of_using_SolitonTAExperimentSelfContained_class
close all;clear;clc;startup;
%% load energy calibration
% calibration still has to be done elsewhere. In the future the calibration
% step
% can be added to this program
x=1:1024; % pass a better energy axis if you wish
%A= 419540         ;E= 13.7190         ;T0=1.7733e3        ; % path to calibration file here
A = 2.7639e+05;
E = 12.4698;
T0 = 1.3204e+03;
E_new_sqrd = E^2 - A./(x + T0);
E_new = real(sqrt(E_new_sqrd));
mask= E_new>0;
% cd(folderPath)
energyAxis=E_new.*mask;
%% Load data
% provide the path to the directory where data are stored - it must end
% with a forward slash "/"
folderPath= '/Users/azizmohammed/Desktop/Bradforth/TransientAbrosption/2026_03_14/Pyrazine_7p8mM_HPLCWater_73nJ_3H_MApump_ODjustunder0p5_20260314_309pm_8000shots_mirrored/';
% provide the name of the file containing scan 1 data
filename= 'scan1_unref_DOD.txt'
% define the number of scans that will be loaded in
numscans=1;
% Is the delay scheme mirrored? if so isMirrored=1, if not mirrored then =0
isMirrored=1;
% provide a time before which TA spectra will be averaged and subtracted
% from all other time points, given in ps
upperBoundForBackgroundSubtraction=-1.85;

% initialize instance of the class
pyrTA=SolitonTAExperimentSelfContained(folderPath,filename,numscans)%,energyAxis,isMirrored,upperBoundForBackgroundSubtraction);
% load and process data (including unmirroring, background subtraction, and
% chirp correction
mainSoliton_SC(pyrTA,energyAxis,isMirrored,upperBoundForBackgroundSubtraction)

% Example of providing chirp fit coefficients  
myChirpCorrCoeffs=pyrTA.chirpCorrectionCoefficients; % taking the coefficients from the above fit
pyrTA_userPassedCCC=SolitonTAExperimentSelfContained(folderPath,filename,numscans)
mainSoliton_SC(pyrTA_userPassedCCC,energyAxis,isMirrored,upperBoundForBackgroundSubtraction,"chirpCorrectionCoefficients",myChirpCorrCoeffs)

%% Load absorbance spectra for sample prior to and post irradiation

absfilenameUnirr="/Users/azizmohammed/Desktop/Bradforth/TransientAbrosption/2026_03_14/UNIRRADIATEDPYRAZINE_DILUTED100ULIN10ML_7P8MM_20260314_1CMPATH.CSV"
absnotes="Diluted 100uL (micropipettor) in 10 mL water (graduated cylinder), initial soln was 7.8 mM pyrazine (aq), 1 cm pathlength, unirradiated."
importAbsSpectrum(pyrTA,absfilenameUnirr,0,"containsVariance",1,"notes",absnotes)

absfilenameIrr="/Users/azizmohammed/Desktop/Bradforth/TransientAbrosption/2026_03_14/IRRADIATEDPYRAZINE_DILUTED100ULIN10ML_7P8MM_20260314_1CMPATH.CSV"
absnotesIrr="Diluted 100uL (micropipettor) in 10 mL water (graduated cylinder), initial soln was 7.8 mM pyrazine (aq), 1 cm pathlength, irradiated."
importAbsSpectrum(pyrTA,absfilenameIrr,1,"containsVariance",1,"notes",absnotesIrr)
% you can also load a fluorescence spectrum to overlap with stimulated
% emission
%% Load other information about scan

% sample information -- need to fix analytes structure


% all of these properties are required to be numeric or logical. if you
% don't know what a value is, set it equal to False


%
pyrTA.pumpWavelength = 267; % nm, central wavelength of the pump beam
pyrTA.pumpBandwidth = 1000; % nm, badnwdith at FWHM of pump spectrum peak
pryTA.pumpEnergy = 75; % nJ
pyrTA.opticalDensity = 0.46; %OD of jet
pyrTA.jetPathLength = 125; % microns --> estimated from OD, must include losses to reflectance (using water jet for comparison)
pyrTA.pumpAngle = 15 ; % degrees the angle between the pump and probe beams (made this value up) 
pyrTA.pumpSpotSize = 107; % microns FWHM of pump spot measured at the sample
pyrTA.probeSpotSize = 1000; % microns FHWM of probe spot measured at the sample
% POLARIZATION -- VERY IMPORTANT, defined in the lab frame (think unit
% circle where 0º is horizontal (parallel to the plane of the table (table
% surface and perpendicular to the surface normal of the table) and 90º is vertical ( perpendicular to the plane of the table
% which is also parallel to the surface normal of the table)

pyrTA.pumpPolarization=35; %º
pyrTA.probePolarization = 90; % soliton is always vertically polarized b/c fiber 1 output has to be vertical for chirped mirrors

% in the future we can add attributes for time and date of scan

%% Plotting



% look at surface map of TA data

surfplotRaw(pyrTA)
surfplotSorted(pyrTA)


% look at surface map of TA data, here BS means background subtracted, and
% CC means chirp corrected
surfplotSorted_BS_CC_SC(pyrTA)

% Inspect background used in background subtraction
figure()
plot(pyrTA.energyAxis,pyrTA.background)

% Inspect standard deviation of spectrum used in background average
figure()
plot(pyrTA.energyAxis,pyrTA.backgroundStdev)
% Inspect time points used in background average to get an idea of the
% noise

[idx,val]=findNearestIV_SC(pyrTA, pyrTA.timesSorted ,pyrTA.backgroundSubtractionUpperBound);
timesInBackground=pyrTA.timesSorted(1:idx)';
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,timesInBackground,"colormap_",colormapz.guppy);

times1=0.1:0.015:0.5;
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,times1,"colormap_",colormapz.turbo);


% examine absorbance spectra of sample before and after irradiation
plotAbsSpectraPreAndPostIrradiation(pyrTA)
plotAbsSpectraPreAndPostIrradiation(pyrTA,"energyUnits","nm")

% plot spectra at various times - single time points
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100])

% plot spectra at various times - averaged spectra in a given time range
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100;2,6,15,60,130])

% plot spectra at various times as above, now include absorbance specrtrum
% of sample for ground state bleach
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100;2,6,15,60,130],"plotUnirradiatedAbsSpec",1)

% plot spectra at various times as above, now include absorbance specrtrum
% of sample for ground state bleach, include a scaling factor to make
% spectrum roughly the same magnitude as the TA signal
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100;2,6,15,60,130],"plotUnirradiatedAbsSpec",1,"plotUnirradiatedAbsSpec_scalfactor",0.003)

% make it a waterfall plot using the keyword "waterfallOffset"
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100],"waterfallOffset",0.001)

% normalize spectra
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100],"normalize_",1)

% change the colormap, marker type, marker size, and add lines between
% points
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100],"colormap_",colormapz.guppy,"mrkr_",'x','ms_',4,'lines_',1)

% Note that you can also drop the quotes and use keyword=value in the
% function
plotSolitonSpectraAtSpecificDelays_SC(pyrTA,[1,5,10,50,100],colormap_=colormapz.guppy,mrkr_='x',ms_=4,lines_=1)

% plot kinetics traces 
plotSolitonKineticsTraces_SC(pyrTA,[2.1,3.1,4.66,5])

% plot kinetics traces integrated between energies or wavelengths
plotSolitonKineticsTraces_SC(pyrTA,[2.1,3.1,4.66,5;2.3,3.3,4.8,5.1])

% as above you can use the same keywords to make a waterfall plot, change
% the colormap, marker type, marker size, and whether or not lines connect
% the markers




