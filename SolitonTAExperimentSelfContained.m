classdef SolitonTAExperimentSelfContained < handle
    % Class for soltion probed TA from 716. Designed to mimic class
    % structure from the high rep-rate TA in 714
    %   Constructor requires three arguments:  
    %   folderPath {string} path to folder with TA Data, must end with "/"
    %   fileName {string} name of file containing data from the first scan
    %   numScans {numeric} number of scans to be combined
    % 

    % future improvements: 
    % 1) plot scan by scan 
    % 2) fix axis labels and titling 
    % 3) load specific scans


    properties
                % json;
        folderPath (1,1) string
        fileName (1,1) string
        energyAxis(1, :) {mustBeNonnegative}
        pumpWavelength (1, 1) {mustBeNumeric} = 0 % Units: nanometers
        pumpBandwidth (1, 1) {mustBeNumeric} = 0 % Units: nanometers
        pumpEnergy (1, 1) {mustBeNumeric} = 0 % Units: nanojoules
        opticalDensity (1, 1) {mustBeNumeric} = 0 % Units: dimensionless
        jetPathLength (1, 1) {mustBeNumeric} = 0 % Units: micrometers
        pumpAngle (1, 1) {mustBeNumeric} = 0 % Units: degrees
        pumpSpotSize (1, 1) {mustBeNumeric} = 0 % Units: micrometers
        probeSpotSize(1, 1) {mustBeNumeric} = 0 % Units: micrometers
        analytes Analyte
        solvent (1, 1) string = ""    
        nTimes (1, 1) {mustBeNumeric} = 0 
        timesRaw (:, 1) {mustBeNumeric}
        timesSorted (:, 1) {mustBeNumeric}
        nScans (1, 1) {mustBeNonnegative} = 0
        TACubeRaw (:, :, :) {mustBeNumeric}
        % TACubeUnmirrored (:,:,:) {mustBeNumeric}
        TACubeSorted (:,:,:) {mustBeNumeric}
        TAMeanRaw (:, :) {mustBeNumeric}
        % TAMeanUnmirrored (:,:) {mustBeNumeric}
        TAMeanSorted (:,:) {mustBeNumeric}
        pumpPolarization (1,1) {mustBeNumeric} = 0 % Units: degrees horizontal (P) is 0 degrees, vertical (S) is 90 degrees, magic angle is 54.7 deg between pump and probe but each polarization must be specified in the lab frame
        probePolarization (1,1) {mustBeNumeric} = 0 % Units: degrees
        backgroundSubtractionUpperBound (1,1) {mustBeNumeric};
        backgroundSubtractionUpperBoundInxed (1,1) {mustBeNumeric};
        background (:,1) {mustBeNumeric}
        backgroundStdev (:,1) {mustBeNumeric}
        TAMeanSortedBackgroundSub (:,:) {mustBeNumeric}
        TAMeanSortedBackgroundSub_CC (:,:) {mustBeNumeric}
        timesSorted_CC (:,1) {mustBeNumeric}
        energyAxis_CC (:,1) {mustBeNumeric}
        UVVis_unirradiated = struct('filePath','wavelength','eV','absorbance','varaince','notes'); %(:,:) {mustBeNumeric} % Two column array, first coulmn is wavelength, second column is intensity
        UVVis_irradiated = struct('filePath','wavelength','eV','absorbance','varaince','notes'); % custom structure with fields for wavelength, absorbance, and (optionally) variance and notes eg dilution details and path length (:,:) {mustBeNumeric} % Two column array, first coulmn is wavelength, second column is intensity
        chirpCorrectionCoefficients (5,1) {mustBeNumeric} 

        pumpSpectrum = struct('filePath','wavelength','eV', 'intensity');
        fluorescenceSpectrum = struct('filePath','wavelength','eV','intensity');
        % solitonSpectrum = struct('filePath','pixel', 'intensity');
        mainGraphTitle (1,1) string;
   end
    methods
        function obj = SolitonTAExperimentSelfContained(folderPath,fileName, nScans) %,pumpBandwidth,pumpEnergy,opticalDensity,jetPathLength,pumpAngle,pumpSpotSize,probeSpotSize,nTimes,nScans,timeAxis,calbiratedEnergyAxis,averagedTA,pumpPol,probePol,analytes, solvent)
            %SolitonTAExperiement Construct an instance of this class
            %   Pass arguments for folderName (must end with "/"), fileName 
            %   (the file of the first scan that will be loaded), nScans
            %   (the number of scans)
            %   pumpBandwidth, pumpEnergy, opticalDensity, jetPathLength,
            %   pumpAngle, pumpSpotSize, probeSpotSize,number of times,
            %   timeAxis, calibratedEnergyAxis, averagedTA, pump
            %   polarization, and probe polarization
            obj.fileName = fileName;
            obj.folderPath = folderPath;
            obj.nScans = nScans;
            % obj.pumpWavelength=pumpWavelength
            % obj.pumpBandwidth=pumpBandwidth
            % obj.pumpEnergy=pumpEnergy
            % obj.opticalDensity=opticalDensity
            % obj.jetPathLength=pumpAngle
            % obj.pumpAngle=jetPathLength
            % obj.pumpSpotSize=pumpSpotSize
            % obj.probeSpotSize=probeSpotSize
            % obj.nTimes=nTimes
            % obj.nScans=nScans
            % obj.times=timeAxis
            % obj.energyAxis=calbiratedEnergyAxis
            % obj.TAMean=averagedTA
            % obj.analytes=analytes
            % obj.solvent=solvent      
            % obj.pumpPolarization=pumpPol
            % obj.probePolarization=probePol
        end
        % functions to read in and process soliton TA data will be included
        % below (named like the ones that are external functions but
        % appended with "_SC" for self-contained

        % functions to be included: main, load, unmirror,
        % process,background subtract, and chirp correct
        function mainSoliton_SC(obj,energyAxis,isMirrored,upperBoundForBS,options)
        % main function to load in and process soliton probed TA data 
            % Arguments are folder path (string), file name (string), number of scans (integer or double), energyAxis(array of 1024 doubles),
            % isMirrored (boolean 0 for not mirrored, 1 for mirrored scan), upper
            % bound for background subtraction (double, pick a time point before any 2PA
            % signal is visible -1.5 is usually safe), and optionally a row of 5
            % coefficeints for a chirp correction
        
            % initialize a instance of SolitonTAExperiement that will hold your
            % data and experiment details. This requires the folder path ending in
            % "/" on mac, the file name of the first scan you want to load in, and
            % the number of scans in this folder. In the future this can be adapted
            % to include selected scans
            arguments
                % folderPath;
                % fileName;
                % numScans {mustBeNumeric};
                obj SolitonTAExperimentSelfContained
                energyAxis {mustBeNumeric};
                isMirrored {mustBeNumericOrLogical};
                upperBoundForBS {mustBeNumeric};
                options.chirpCorrectionCoefficients (1,5) {mustBeNumeric} = [0 0 0 0 0];
            end
            %DS=SolitonTAExperimentSelfContained(folderPath,fileName,numScans);
            loadSolitonScans_SC(obj); % read in data
            if isMirrored
                unmirrorScans_SC(obj)
            else
                processSolitonTimes_SC(obj) 
            end
            % need to add a calibrated axis or calibration step
            obj.energyAxis=energyAxis;
            % background subtract
            backgroundSubtractSoliton_SC(obj,upperBoundForBS);
            % chirp correct 
            chirpCorrectionSoliton_SC(obj,options.chirpCorrectionCoefficients)
        end
        function loadSolitonScans_SC(obj)
            % obj must be a SolitonTAExperimentSelfContained, isMirrored is a boolean (0 for
            % scans that are not mirrored 1 for scans that are mirrored) 
            scanName='scan' ; % hard-coded in, maybe add it as an option in the future
            % dataFolder = DS.folderPath;
            datapath = strcat(obj.folderPath, obj.fileName);
            cd(obj.folderPath);
            sizing = readmatrix(datapath);
            n_t = size(sizing,1); % number of time points + 1 (includes 0 at (1,1)) 
            n_p = size(sizing,2); % number of pixels + 1
            clear sizing;
        
            datCube = zeros(n_t,n_p,obj.nScans);
        
            for j=1:obj.nScans
                file_n = strrep(datapath,[scanName '1'],[scanName num2str(j)]);
                datCube(:,:,j) = readmatrix(file_n);
            end
        
            datCube(isinf(datCube)) = 0 ; %removes any cells that have infinite values
            datCube(isnan(datCube)) = 0 ; %removes any NaNs (not a number)
            timeAxis = mean(datCube(:,1,:),3);
            obj.timesRaw=timeAxis;
            datCube = datCube(:,2:end,:);
            obj.TACubeRaw = datCube;
            obj.TAMeanRaw=mean(datCube,3);
        end

        function unmirrorScans_SC(obj)
            tempcube = obj.TACubeRaw;
            [a, b, c]=size(tempcube);
            unmirroredCube= zeros(a/2,b,obj.nScans*2);
            for i=1:c
                A = tempcube(1:a/2,:,i);
                size(A);
                B = flipud(tempcube(a/2+1:end,:,i));
                size(B);
                % C = cat(3,A,B);
                unmirroredCube(:,:,2*i-1)=A;
                unmirroredCube(:,:,2*i)=B;
                % clear A, B
            end
            unmirroredTimes=obj.timesRaw(1:a/2);
            % Now that the unmirrored cube has been populated, search for
            % duplicates in the time axis and remove them/confirm that the time
            % axis is appropriately sorted
            [B,I,~] = unique(unmirroredTimes); %examines time delays and sorts into ascending pump-probe time delays
            t = B; %resorts time delays
            sortedUnmirroredTACube= unmirroredCube(I,:,:); %resorts data matrix into ascending order
            size(sortedUnmirroredTACube)
            clear B I;
        
            obj.TACubeSorted = sortedUnmirroredTACube;
            obj.TAMeanSorted=mean(sortedUnmirroredTACube,3);
            obj.timesSorted=t;
            % clear tempcube, a, b, c,t,sortedUnmirroredTACube,unmirroredTimes,unmirroredCube;
        end
        function processSolitonTimes_SC(obj)
            tempcube = obj.TACubeRaw;
            [a, b, c]=size(tempcube);
            % unmirroredCube= zeros(a/2,b,obj.nScans*2);
            % for i=1:c
            %     A = tempcube(1:a/2,:,i);
            %     size(A);
            %     B = flipud(tempcube(a/2+1:end,:,i));
            %     size(B);
            %     % C = cat(3,A,B);
            %     unmirroredCube(:,:,2*i-1)=A;
            %     unmirroredCube(:,:,2*i)=B;
            %     % clear A, B
            % end
            rawTimes=obj.timesRaw;
            % Now that the unmirrored cube has been populated, we must search for
            % duplicates in the time axis and remove them/confirm that the time
            % axis is appropriately sorted
            [B,I,~] = unique(rawTimes); %examines time delays and sorts into ascending pump-probe time delays
            t = B; %resorts time delays
            sortedTACube= tempcube(I,:,:); %resorts data matrix into ascending order
            size(sortedTACube)
            clear B I;
        
            obj.TACubeSorted = sortedTACube;
            obj.TAMeanSorted=mean(sortedTACube,3);
            obj.timesSorted=t;
            % clear tempcube, a, b, c,t,sortedUnmirroredTACube,unmirroredTimes,unmirroredCube;
        end
        function cmap = divergingBlueWhiteRed_SC(~, n)
            % MATLAB-only diverging colormap for TA plots
            % Blue -> White -> Red
            if nargin < 2 || isempty(n)
                n = 256;
            end

            if mod(n,2) ~= 0
                n = n + 1; % keep halves even
            end

            n1 = n/2;
            n2 = n/2;

            blue  = [33 102 172] / 255;
            white = [1 1 1];
            red   = [178 24 43] / 255;

            cmap1 = [ ...
                linspace(blue(1),  white(1), n1)' ...
                linspace(blue(2),  white(2), n1)' ...
                linspace(blue(3),  white(3), n1)' ];

            cmap2 = [ ...
                linspace(white(1), red(1), n2)' ...
                linspace(white(2), red(2), n2)' ...
                linspace(white(3), red(3), n2)' ];

            cmap = [cmap1; cmap2];
        end
        function backgroundSubtractSoliton_SC(obj, upperBoundTime)
            [UBT_idx,UBT_val]=obj.findNearestIV_SC(obj.timesSorted,upperBoundTime);
            backgroundAvg=mean(obj.TAMeanSorted(1:UBT_idx,:),1);
            bckgrndSubtracted=obj.TAMeanSorted-backgroundAvg;
            % figure()
            % surf(bckgrndSubtracted,"EdgeColor","none");clim([-0.005,0.005]);view(2);xlim([400,920])
            % figure()
            % surf(bckgrndSubtracted./DS.TAMeanSorted,"EdgeColor","none");clim([-0.005,0.005]);view(2);xlim([400,920])
            obj.background=backgroundAvg;
            obj.TAMeanSortedBackgroundSub=bckgrndSubtracted;
            obj.backgroundSubtractionUpperBound=UBT_val;
            obj.backgroundStdev=std(obj.TAMeanSorted(1:UBT_idx,:),1);
            obj.backgroundSubtractionUpperBoundInxed=UBT_idx;

        end
        function chirpCorrectionSoliton_SC(obj,CCcoeffs)
            % chirpCorrectionSoliton_SC Corrects for dispersion in probe by
            % 
            arguments
                obj SolitonTAExperimentSelfContained
                CCcoeffs (1,5) {mustBeNumeric}; %= [0,0,0,0,0];
            end
             
            timeAx=obj.timesSorted;
            newTimeAx=timeAx(1:end);
            data=obj.TAMeanSortedBackgroundSub;
            energyAx=obj.energyAxis;
            CCenergyAx=1:size(energyAx,2);
            timeWindow=[-2,2]; % window for polynomial fit
            [~,TimeInx] = min(abs(timeAx-timeWindow));
            if sum(CCcoeffs)==0
            % plot data
                figure();
                [~,h]=contourf(CCenergyAx,newTimeAx(TimeInx(1):TimeInx(2)),data(TimeInx(1):TimeInx(2),:),100);set(h,'Linestyle','none');colormap(turbo)
                colorbar()
                % % fixes white to zero
                cmax = .002;%max(data1(TimeInx(1):TimeInx(2),:),[],'all');
                cmin = -.002;%min(data1(TimeInx(1):TimeInx(2),:),[],'all');
                xlim([300,1000]);
                colormap(obj.divergingBlueWhiteRed_SC(200));
                clim([cmin cmax])
                % User selects points
                [x,y] = getpts;
                % fourth order polynomial is used to fit the chirp
                f = polyfit(x,y,4);
                % coeffs = coeffvalues(f);
                coeffs = f;
                p1 = coeffs(1);
                p2 = coeffs(2);
                p3 = coeffs(3);
                p4 = coeffs(4);
                p5 = coeffs(5);
                obj.chirpCorrectionCoefficients=coeffs;
                fitted = p1*CCenergyAx.^4 + p2*CCenergyAx.^3 + p3*CCenergyAx.^2 + p4*CCenergyAx + p5;
                figure();
                % plot(x,y);hold on;plot(f);plot(newprobeax,fitted);
                plot(x,y);hold on;plot(CCenergyAx,fitted); 
                %clear p1 p2 p3 p4 p5 coeffs x y f
            else
                fitted=CCcoeffs(1)*CCenergyAx.^4 + CCcoeffs(2)*CCenergyAx.^3 + CCcoeffs(3)*CCenergyAx.^2 + CCcoeffs(4)*CCenergyAx + CCcoeffs(5);
            end
        
            
            %Creates empty matrix with extra time delay coloumns- to allow data to be
            %shifted according to chirp?
            %new_t=t2(1:end);
            newmat = zeros(size(data));
            
            % Loops through each coloumn of data (list time delays for each pixel)
            for i =1:size(data,2)
                % Transpose data into a row of delays
                tmpdat = data(:,i)';
                % Shifts delays according to fited polynomial 
                t_i = newTimeAx - fitted(i);
                % Interpolate in time. Finds the data values only at the existing time
                % 
                vq = interp1(t_i,tmpdat,newTimeAx,'linear',0); 
                newmat(:,i) = vq';
                % clear tmpdat t_i vq 
            end
            
            chirpcorrdat = newmat;
            figure;[h]=surf(energyAx,newTimeAx,chirpcorrdat,EdgeColor="none");set(h,'Linestyle','none');ylabel('t / ps');xlabel('Probe Energy /eV');title('Chirp Corrected');ylim([-1 1]);colormap(obj.divergingBlueWhiteRed_SC(200));
            clim([-0.005,0.005]);colorbar();view(2);%xlim([1,6]);
            % TAdat_cc = chirpcorrdat;
            % t = newTimeAx;
            obj.TAMeanSortedBackgroundSub_CC=chirpcorrdat;
            obj.timesSorted_CC=newTimeAx;
            obj.energyAxis_CC=CCenergyAx;
        end
        

        % Plotting Functions

        % function  contourplot(obj) issue with multiple 0's in the x axis
        % - need to clean this up by trimming the matrix 
        %     %METHOD1 Summary of this method goes here
        %     %   Detailed explanation goes here
        %     fig=figure();view(2);
        %     colormap(obj.divergingBlueWhiteRed_SC(200))
        %     contourf(obj.energyAxis,obj.times,obj.TAMean,200,'EdgeColor','none');
        %     clim([-0.005 0.005]);colorbar();fontsize(14,"points");
        %     xlabel("Probe Energy /eV");ylabel("Time Delay /ps");
        %     xlim([0.8,6.5]);
        % end

        function  surfplotRaw(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            fig=figure();
            colormap(obj.divergingBlueWhiteRed_SC());
            surf(obj.energyAxis,obj.timesRaw,obj.TAMeanRaw,'EdgeColor','none');view(2);
            clim([-0.005 0.005]);colorbar();fontsize(14,"points");
            xlabel("Probe Energy /eV");ylabel("Time Delay /ps");
            % xlim([0.8,6.5]);
        end
        function  surfplotRawSingleScan(obj,scanNum)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            fig=figure();
            colormap(obj.divergingBlueWhiteRed_SC());
            surf(obj.energyAxis,obj.timesRaw,obj.TACubeRaw(:,:,scanNum),'EdgeColor','none');view(2);
            clim([-0.005 0.005]);colorbar();fontsize(14,"points");
            xlabel("Probe Energy /eV");ylabel("Time Delay /ps");
            % xlim([0.8,6.5]);
        end
        function  surfplotSorted(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            fig=figure();
            colormap(obj.divergingBlueWhiteRed_SC());
            surf(obj.energyAxis,obj.timesSorted_CC,obj.TAMeanSorted,'EdgeColor','none');view(2);
            clim([-0.005 0.005]);colorbar();fontsize(14,"points");
            xlabel("Probe Energy /eV");ylabel("Time Delay /ps");
            % xlim([0.8,6.5]);
        end
        function  surfplotSorted_BS_CC_SC(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            fig=figure();
            colormap(obj.divergingBlueWhiteRed_SC());
            surf(obj.energyAxis,obj.timesSorted_CC,obj.TAMeanSortedBackgroundSub_CC,'EdgeColor','none');view(2);
            clim([-0.005 0.005]);colorbar();fontsize(14,"points");
            xlabel("Probe Energy /eV");ylabel("Time Delay /ps");
            % xlim([0.8,6.5]);
            fontsize(14,'points');
        end
        

        function plotSolitonKineticsTraces_SC(obj,wvlnbndz,options)
        % Average and plot a range of wavelengths, wavelength bounds should be a 2D
        % array with start bound in row 1 and corresponding end bound in row two of
        % a given column
            %parameters
            arguments
                obj SolitonTAExperimentSelfContained;
                wvlnbndz (:,:) {mustBeNumeric};
                options.lines_ = 0;
                options.colormap_= nebula;
                options.numLegCols_ = 1;
                options.coloffset_ = 0;
                options.normalize_ = 0;
                options.waterfallOffset {mustBeNumeric} = 0;%0.001;
                options.ms_ {mustBeNumeric} = 12;
                options.mrkr_ =".";
            end
            colmap=options.colormap_;
            ms=options.ms_;
            mrkr=options.mrkr_;
            if options.lines_
                myLineStyle="-";
            else
                myLineStyle="none";
            end
            myOffset=options.waterfallOffset;
            coloffset=options.coloffset_;
            normalize=options.normalize_;
            figure()
            hold on
            for i = 1:size(wvlnbndz,2)
                if size(wvlnbndz,1)==2
                    [w1_idx,w1_val]=findNearestIV_SC(obj,obj.energyAxis,wvlnbndz(1,i));
                    [w2_idx,w2_val]=findNearestIV_SC(obj,obj.energyAxis,wvlnbndz(2,i));
                    idxs=sort([w1_idx,w2_idx]); % must be sorted so that lower index is called first in the next line
                    toplot=mean(obj.TAMeanSortedBackgroundSub_CC(:,idxs(1):idxs(2)),2);      
                    tempname = [num2str(round(w1_val,2)), ' to ', num2str(round(w2_val,2)),' eV'];
                else
                    [w1_idx,w1_val]=findNearestIV_SC(obj,obj.energyAxis,wvlnbndz(1,i));
                    toplot=obj.TAMeanSortedBackgroundSub_CC(:,w1_idx);
                    tempname = [num2str(round(w1_val,2)),' eV'];
                end
                if normalize
                    plot(obj.timesSorted,toplot/max(toplot) + (i-1)*myOffset,'DisplayName',tempname,Marker=mrkr,MarkerSize=ms,LineStyle=myLineStyle, Color=colmap(i*floor(255/(size(wvlnbndz,2)+coloffset)),:));
                else
                    plot(obj.timesSorted,toplot + (i-1)*myOffset,'DisplayName',tempname,Marker=mrkr,MarkerSize=ms,LineStyle=myLineStyle,Color=colmap(i*floor(255/(size(wvlnbndz,2)+coloffset)),:));
                end
            end
            if normalize
                title("Normalized Kinetics Traces");ylabel("Normalized Signal");%ylim([-1.2,1.2]);
            else
                ylabel("\Delta OD");ylim([-0.005,0.005]);
            end
            legend("show",Location="best",NumColumns=options.numLegCols_);xlabel("Time Delay /ps");
            xlim([-10,500]);
            fontsize(14,'points');
            if strlength(obj.mainGraphTitle)>0
                title(obj.mainGraphTitle)
            end
            subtitle("Kinetics Traces")
        end

        function plotSolitonSpectraAtSpecificDelays_SC(obj,timebndz,options)
            % Average and plot a range of wavelengths, wavelength bounds should be a 2D
            % array with start bound in row 1 and corresponding end bound in row two of
            % a given column
            %parameters
            arguments
                obj SolitonTAExperimentSelfContained; 
                timebndz; 
                options.lines_ {mustBeNumericOrLogical} = 0; 
                options.colormap_ = nebula; 
                options.ms_ = 12;
                options.mrkr_ = '.';
                options.waterfallOffset {mustBeNumeric} = 0;%0.001;
                options.normalize_ = 0;
                options.plotUnirradiatedAbsSpec = 0;
                options.plotUnirradiatedAbsSpec_scalfactor = 0.02;
                options.numLegCols_ {mustBeNumeric} =1;
            end
            colmap=options.colormap_;
            ms=options.ms_;
            mrkr=options.mrkr_;
            if options.lines_
                myLineStyle="-";
            else
                myLineStyle="none";
            end
            myOffset=options.waterfallOffset;
            figure()
            hold on
            for i = 1:size(timebndz,2)
                if size(timebndz,1)==2
                    [w1_idx,w1_val]=findNearestIV_SC(obj,obj.timesSorted_CC,timebndz(1,i));
                    [w2_idx,w2_val]=findNearestIV_SC(obj,obj.timesSorted_CC,timebndz(2,i));
                    numPts=w2_idx-w1_idx+1;
                    sortedIdx=sort([w1_idx,w2_idx]);
                    toplot=mean(obj.TAMeanSortedBackgroundSub_CC(sortedIdx(1):sortedIdx(2),:),1);
                    tempname = [num2str(round(w1_val,2)), ' to ', num2str(round(w2_val,2)),' ps, ',num2str(numPts),' points'];
                else
                    [w1_idx,w1_val]=findNearestIV_SC(obj,obj.timesSorted_CC,timebndz(1,i));
                    toplot=obj.TAMeanSortedBackgroundSub_CC(w1_idx,:);
                    tempname = [num2str(round(w1_val,2)),' ps'];
                end
                if options.normalize_
                    plot(obj.energyAxis,toplot/max(toplot)+(i-1)*myOffset,'DisplayName',tempname,Marker=mrkr,MarkerSize=ms,LineStyle=myLineStyle,Color=colmap(i*floor(255/size(timebndz,2)),:));
                else
                    plot(obj.energyAxis,toplot+(i-1)*myOffset,'DisplayName',tempname,Marker=mrkr,MarkerSize=ms,LineStyle=myLineStyle,Color=colmap(i*floor(255/size(timebndz,2)),:));
                end
                

            end
            if options.plotUnirradiatedAbsSpec
                tempname=strcat("GroundStateUnirradiatedSample * -", num2str(options.plotUnirradiatedAbsSpec_scalfactor));
                plot(obj.UVVis_unirradiated.eV,-1*options.plotUnirradiatedAbsSpec_scalfactor*obj.UVVis_unirradiated.absorbance,DisplayName=tempname); 
            end
            legend("show",Location="best",NumColumns=options.numLegCols_);xlabel("Probe Energy /eV");
            if options.normalize_
                ylabel("Normalized Signal");     
            else
                ylabel("\DeltaOD");
            end
            fontsize(14,'points');% xlim([1.1,6.2]);
            if strlength(obj.mainGraphTitle)>0
                title(obj.mainGraphTitle)
            end
            subtitle("Spectral Slices at Specific Delays")
        end

        
        function plotAbsSpectraPreAndPostIrradiation(obj,options)
            arguments
                obj SolitonTAExperimentSelfContained
                options.energyUnits = "eV";
            end
            if options.energyUnits == "eV"
                figure()
                plot(obj.UVVis_unirradiated.eV,obj.UVVis_unirradiated.absorbance,'k',obj.UVVis_irradiated.eV,obj.UVVis_irradiated.absorbance,'r')
                legend(["unirradiated sample" "irradiated sample"]);
                xlabel("Energy /eV");ylabel("Absorbance");
            else
                figure()
                plot(obj.UVVis_unirradiated.wavelengths,obj.UVVis_unirradiated.absorbance,'k',obj.UVVis_irradiated.wavelengths,obj.UVVis_irradiated.absorbance,'r')
                legend(["unirradiated sample" "irradiated sample"]);
                xlabel("Wavelength /nm");ylabel("Absorbance");
            end
        end

        % Saving functions
        function myBinnedTrace = getBinnedTimeTraceForFitting_SC(obj,wavelengthbounds)
            [Eidx1,~]=findNearestIV_SC(obj,obj.energyAxis,wavelengthbounds(1,1));
            [Eidx2,~]=findNearestIV_SC(obj,obj.energyAxis,wavelengthbounds(2,1));
            idxs=sort([Eidx1,Eidx2]);
            myBinnedTrace = mean(obj.TAMeanSortedBackgroundSub_CC(:,idxs(1):idxs(2)),2);
        end
        function [myMat]=saveMyTAForGlotaran_SC(obj,options)%toPath,asFileName, saveChirpCorrected)
            % take the TAMean matrix from the data structure and add a coulmn on
            % the left for wavelength and a row on the top time delay. Then save
            % that matrix as a .csv to be loaded into pyglotaran
            arguments
                obj;
                options.toPath= obj.folderPath;
                options.asFileName = "TAMeans_BS_for_PyGloTarAn"; 
                options.saveChirpCorrected = 0; % defaults to saving background subtracted sorted TA Means matrix, but NOT chirp corrected
            end
      
            if options.saveChirpCorrected
                myTA=obj.TAMeanSortedBackgroundSub_CC';
                options.asFileName=strcat(options.asFileName,"_CC"); %% appends "_CC" to file name to indicate that the chirp corrected version is saved
            else
                myTA=obj.TAMeanSortedBackgroundSub'; % Transpose my matrtix so that each column is a delay and each row is a wavelength
            end  
            myDelays=obj.timesSorted'; % make the sorted time delays a row that will be added at the top of the matrix
            myWavelengths=[NaN,obj.energyAxis]'; % Make energy axis that will be a coulmn on the 
            tempMat=cat(1,myDelays,myTA); % add a row for the delays at the top of the matrix
            myMat=cat(2,myWavelengths,tempMat);
            cd(options.toPath);
            writematrix(myMat,options.asFileName);

        end

        % Helper functions 
        function [myIndex,myActualValue] = findNearestIV_SC(~,myArray,myValue)
            % Finds the index of the elemenet nearest in value to a passed target
            % value. Returns [idx,value]as indicated by the "IV" appended to
            % the function name
            differences=abs(myArray-myValue);
            [~,idx]=min(differences);
            myIndex=idx;
            myActualValue=myArray(idx);
        end

        function importAbsSpectrum(obj, filePath, isIrradiated, options)
            % function to load in absorbance spectra for unirradiated and
            % irradiated samples
            %   
            arguments
                obj SolitonTAExperimentSelfContained;
                filePath;
                isIrradiated {mustBeNumericOrLogical};
                options.notes;
                options.containsVariance =1; %assumes that the file contains variances as the third coulmn - this is true for .csv's from 717
            end
            tempData=readmatrix(filePath);
            if isIrradiated
                obj.UVVis_irradiated.filePath=filePath;
                obj.UVVis_irradiated.notes=options.notes;
                obj.UVVis_irradiated.wavelengths = tempData(:,1);
                obj.UVVis_irradiated.eV = 1240./tempData(:,1);
                obj.UVVis_irradiated.absorbance = tempData(:,2);
                if options.containsVariance
                    obj.UVVis_irradiated.variance = tempData(:,3);
                end
            else
                obj.UVVis_unirradiated.filePath=filePath;
                obj.UVVis_unirradiated.notes=options.notes;
                obj.UVVis_unirradiated.wavelengths = tempData(:,1);
                obj.UVVis_unirradiated.eV = 1240./tempData(:,1);
                obj.UVVis_unirradiated.absorbance = tempData(:,2);
                if options.containsVariance
                    obj.UVVis_unirradiated.variance = tempData(:,3);
                end
            end
        end
        function importPumpSpectrum(obj, filePath) % add notes section
            tempdata=readmatrix(filePath);
            obj.pumpSpectrum.wavelength=tempdata(:,1);
            obj.pumpSpectrum.eV=1240./obj.pumpSpectrum.wavelength;
            obj.pumpSpectrum.intensity=tempdata(:,2);
        end
        function importFluorescenceSpectrum(obj, filePath) % add notes section for solution details
            tempdata=readmatrix(filePath);
            obj.fluorescenceSpectrum.wavelength=tempdata(:,1);
            obj.fluorescenceSpectrum.eV=1240./obj.fluorescenceSpectrum.wavelength;
            obj.fluorescenceSpectrum.intensity=tempdata(:,2);
        end
        
        function saveSolitonTAExperiment_SC(obj,filename,options)
            arguments
                obj SolitonTAExperimentSelfContained;
                filename ;
                options.dirToSaveto;
            end
        end

        function makeMainGraphTitle_SC(obj,options)
            arguments
                obj SolitonTAExperimentSelfContained;
                options.Sample;
                options.Solvent;
                options.Concentration;
                options.pumpEnergy = obj.pumpEnergy;
                options.pumpWavelength = obj.pumpWavelength;
            end
            tempStr=strcat(options.Sample, " ", num2str(options.Concentration), " mM, in ",options.Solvent, " excited by ", num2str(options.pumpEnergy), " nJ ",num2str(options.pumpWavelength));
            obj.mainGraphTitle=tempStr;
        end
        % function jacobian_nm_to_eV(obj,intensity_nm,nmAx)
        %     % takes a spectrum (histogram) in wavelength space and
        %     % transforms it into eV space
        %     % This requires using the jacobian b/c the relationship between
        %     % energy and wavelength is reciprocal 
        %     eVAx=1239.8/
        %     intensity_eV=intensity_nm *1239.8/
        % end
        % function jacobian_eV_to_nm(obj,intensity,eVax)
        % end
        % Setter functions
        function set.timesRaw(obj,timeAxis)
            obj.timesRaw=timeAxis;
        end
        function set.TACubeRaw(obj,TACube)
            obj.TACubeRaw=TACube;
        end
        function set.TAMeanRaw(obj,TAMean)
            obj.TAMeanRaw=TAMean;
        end
        function set.energyAxis(obj,energyAxis)
            obj.energyAxis=energyAxis;
        end
        % function set.TACubeUnmirrored(obj,TACubeUnmirrored)
        %     obj.TACubeUnmirrored=TACubeUnmirrored;
        % end
        % function set.TAMeanUnmirrored(obj,TAMeanUnmirrored)
        %     obj.TAMeanUnmirrored=TAMeanUnmirrored;
        % end
        function set.timesSorted_CC(obj,timesUnmirrored)
            obj.timesSorted_CC=timesUnmirrored;
        end
        function set.TAMeanSorted(obj,TAMeanSorted)
            obj.TAMeanSorted=TAMeanSorted;
        end
        function set.TACubeSorted(obj,TACubeSorted)
            obj.TACubeSorted=TACubeSorted;
        end
        function set.background(obj,background)
            obj.background=background;
        end
        function set.TAMeanSortedBackgroundSub(obj,TAMeanSortedBackgroundSub)
            obj.TAMeanSortedBackgroundSub=TAMeanSortedBackgroundSub;
        end
        function set.TAMeanSortedBackgroundSub_CC(obj,TAMeanSortedBackgroundSub_CC)
            obj.TAMeanSortedBackgroundSub_CC=TAMeanSortedBackgroundSub_CC;
        end
        function set.energyAxis_CC(obj,energyAxis_CC)
            obj.energyAxis_CC=energyAxis_CC;
        end
        function set.chirpCorrectionCoefficients(obj,chirpCorrectionCoefficients)
            obj.chirpCorrectionCoefficients=chirpCorrectionCoefficients;
        end
        function set.backgroundSubtractionUpperBound(obj,backgroundSubtractionUpperBound)
            obj.backgroundSubtractionUpperBound=backgroundSubtractionUpperBound;
        end
        function set.pumpEnergy(obj,pumpE)
            obj.pumpEnergy=pumpE;
        end
        function set.pumpWavelength(obj,pumpWavelength)
            obj.pumpWavelength=pumpWavelength;
        end
        function set.mainGraphTitle(obj,mainGraphTitle_)
            obj.mainGraphTitle=mainGraphTitle_;
        end
    end
end