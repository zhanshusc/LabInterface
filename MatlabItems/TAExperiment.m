classdef TAExperiment < handle
    % Class contains the information contained within a single
    % transient absorption experiment

    properties
       
        json;

        fileName (1,1) string

        nPixels (1, 1) {mustBeNumeric}
        pixels (1, :) {mustBeNonnegative}
        wavelengths(1, :) {mustBeNonnegative}

        Calibration (1, 2) {mustBeNumeric}
        CalibrationPoints (:, 1) CalibrationPoint

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
        
        nTimes (1, 1) {mustBeNumeric}
        times (:, 1) {mustBeNumeric}

        nScans (1, 1) {mustBeNonnegative}

        scans TAScan

        TAMean (:, :) {mustBeNumeric}
        TAVariance (:, :) {mustBeNumeric}
        TANShots (:, :) {mustBeNumeric}

        pumpOnMean (:, :) {mustBeNumeric}
        pumpOnVariance (:, :) {mustBeNumeric}
        pumpOnNShots (:, :) {mustBeNumeric}

        pumpOffMean (:, :) {mustBeNumeric}
        pumpOffVariance (:, :) {mustBeNumeric}
        pumpOffNShots (:, :) {mustBeNumeric}

        dispersionFitCoefficients (1, :) {mustBeNumeric}
        dispersionFit (1, :) {mustBeNumeric}
        
        % Temp
        temp
    end

    methods (Access = public)
        
        function obj = TAExperiment(fileName, options)
            
            arguments
                fileName (1, 1) string
                options.Experiment = 'Burritos';
            end

            obj.fileName = fileName;

            
            
            if isfile(obj.fileName)
                switch lower(options.Experiment)
                    case 'burritos'
                        obj.loadFile();
                    case 'specter'
                        obj.loadFile716();
                end
            else
                error('File not found');
            end

        end

        function subtractTABackground(obj, upperLimit)

            index = obj.findTime(upperLimit);

            for i = 1:obj.nScans

                background = mean(obj.scans(i).TAMean(1:index, :), 1);

                for j = 1:obj.nTimes

                    obj.scans(i).TAMean(j, :) = obj.scans(i).TAMean(j, :) - background;

                end
                
            end

            obj.combine_scans();
        end

        function correctDispersion(obj)

            correctionApp = DispersionCorrection(obj.pixels, obj.times, obj.TAMean);
            uiwait(correctionApp.UIFigure);
            obj.dispersionFitCoefficients = correctionApp.dispersionFitCoefficients;
            obj.buildDispersionFit();
            obj.dispersionCorrectData()
           
        end

        function [figTA, axTA, cbar] = peekTA(obj, options)
            arguments
                obj TAExperiment
                options.Units = "Wavelengths"
                options.TimeLimits = [min(obj.times), max(obj.times)];
            end

            figTA = figure();
            axTA = axes('Parent',figTA);
            if options.Units == "Wavelengths"
                surf(axTA, obj.wavelengths, obj.times, obj.TAMean, 'EdgeColor', 'none');
                xlabel(axTA, 'Wavelengths / nm');
                xlim(axTA, [min(obj.wavelengths), max(obj.wavelengths)]);
            else
                surf(axTA, obj.pixels, obj.times, obj.TAMean, 'EdgeColor', 'none');
                xlabel(axTA, 'Pixels');
                xlim(axTA, [min(obj.pixels), max(obj.pixels)]);
            end
            ylabel(axTA, 'Times / ps');
            ylim(axTA, options.TimeLimits);
            set(axTA, 'FontSize', 12);
            set(axTA, 'Box', 'On', 'LineWidth', 2);
            set(axTA, 'TickLength', [0.015, 0.015]);
            set(axTA, 'Layer', 'top', 'XGrid', 'off', 'YGrid', 'off');
            view(axTA, [0, 90]);
            title(axTA, obj.makeTitle());
            subtitle(axTA, obj.makeSubtitle());
            box(axTA, 'On');
            defaultColourLims = clim(axTA);
            maxColourLim = min(abs(defaultColourLims));
            clim(axTA, [-maxColourLim, maxColourLim]);
            if exist('brewermap', 'file') == 2
                colormap(axTA, brewermap([],'RdBu'));
            else
                colormap(axTA, 'parula');
            end
            cbar = colorbar(axTA);
            cbar.Label.String = "Transient Absorption / OD";
            cbar.Label.FontSize = 12;
            

        end

        function [figSpectra, axSpectra] = peekSpectra(obj, times, options)
            arguments
                obj TAExperiment
                times (:, 1) {mustBeNumeric}
                options.ShowConfidenceIntervals = false;
                options.Units string = "Wavelengths";
                options.dataType = 'TA';
            end
            
            scatterplots = gobjects(length(times), 1);
            legendlabels = strings(length(times), 1);
            colours = parula(length(times));

            figSpectra = figure();
            axSpectra = axes('Parent',figSpectra);
            hold(axSpectra, 'on');

            for i = 1:length(times)
                [spectrum, var] = obj.getWavelengthTrace(times(i), 'dataType', options.dataType);
                if options.Units == "Wavelengths"
                    scatterplots(i) = scatter(axSpectra, obj.wavelengths, spectrum, 10, colours(i, :), 'filled');
                    if options.ShowConfidenceIntervals
                        patch(axSpectra, [obj.wavelengths, flip(obj.wavelengths)], [spectrum + 2*(sqrt(var)), flip(spectrum - 2*(sqrt(var)))], 1, 'FaceColor', colours(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
                    end
                else
                    scatterplots(i) = scatter(axSpectra, obj.pixels, spectrum, 10, colours(i, :), 'filled');
                    if options.ShowConfidenceIntervals
                        patch(axSpectra, [obj.pixels, flip(obj.pixels)], [spectrum - 2*(sqrt(var)), flip(spectrum + 2*(sqrt(var)))], 1, 'FaceColor', colours(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
                    end
                end
                legendlabels(i) = strcat(sprintf('%.2f', times(i)), " / ps");
            end
            if options.Units == "Wavelengths"
                xlim(axSpectra, [min(obj.wavelengths), max(obj.wavelengths)]);
                xlabel(axSpectra, 'Wavelength / nm');
            else
                xlim(axSpectra, [min(obj.pixels), max(obj.pixels)]);
                xlabel(axSpectra, 'Pixels');
            end
            if options.dataType == "TA"
                ylabel(axSpectra, 'Transient Absorption / OD');
            else
                ylabel(axSpectra, 'Intensity');
            end

            set(axSpectra, 'FontSize', 12);
            set(axSpectra, 'Box', 'On', 'LineWidth', 2);
            set(axSpectra, 'TickLength', [0.015, 0.015]);
            set(axSpectra, 'Layer', 'top', 'XGrid', 'off', 'YGrid', 'off');
            title(axSpectra, obj.makeTitle());
            subtitle(axSpectra, obj.makeSubtitle());
            box(axSpectra, 'On');

            spectraLegend = legend(axSpectra, scatterplots, legendlabels);
            set(spectraLegend, 'FontSize', 12);

        end

        function [figTimeTraces, axTimeTraces] = peekTimeTraces(obj, xvalues, options)
            arguments
                obj TAExperiment
                xvalues (:, 1) {mustBeNumeric}
                options.ShowConfidenceIntervals = false;
                options.Units string = "Wavelengths";
                options.dataType = 'TA';
            end

            scatterplots = gobjects(length(xvalues), 1);
            legendlabels = strings(length(xvalues), 1);
            colours = parula(length(xvalues));

            figTimeTraces = figure();
            axTimeTraces = axes('Parent',figTimeTraces);
            hold(axTimeTraces, 'on');

            for i = 1:length(xvalues)
                [timeTrace, var] = obj.getTimeTrace(xvalues(i), 'Units', options.Units, 'dataType', options.dataType);
                scatterplots(i) = scatter(axTimeTraces, obj.times, timeTrace, 10, colours(i, :), 'filled');
                if options.ShowConfidenceIntervals
                    patch(axTimeTraces, [obj.times; flip(obj.times)], [timeTrace - 2*(sqrt(var)); flip(timeTrace + 2*(sqrt(var)))], 1, 'FaceColor', colours(i, :), 'EdgeColor', 'none', 'FaceAlpha', 0.5); 
                end
                
                if lower(options.Units) == "wavelengths"
                    unit = " / nm";
                else
                    unit = " / pixels";
                end
                legendlabels(i) = strcat(sprintf('%.2f', xvalues(i)), unit);
            end
            
            xlabel(axTimeTraces, 'Time / ps');
            if options.dataType == "TA"
                ylabel(axTimeTraces, 'Transient Absorption / OD');
            else
                ylabel(axTimeTraces, 'Intensity');
            end

            set(axTimeTraces, 'FontSize', 12);
            set(axTimeTraces, 'Box', 'On', 'LineWidth', 2);
            set(axTimeTraces, 'TickLength', [0.015, 0.015]);
            set(axTimeTraces, 'Layer', 'top', 'XGrid', 'off', 'YGrid', 'off');
            title(axTimeTraces, obj.makeTitle());
            subtitle(axTimeTraces, obj.makeSubtitle());
            box(axTimeTraces, 'On');

            spectraLegend = legend(axTimeTraces, scatterplots, legendlabels);
            set(spectraLegend, 'FontSize', 12);
        end

        function [timeTraceMean, timeTraceVar, timeTraceNshots] = getTimeTrace(obj, value, options)
            arguments
                obj TAExperiment
                value (1,:)
                options.Units = 'Wavelengths';
                options.dataType = 'TA';
            end

            switch options.dataType
                
                case 'TA'
                    datamean = obj.TAMean;
                    datavar = obj.TAVariance;
                    dataNShots = obj.TANShots;
                case 'Pump On'
                    datamean = obj.pumpOnMean;
                    datavar = obj.pumpOnVariance;
                    dataNShots = obj.pumpOnNShots;
                case 'Pump Off'
                    datamean = obj.pumpOffMean;
                    datavar = obj.pumpOffVariance;
                    dataNShots = obj.pumpOffNShots;
            end
            
            if length(value) == 1
                % Return single wavelength trace
                if lower(options.Units) == "wavelengths"
                    index = obj.findWavelength(value);
                else
                    index = obj.findPixel(value);
                end

                timeTraceMean = datamean(:, index);
                timeTraceVar = datavar(:, index);
                timeTraceNshots = dataNShots(:, index);
                
            elseif length(value) == 2
                if lower(options.Units) == "wavelengths"
                    indices = [obj.findWavelength(value(1,1)), obj.findWavelength(value(1,2))];
                else
                    indices = [obj.findPixel(value(1,1)), obj.findPixel(value(1,2))];
                end

                minIndex = min(indices);
                maxIndex = max(indices);

                timeTraceNshots = sum(dataNShots(:, minIndex:maxIndex), 2);
                timeTraceMean = sum(dataNShots(:, minIndex:maxIndex).*datamean(:, minIndex:maxIndex), 2, 'omitnan')./timeTraceNshots;
                timeTraceVar = sum(dataNShots(:, minIndex:maxIndex).*datavar(:, minIndex:maxIndex), 2, 'omitnan')./timeTraceNshots;
                
            else
                error('Wavelength argument must be a single wavelength or an array of two represneing upper and lower limits');
            end


        end

        function [wavelengthTraceMean, wavelengthTraceVar, wavelengthTraceNshots] = getWavelengthTrace(obj, time, options)
            arguments
                obj TAExperiment
                time (1,:)
                options.dataType = 'TA';
            end

            switch options.dataType
                % Use function for all data types default to TA
                case 'TA'
                    datamean = obj.TAMean;
                    datavar = obj.TAVariance;
                    dataNShots = obj.TANShots;
                case 'Pump On'
                    datamean = obj.pumpOnMean;
                    datavar = obj.pumpOnVariance;
                    dataNShots = obj.pumpOnNShots;
                case 'Pump Off'
                    datamean = obj.pumpOffMean;
                    datavar = obj.pumpOffVariance;
                    dataNShots = obj.pumpOffNShots;
            end
            
            if length(time) == 1
                % Return single wavelength trace
                index = obj.findTime(time);
                wavelengthTraceMean = datamean(index, :);
                wavelengthTraceVar = datavar(index, :);
                wavelengthTraceNshots = dataNShots(index, :);
                
            elseif length(time) == 2
                indices = [obj.findTime(time(1,1)), obj.findTime(time(1,2))];

                minIndex = min(indices);
                maxIndex = max(indices);

                wavelengthTraceNshots = sum(dataNShots(minIndex:maxIndex, :), 1);
                wavelengthTraceMean = sum(dataNShots(minIndex:maxIndex, :).*datamean(minIndex:maxIndex, :), 1, 'omitnan')./wavelengthTraceNshots;
                wavelengthTraceVar = sum(dataNShots(minIndex:maxIndex, :).*datavar(minIndex:maxIndex, :), 1, 'omitnan')./wavelengthTraceNshots;
                
            else
                error('Wavelength argument must be a single wavelength or an array of two represneing upper and lower limits');
            end


        end

        function applyExternalDispersionCorrection(obj, dispersionFitCoefficients)
            
            arguments
                obj TAExperiment
                dispersionFitCoefficients (1, :) {mustBeNumeric}
            end

            obj.dispersionFitCoefficients = dispersionFitCoefficients;
            obj.buildDispersionFit();
            obj.dispersionCorrectData();
            
        end

        function concExcitedMolecules = concExcitedMolecules(obj)
            concExcitedMolecules = (obj.numExcitedMolecules./6.02E23) ./ obj.pumpCylinderVolume();
        end

        function numExcitedMolecules = numExcitedMolecules(obj)
            if obj.opticalDensity > 0
                numExcitedMolecules = obj.numPumpPhotons - (obj.numPumpPhotons./10.^obj.opticalDensity);
            else
                error('Optical density value not valid');
            end
        end

        function numPhotons = numPumpPhotons(obj)
            if obj.pumpEnergy > 0
                numPhotons = (obj.pumpEnergy*1E-9)/obj.pumpPhotonEnergy;
            else
                error('Pump energy value not valid')
            end
                
        end

        function energy = pumpPhotonEnergy(obj)
            if obj.pumpWavelength > 0
                energy = (6.636E-34 * 2.998E8)./(obj.pumpWavelength*1E-9);
            else
                error('Pump wavelength value not valid');
            end
        end

        function volume = pumpCylinderVolume(obj)
            if obj.jetPathLength > 0
                pump_volume_m3 = obj.pumpArea() * obj.jetPathLength*1E-6; % Meters cubed
                volume = pump_volume_m3 * 1000; % volume in Litres
            else
                error('Jet path length value not valid')
            end
        end

        function area = pumpArea(obj)
            if obj.pumpSpotSize > 0
                area = 3.141*((obj.pumpSpotSize*1E-6)/2).^2; % Meters squared
            else
                error('Pump spot size value not valid')
            end
        end

        function updateCalibration(obj, calibrationPolynomial)
            obj.Calibration = calibrationPolynomial;
            obj.buildWavelengthVector();
        end
        
    end

    methods (Access = private)

        function buildDispersionFit(obj)
            obj.dispersionFit = zeros(size(obj.pixels));
            coeffs = flip(obj.dispersionFitCoefficients);
            
            for i = 1:length(coeffs)
                obj.dispersionFit = obj.dispersionFit + coeffs(i).*obj.pixels.^(i-1);
            end
            
        end

        function index = findWavelength(obj, wavelength)
            [~, index] = min(abs(obj.wavelengths - wavelength));
        end

        function index = findPixel(obj, pixel)
            [~, index] = min(abs(obj.pixels - pixel));
        end

        function index = findTime(obj, time)
            [~, index] = min(abs(obj.times - time));
        end

        function loadFile(obj)
            disp('Loading Burritos File...');

            data = jsondecode(fileread(obj.fileName));
            obj.json = data;
            
            obj.nPixels = data.Spectrometer.num_pixels;
            
            obj.buildPixelsVector();
            
            obj.Calibration(1, 1) = data.Spectrometer.wavelength_calibration.slope;
            obj.Calibration(1, 2) = data.Spectrometer.wavelength_calibration.intercept;

            obj.buildWavelengthVector;

            if isempty(data.Spectrometer.calibration_points)
                warning('No calibration points in data file.');
            else
                for i = 1:length(data.Spectrometer.calibration_points)
                    obj.CalibrationPoints = [obj.CalibrationPoints; CalibrationPoint(data.Spectrometer.calibration_points(i).wavelength, data.Spectrometer.calibration_points(i).pixel, data.Spectrometer.calibration_points(i).uncertainty)];
                end
            end

            obj.pumpWavelength = data.TransientAbsorption.pump_wavelength;
            obj.pumpBandwidth = data.TransientAbsorption.pump_bandwidth;
            obj.pumpEnergy = data.TransientAbsorption.pump_energy;
            obj.solvent = data.TransientAbsorption.solvent;
            try
                obj.opticalDensity = data.TransientAbsorption.optical_density;
            end

            try
                obj.jetPathLength = data.TransientAbsorption.jet_path_length;
            end

            try 
                obj.pumpAngle = data.TransientAbsorption.pump_angle;
            end

            try 
                obj.pumpSpotSize = data.TransientAbsorption.pump_spot_size;
            end

            try
                obj.probeSpotSize = data.TransientAbsorption.probe_spot_size;
            end

            obj.nTimes = length(data.TransientAbsorption.time_delays);
            obj.times = data.TransientAbsorption.time_delays;

            obj.nScans = length(data.TransientAbsorption.scans);
            obj.scans = TAScan.empty(obj.nScans, 0);

            for i = 1:obj.nScans
                obj.scans(i) = TAScan(obj.nTimes, obj.nPixels);
                obj.scans(i).populate714(data.TransientAbsorption.scans(i));
            end
            
            if ~isempty(data.TransientAbsorption.analytes)
                obj.analytes = Analyte.empty(length(data.TransientAbsorption.analytes), 0);
                for i = 1:length(data.TransientAbsorption.analytes)
                    mol = data.TransientAbsorption.analytes(i).analyte;
                    conc = data.TransientAbsorption.analytes(i).concentration;%sprintf('%.3f', data.TransientAbsorption.analytes(i).concentration);
                    obj.analytes(i) = Analyte(mol, conc);
                end
            end

            [obj.TAMean, obj.TAVariance, obj.TANShots, obj.pumpOnMean, obj.pumpOnVariance, obj.pumpOnNShots, obj.pumpOffMean, obj.pumpOffVariance, obj.pumpOffNShots] = deal(zeros(obj.nTimes, obj.nPixels));
            
            obj.combine_scans();

        end

        function loadFile716(obj)

            listTimes = [];

            lines = readlines(obj.fileName);
            
            obj.nPixels = 256;
            obj.buildPixelsVector();

            listNShots = [];
            listPumpOff = [];
            listPumpOn = [];
            listTAMean = [];
            listTAVar = [];

            for i = 1:length(lines)

                line = lines(i);

                switch line
                    case "# Spectrograph:  Slope and intercept"
                        slopeAndIntercept = split(lines(i+1));
                        obj.Calibration(1, 1) = str2double(slopeAndIntercept(1));
                        obj.Calibration(1, 2) = str2double(slopeAndIntercept(2));
                        obj.buildWavelengthVector();
                    case "# Pump-Probe Delay"
                        listTimes = [listTimes; str2double(lines(i+1))];
                        listNShots = [listNShots; str2double(split(lines(i+4)))'];
                        listPumpOff = [listPumpOff; str2double(split(lines(i+7)))'];
                        listPumpOn = [listPumpOn; str2double(split(lines(i+13)))'];
                        listTAMean = [listTAMean; 1E-3*str2double(split(lines(i+19)))'];
                        listTAVar = [listTAVar; (1E-3*str2double(split(lines(i+22)))).^2'];
                end

                
            end
            
            obj.times = unique(listTimes);
            obj.nTimes = length(obj.times);

            obj.nScans = ceil(length(listTimes) / length(obj.times));

            pad = zeros(obj.nScans*length(obj.times) - length(listTimes), obj.nPixels);

            listNShots = [listNShots; pad];
            listPumpOff = [listPumpOff; pad];
            listPumpOn = [listPumpOn; pad];
            listTAMean = [listTAMean; pad];
            listTAVar = [listTAVar; pad];

            cubeNShots = reshape(listNShots, obj.nTimes, obj.nScans, obj.nPixels);
            cubePumpOff = reshape(listPumpOff, obj.nTimes, obj.nScans, obj.nPixels);
            cubePumpOn = reshape(listPumpOn, obj.nTimes, obj.nScans, obj.nPixels);
            cubeTAMean = reshape(listTAMean, obj.nTimes, obj.nScans, obj.nPixels);
            cubeTAVar = reshape(listTAVar, obj.nTimes, obj.nScans, obj.nPixels);
            
            obj.scans = TAScan.empty(obj.nScans, 0);
            
            for i = 1:obj.nScans
                obj.scans(i) = TAScan(obj.nTimes, obj.nPixels);
                obj.scans(i).populate716(cubeTAMean(:, i, :), cubeTAVar(:, i, :), cubeNShots(:, i, :), cubePumpOff(:, i, :), cubePumpOn(:, i, :));
                
            end

            [obj.TAMean, obj.TAVariance, obj.TANShots, obj.pumpOnMean, obj.pumpOnVariance, obj.pumpOnNShots, obj.pumpOffMean, obj.pumpOffVariance, obj.pumpOffNShots] = deal(zeros(obj.nTimes, obj.nPixels));
            
            obj.combine_scans();

        end

        function buildPixelsVector(obj)
            obj.pixels = linspace(1, obj.nPixels, obj.nPixels);
        end

        function buildWavelengthVector(obj)
            obj.wavelengths = obj.pixels * obj.Calibration(1,1) + obj.Calibration(1, 2);
        end

        function combine_scans(obj)
            % Compute a weightd mean of the data from the different scans.
            % Check both the mean and variance in each scan for nan values
            % and exclude that point in that scan for the averaging
            % Pooled variance = weighted average of variance

            [TAMeanCube, TAVarCube, TANShotsCube, pumpOnMeanCube, pumpOnVarCube, pumpOnNShotsCube, pumpOffMeanCube, pumpOffVarCube, pumpOffNShotsCube] = deal(zeros(obj.nScans, obj.nTimes, obj.nPixels));

            TANShotsCube(isnan(TAMeanCube + TAVarCube)) = 0;
            pumpOnNShotsCube(isnan(pumpOnMeanCube + pumpOnVarCube)) = 0;
            pumpOffNShotsCube(isnan(pumpOffMeanCube + pumpOffVarCube)) = 0;

            for i = 1:obj.nScans
                TAMeanCube(i, :, :) = obj.scans(i).TAMean;
                TAVarCube(i, :, :) = obj.scans(i).TAVariance;
                TANShotsCube(i, :, :) = obj.scans(i).TANShots;
                pumpOnMeanCube(i, :, :) = obj.scans(i).pumpOnMean;
                pumpOnVarCube(i, :, :) = obj.scans(i).pumpOnVariance;
                pumpOnNShotsCube(i, :, :) = obj.scans(i).pumpOnNShots;
                pumpOffMeanCube(i, :, :) = obj.scans(i).pumpOffMean;
                pumpOffVarCube(i, :, :) = obj.scans(i).pumpOffVariance;
                pumpOffNShotsCube(i, :, :) = obj.scans(i).pumpOffNShots;
            end

            obj.TAMean = reshape(sum(TANShotsCube.*TAMeanCube, 1, 'omitnan')./sum(TANShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.TAVariance = reshape(sum(TANShotsCube.*TAVarCube, 1, 'omitnan')./sum(TANShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.TANShots = reshape(sum(TANShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOnMean = reshape(sum(pumpOnNShotsCube.*pumpOnMeanCube, 1, 'omitnan')./sum(pumpOnNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOnVariance = reshape(sum(pumpOnNShotsCube.*pumpOnVarCube, 1, 'omitnan')./sum(pumpOnNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOnNShots = reshape(sum(pumpOnNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOffMean = reshape(sum(pumpOffNShotsCube.*pumpOffMeanCube, 1, 'omitnan')./sum(pumpOffNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOffVariance = reshape(sum(pumpOffNShotsCube.*pumpOffVarCube, 1, 'omitnan')./sum(pumpOffNShotsCube, 1), obj.nTimes, obj.nPixels);
            obj.pumpOffNShots = reshape(sum(pumpOffNShotsCube, 1), obj.nTimes, obj.nPixels);

            % Perform a final check for any points which were a nan across
            % al scans (typical of a portion of the spectrometer upon which
            % no white light falls and set them to zero for future
            % algorithms.
            %obj.TAMean(isnan(obj.TAMean)) = 0;
            obj.TAVariance(isnan(obj.TAVariance)) = 0;
            obj.pumpOnMean(isnan(obj.pumpOnMean)) = 0;
            obj.pumpOnVariance(isnan(obj.pumpOnVariance)) = 0;
            obj.pumpOffMean(isnan(obj.pumpOffMean)) = 0;
            obj.pumpOffVariance(isnan(obj.pumpOffVariance)) = 0;

        end

        function titlestring = makeTitle(obj)
            titlestring = "";
            if ~isempty(obj.analytes)
                for i = 1:length(obj.analytes)
                    analytestring = strcat(obj.analytes.name, " (", sprintf('%.2f', obj.analytes.concentration), " \muM)  ");
                    titlestring = strcat(titlestring, analytestring);
                end
            else
                titlestring = "Analyte Unknown";
            end

            titlestring = strcat(titlestring, '/ ', obj.solvent);

        end

        function subtitlestring = makeSubtitle(obj)

            wavelength = sprintf('%.0f', obj.pumpWavelength);
            bandwidth = sprintf('%.0f', obj.pumpBandwidth);
            energy = sprintf('%.1f', obj.pumpEnergy);
            subtitlestring = strcat("\lambda_{exc} = ", wavelength, " nm \Delta_{\lambda} = ", bandwidth, " nm E = ", energy, " nJ");
            
        end

        function dispersionCorrectData(obj)

            [~, timeGrid] = meshgrid(obj.pixels, obj.times);
            time_shifted_grid = timeGrid - obj.dispersionFit; % Care must be taken that the dispersion curve is a row vec while the timegrid has columns as pixels amnd the rows as times
            
            %for i = 1:obj.nScans
            %    for j = 1:obj.nPixels
            %        obj.scans(i).TAMean(:, j) = interp1(time_shifted_grid(:, j), obj.scans(i).TAMean(:, j), obj.times, 'linear', 'extrap');
            %    end
            %end

            for i = 1:obj.nPixels
                   obj.TAMean(:, i) = interp1(time_shifted_grid(:, i), obj.TAMean(:, i), obj.times, 'linear', 'extrap');
            end

            %obj.combine_scans();

        end

    end

end
