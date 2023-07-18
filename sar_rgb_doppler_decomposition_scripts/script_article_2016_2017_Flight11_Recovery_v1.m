% Display Recovery Ice Stream SAR images and SAR-derived products
%
%
%
%

fileInSingleCell = {...
                    '...\2016_2017_Flight11_SAR_ktrace33to35_TxPortC13L4PortJ13L4_RxSC_RangeCompression.nc',...
                    };
tracesPRI = 0.2; % (s)
% Labels
outputNamesCell = {'recSARimage_HV{1}'};
titleCell = cell(1,length(outputNamesCell));
% Image thresholds
dBLim = [0 60]; % (dB), maximum lower limit below the maximum
% Palette parameters for universal readability
palette_codes = [-1 0 -2 -3 3];
blindness_type = 'deuteranopia'; % 'protanopia', 'deuteranopia', 'tritanopia' and 'achromatopsia'
% Display parameters
aslFlag = true; % 0, thickness; 1, height above sea level
if ~aslFlag % thickness
    vertical_limits_shown = [-180 1600]/1e3; % [km]
else
    vertical_limits_shown = [-1500 250]/1e3; % [km]
end
% along_track_limits = [3300 3600]/1e3; % [kilotraces]
along_track_limits = [3340 3520]/1e3; % [kilotraces]
decomposition_band_from_aperture_flag = false;
decomposition_processed_band = 50; % [Hz]
% Change how to display the along-track axes
display_track_distance_flag = true; % false: time; true: distance
number_display_intervals = 4;
precission_m = 10; % [m]
% Image mosaicking type
mosaic_type = 2; % 0: no mosaicking; 1: 2x2; 2: 3x1
% Add box labels in the images of the mosaic
add_box_labels_in_mosaic_flag = true; % false: not to add; true: to add
% Figure sizes
figure_units = 'normalized';
figure_position_norm = [0    0.0581    1    0.8686]; % normalised units
% Figures to be displayed
figure_objects = gobjects(0);
SAR_full_basic_disp_flag = true;
RGB_full_basic_disp_flag = true;
RGB_unique_basic_disp_flag = true;
RGB_max_basic_disp_flag = false;
RGB_full_universal_disp_flag = true;
RGB_unique_universal_disp_flag = true;
RGB_max_universal_disp_flag = false;
RGB_over_full_universal_disp_flag = false;
RGB_over_unique_universal_disp_flag = false;
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% Read SAR images and display basic representations.
%
% The SAR image is represented in grey scale.
% The RGB Doppler-Decompositions are represented with the pure Red, Green
%   and Blue colours.
%
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% Get maximum of each SAR image
for ff = 1:length(fileInSingleCell)
    out_file = fileInSingleCell{ff};
    %--------------------------------
    % Retrieve data from netCDF file
    %--------------------------------
    % Open file
    ncID = netcdf.open(out_file);
    disp(['out_file (netCDF) = ' out_file])
    % Retrieve 'intervalAbsoluteSAR' (-)
    recVarName = 'intervalAbsoluteSAR';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recIntervalAbsoluteSAR = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (-)
    recMaxAbsoluteSAR_pol(ff) = recIntervalAbsoluteSAR(2);
    % Close file
    netcdf.close(ncID);
end
polarimetricCal = max(recMaxAbsoluteSAR_pol); % (-)
% Get and represent SAR images and RGB Doppler decompositions
for ff = 1:length(fileInSingleCell)
    out_file = fileInSingleCell{ff};
    %--------------------------------
    % Retrieve data from netCDF file
    %--------------------------------
    % Open file
    ncID = netcdf.open(out_file);
    disp(['out_file (netCDF) = ' out_file])
    %--------------------------------
    % Global attributes
    %--------------------------------
    [ndims,nvars,ngatts,unlimdimid] = netcdf.inq(ncID);
    fprintf(' global attributes:');
    maxAttlen = 320; % maximum lenght of chars to show
    maxDispLen = 80;
    for k = 0:ngatts-1
        % Read parameters
        attname = netcdf.inqAttName(ncID, netcdf.getConstant('NC_GLOBAL'),k);
        [xtype,attlen] = netcdf.inqAtt(ncID, netcdf.getConstant('NC_GLOBAL'), attname);
        gattval = netcdf.getAtt(ncID,netcdf.getConstant('NC_GLOBAL'),attname);
        % Limit length
        gattval = gattval(1:min(maxAttlen,attlen));
        % Replace 'newline' and 'carriage return' character
        gattval = strrep(gattval, newline, '');
        % Replace 'carriage return' character
        gattval = strrep(gattval, char(13), '');
        attlen = length(gattval);
        % Display
        fprintf('\n\t''%s'': %s', attname, gattval(1:min(length(gattval),maxDispLen)));
        currentAttLenPos = min(maxAttlen,maxDispLen);
        while currentAttLenPos < length(gattval)
            fprintf('...')
            lastCurrentPosition = min(currentAttLenPos+maxDispLen, length(gattval));
            fprintf('\n\t\t %s', gattval(currentAttLenPos+1 : lastCurrentPosition));
            currentAttLenPos = lastCurrentPosition;
        end
    end
    fprintf(newline);
    %--------------------------------
    % Variables
    %--------------------------------
    % Retrieve 'alongTrackCentralAngle' (deg)
    recAlongTrackCentralAngle = netcdf.getVar(ncID, netcdf.inqVarID(ncID,'alongTrackCentralAngle')); % (deg)
    % Retrieve 'alongTrackApertureAngle' (deg)
    recAlongTrackApertureAngle = netcdf.getVar(ncID, netcdf.inqVarID(ncID,'alongTrackApertureAngle')); % (deg)
    % Retrieve 'pulse'
    recVarName = 'pulse';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recPulse = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths);
    recPulse_fillValue = netcdf.getAtt(ncID, netcdf.inqVarID(ncID,recVarName), '_FillValue');
    recValidPulses = recPulse ~= recPulse_fillValue;
    recPulse = recPulse(recValidPulses);
    % Retrieve 'pulseTime' (according to 'timeFormatFlag')
    recVarName = 'pulseTime';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recPulseTime = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (according to 'timeFormatFlag')
    recPulseTime = recPulseTime(recValidPulses); % (according to 'timeFormatFlag')
    % Retrieve 'SARimageReal'
    recVarName = 'SARimageReal';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recSARimage = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); %
    recSARimage = recSARimage(:,recValidPulses); %
    % Retrieve 'SARimageImag'
    recVarName = 'SARimageImag';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recSARimageImag = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); %
    recSARimageImag = recSARimageImag(:,recValidPulses); %
    recSARimage = recSARimage + 1j*recSARimageImag;
    clear recSARimageImag;
    % Retrieve 'intervalAbsoluteSAR' (-)
    recVarName = 'intervalAbsoluteSAR';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recIntervalAbsoluteSAR = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (-)
    recMaxAbsoluteSAR = recIntervalAbsoluteSAR(2);
    % Retrieve 'radialDepthGrid' (m)
    recVarName = 'radialDepthGrid';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recRadialDepthGrid = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (m)
    % retrieve 'heightAircraft' (m)
    recVarName = 'heightAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recHeightAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recHeightAircraft = recHeightAircraft(recValidPulses); % (deg)
    % Retrieve 'terrainClearanceAircraft' (m)
    recVarName = 'terrainClearanceAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recTerrainClearanceAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recTerrainClearanceAircraft = recTerrainClearanceAircraft(recValidPulses); % (deg)
    % Retrieve 'speedAircraft' (m/s)
    recVarName = 'speedAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recSpeedAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recSpeedAircraft = recSpeedAircraft(recValidPulses); % (deg)
    % Retrieve 'latitudeAircraft' (deg)
    recVarName = 'latitudeAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recLatitudeAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recLatitudeAircraft = recLatitudeAircraft(recValidPulses); % (deg)
    % Retrieve 'longitudeAircraft' (deg)
    recVarName = 'longitudeAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recLongitudeAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recLongitudeAircraft = recLongitudeAircraft(recValidPulses); % (deg)
    % Retrieve 'rollAircraft' (deg)
    recVarName = 'rollAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recRollAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recRollAircraft = recRollAircraft(recValidPulses); % (deg)
    % Retrieve 'pitchAircraft' (deg)
    recVarName = 'pitchAircraft';
    [varname,xtype,dimids,natts] = netcdf.inqVar(ncID,netcdf.inqVarID(ncID,recVarName));
    fullLengths = zeros(1,length(dimids));
    for kd = 1:length(dimids)
        [dimname, dimLength] = netcdf.inqDim(ncID, dimids(kd)); % unlimdimid = dimids(end)
        fullLengths(kd) = dimLength;
    end
    recPitchAircraft = netcdf.getVar(ncID, netcdf.inqVarID(ncID,recVarName), zeros(1,length(dimids)), fullLengths); % (deg)
    recPitchAircraft = recPitchAircraft(recValidPulses); % (deg)
    % Retrieve 'flight' global attribute
    recFlight = netcdf.getAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'flight');
    % Retrieve 'transmittersAndWaveforms' global attribute
    recTransmittersAndWaveformsString = netcdf.getAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'transmittersAndWaveforms');
    % Retrieve 'receivers' global attribute
    recReceiversString = netcdf.getAtt(ncID, netcdf.getConstant('NC_GLOBAL'), 'receivers');
    % Close file
    netcdf.close(ncID);
    %--------------------------------
    % strings for titles
    %--------------------------------
    % Transmitters label in format 'TxS9ABC_C13L4; TxS9ABC_J13L4; ...'
    TXWVstring = recTransmittersAndWaveformsString;
    TXWVstring = strrep(TXWVstring, '_', '');
    TXWVstring = strrep(TXWVstring, ',', '-');
    % Replace P1234 with Port
    TXWVstring = strrep(TXWVstring,'P1234','Port');
    % Replace S9ABC with Star
    TXWVstring = strrep(TXWVstring,'S9ABC','Star');
    % Receivers label in format 'RxP1; RxP2; ...'
    RXstring = recReceiversString;
    RXstring = strrep(RXstring, 'Rx', '');
    RXstring = strrep(RXstring, ',', '');
    RXstring = ['Rx' sprintf('%s',RXstring)];
    % Replace P1P2P3P4 with Port
    RXstring = strrep(RXstring,'P1P2P3P4','Port');
    % Replace P1P2P3P4 with Belly
    RXstring = strrep(RXstring,'B5B6B7B8','Belly');
    % Replace P1P2P3P4 with Star
    RXstring = strrep(RXstring,'S9SASBSC','Star');
    % Replace PortBellyStar with All
    RXstring = strrep(RXstring,'PortBellyStar','All');
    %--------------------------------
    % Reformatting
    %--------------------------------
    recPulseDouble = double(recPulse);
    %--------------------------------
    % Calculated parameters
    %--------------------------------
    recPulseDiff = recPulseDouble(2) - recPulseDouble(1);
    foundPRF = recPulseDiff/(mean(diff(recPulseTime))*24*3600); % (Hz)
    foundResDepth = mean(diff(recRadialDepthGrid)); % (m/sample)
    foundHeightSurface = recHeightAircraft - recTerrainClearanceAircraft; %(m)
    traces = recPulseDouble/(foundPRF*tracesPRI);
    track_distance = [0; cumsum(0.5*(recSpeedAircraft(1:end-1) + recSpeedAircraft(2:end)))./(foundPRF./diff(recPulseDouble))]; % (m)
    useTracesAxis = true; % 0, use radar pulses; 1, use traces
    if useTracesAxis
        alongTrackAxis = traces/1e3; % (kilotraces)
        xLabelString = ['Along-track (x10^3 ' num2str(tracesPRI*1e3) ' ms traces)'];
    else
        alongTrackAxis = recPulseDouble;
        xLabelString = ['Along-track (radar data, prf ' num2str(round(foundPRF*100)/100) 'Hz)'];
    end
    recMaxAbsoluteSAR = polarimetricCal;
    %--------------------------------
    % RGB Doppler decomposition
    %--------------------------------
    if decomposition_band_from_aperture_flag
        approxApBandProcessed = recAlongTrackApertureAngle; % (Hz), in PASIN2 the aperture almost equals the Doppler band
    else
        approxApBandProcessed = decomposition_processed_band; % [Hz]
    end
    ratioAzDop = foundPRF/recPulseDiff/(approxApBandProcessed/3); % >1, of each doppler bandwidth over the full processing band.
    percentageShift = 100*(approxApBandProcessed/3)/(foundPRF/recPulseDiff); % shift of each central frequency, in percentage over the full processing band.
    equalization_flag = true;
    azFilteredRGB = rgb_doppler_decomposition_dB_min_max_v1(recSARimage/recMaxAbsoluteSAR, ratioAzDop, percentageShift, recPulseDouble, recRadialDepthGrid/1e3, dBLim, equalization_flag);
    close(gcf); close(gcf); close(gcf) % 8,32; 25,15;
    numberOfColors = 64;
    RGBFreqLimits = [-percentageShift+(100/ratioAzDop)*[-0.5 0.5] (100/ratioAzDop)*[-0.5 0.5] percentageShift+(100/ratioAzDop)*[-0.5 0.5]];
    RGBFreqLimitsNorm64 = round((RGBFreqLimits-min(RGBFreqLimits))/(max(RGBFreqLimits)-min(RGBFreqLimits))*(numberOfColors-1)+1);
    RGBcolorbar = zeros(numberOfColors,3);
    RGBcolorbar(RGBFreqLimitsNorm64(1):RGBFreqLimitsNorm64(2),:) = repmat([1 0 0],[RGBFreqLimitsNorm64(2)-RGBFreqLimitsNorm64(1)+1 1]);
    RGBcolorbar(RGBFreqLimitsNorm64(3):RGBFreqLimitsNorm64(4),:) = repmat([0 1 0],[RGBFreqLimitsNorm64(4)-RGBFreqLimitsNorm64(3)+1 1]);
    RGBcolorbar(RGBFreqLimitsNorm64(5):RGBFreqLimitsNorm64(6),:) = repmat([0 0 1],[RGBFreqLimitsNorm64(6)-RGBFreqLimitsNorm64(5)+1 1]);
    %--------------------------------
    % Figures
    %--------------------------------
    % Above sea level
    if ~aslFlag % thickness
        % SAR image
        if SAR_full_basic_disp_flag
            current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
            figure_objects = [figure_objects current_figure];
            imagesc(alongTrackAxis, recRadialDepthGrid/1e3, 20*log10(abs(recSARimage)/recMaxAbsoluteSAR)), grid on, colorbar, caxis(-dBLim(end:-1:1)),
            set(gca,'FontSize',18), xlabel(xLabelString), ylabel('Depth (km)'),
            title(['SAR image, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
            colormap('gray'), auxCBar=colorbar('FontSize',18); auxCBarLabel = auxCBar.Label;
            set(auxCBarLabel,'String','Amplitude(dB)','Rotation',-90,'VerticalAlignment','bottom'); clear auxCBar,
        end
        % RGB Doppler decomposition
        if RGB_full_basic_disp_flag
            current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
            figure_objects = [figure_objects current_figure];
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGB),grid on,set(gca,'FontSize',18),
            xlabel(xLabelString), ylabel('Depth (km)'),
            title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
            colormap(RGBcolorbar),
            auxCBar=colorbar('FontSize',18,'Direction','Reverse');
            auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
            set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
        end
    else % Height above sea level
        foundHeightSurface_samples = round(foundHeightSurface/foundResDepth); % (samples)
        foundHeightSurface_shiftSamples = max(foundHeightSurface_samples) - foundHeightSurface_samples; % (samples)
        recSARimage_SL = recSARimage;
        for azkk = 1:size(recSARimage_SL,2)
            recSARimage_SL(:,azkk) = circshift(recSARimage(:,azkk),[foundHeightSurface_shiftSamples(azkk) 0]);
            azFilteredRGB(:,azkk,:) = circshift(azFilteredRGB(:,azkk,:),[foundHeightSurface_shiftSamples(azkk) 0 0]);
        end
        heightASLrecRadialDepthGrid = -recRadialDepthGrid+max(foundHeightSurface_samples)*foundResDepth; % (m)
        % SAR image
        if SAR_full_basic_disp_flag
            current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
            figure_objects = [figure_objects current_figure];
            imagesc(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, 20*log10(abs(recSARimage_SL)/recMaxAbsoluteSAR)), axis xy, grid on, colorbar, caxis(-dBLim(end:-1:1)),
            set(gca,'FontSize',18), xlabel(xLabelString), ylabel('Height a.s.l. (km)'),
            title(['SAR image, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
            colormap('gray'), auxCBar=colorbar('FontSize',18); auxCBarLabel = auxCBar.Label;
            set(auxCBarLabel,'String','Amplitude(dB)','Rotation',-90,'VerticalAlignment','bottom'); clear auxCBar,
        end
        % RGB Doppler decomposition
        if RGB_full_basic_disp_flag
            current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
            figure_objects = [figure_objects current_figure];
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGB),axis xy,grid on,set(gca,'FontSize',18),
            xlabel(xLabelString), ylabel('Height a.s.l. (km)'),
            title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
            colormap(RGBcolorbar),
            auxCBar=colorbar('FontSize',18,'Direction','Reverse');
            auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
            set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
        end
    end
    eval([outputNamesCell{ff} ' = recSARimage;'])
    titleCell{ff} = [TXWVstring(1:6) RXstring];
end
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% RGB Doppler-Decompositions representations with non-overlapping bands.
%   Displays:
%   - image with the intensity of the most intense band in each SAR sample
%   - image with the the most intense band in each SAR sample
%
% The RGB Doppler-Decompositions are represented with the pure Red, Green
%   and Blue colours.
%
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%------------------------------------------------
% Image with the intensity of the most intense band in each SAR sample
%------------------------------------------------
NN = size(azFilteredRGB,1);
MM = size(azFilteredRGB,2);
azFilteredRGBunique = zeros(NN, MM, 3, 'uint8');
% Get most intense band
[~,indColourMax] = max(azFilteredRGB,[],3); 
% Index in 3D matrix
indexMax = (1:NN*MM).'+(indColourMax(:)-1)*NN*MM;
% Assign intensity
azFilteredRGBunique(indexMax) = azFilteredRGB(indexMax);
% Indexed image
indexed_azFilteredRGB = zeros(NN, MM, 'uint8');
wanted_colours_per_band = 8;
for index_rgb = 1:3
    index_colour_offset = (index_rgb-1)*wanted_colours_per_band + 1;
    indexed_azFilteredRGB = indexed_azFilteredRGB + ...
        uint8(...
            index_colour_offset*(azFilteredRGBunique(:,:,index_rgb) > 0) + ...
            floor(double(azFilteredRGBunique(:,:,index_rgb))*(wanted_colours_per_band/255))...
        );
end
% Display image with the intensity of the most intense band in each SAR sample
if RGB_unique_basic_disp_flag
    current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
    figure_objects = [figure_objects current_figure];
    if ~aslFlag % thickness
        image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBunique),grid on,
        ylabel('Depth (km)')
    else
        image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBunique),grid on,axis xy
        ylabel('Height a.s.l. (km)'),
    end
    set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
    xlabel(xLabelString),
    title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
    colormap(RGBcolorbar),
    auxCBar=colorbar('FontSize',18,'Direction','Reverse');
    auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
    set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
end
% Display image with the the most intense band in each SAR sample
if RGB_max_basic_disp_flag
    current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
    figure_objects = [figure_objects current_figure];
    if ~aslFlag % thickness
        image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBunique*255),grid on,
        ylabel('Depth (km)')
    else
        image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBunique*255),grid on,axis xy
        ylabel('Height a.s.l. (km)'),
    end
    set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
    xlabel(xLabelString),
    title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
    colormap(RGBcolorbar),
    auxCBar=colorbar('FontSize',18,'Direction','Reverse');
    auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
    set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
end
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% RGB Doppler-Decompositions representations with non-overlapping bands.
%   Displays:
%   - image with RGB Doppler-Decomposition
%   - image with the intensity of the most intense band in each SAR sample
%   - image with the the most intense band in each SAR sample
%
% The RGB Doppler-Decompositions are represented with other than the pure 
%   Red, Green and Blue colours, depending on the parameters chosen in the
%   function 'universally_readable_colourmap.m', according to the integer
%   parameter 'palette_code':
%   
%   Primary colours
%       red = RGB(255,0,0)
%       green = RGB(0,255,0)
%       blue = RGB(0,0,255)
%   Secondary colours
%       olive = RGB(128,128,0)
%       teal = RGB(0,128,128)
%       purple = RGB(128,0,128)
%
%   palette_code (integer):
%          -3: Interpolation of secondary (olive, teal, purple) and
%              modified-primary (red, green, blue) colours, with the total
%              summation of colours equal to pure white ([1,1,1]). The
%              weight of the secondary colours is 75% and the weight of the
%              primary 25%. These weights are given by the fixed parameter:
%              secondary_weight = 0.75;
%          -2: Secondary colours (olive, teal, purple), with the total
%              summation of colours equal to pure white ([1,1,1]).
%          -1: Modified Red, Green, Blue for uniform perception, still with
%              the total summation of colours equal to pure white ([1,1,1]).
%              Reference: Kovesi, P. Good colour maps: how to design them. CoRR abs/1509.03700 (2015).
%           0: Pure Red, Green, Blue
%           1: Universally readable, v1:
%              [goldenish, light brownish, violetish-bluish], the total
%              summation of colours is not equal to pure white.
%           2: Universally readable, v2:
%              [goldenish, light brownish, violetish-bluish], the total
%              summation of colours is not equal to pure white.
%           3: Universally readable, v3:
%              [yellowish, dark_grey, bluish], the total
%              summation of colours equals to pure white ([1,1,1]).
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% Palette parameters
blind_conversion_method = 4;
wanted_colours_per_band = 8;
kovesi_modification_flag = false;
display_single_palette_flag = false;
display_multiple_palette_flag = false;
% Loop all palette codes
for palette_code = palette_codes
    [readable_colourmap,...
    first_band_template, second_band_template, third_band_template] = universally_readable_colourmap(...
        palette_code,...
        blind_conversion_method,...
        wanted_colours_per_band,...
        kovesi_modification_flag,...
        display_single_palette_flag,...
        display_multiple_palette_flag);
    % Colours representing each of the three Frequency sub-bands
    colour_1st_band = first_band_template(1,:); % [0...255]
    colour_2nd_band = second_band_template(1,:); % [0...255]
    colour_3rd_band = third_band_template(1,:); % [0...255]
    colour_1st_band = rgb_to_colour_blindness(colour_1st_band, blindness_type, blind_conversion_method); % [0...255]
    colour_2nd_band = rgb_to_colour_blindness(colour_2nd_band, blindness_type, blind_conversion_method); % [0...255]
    colour_3rd_band = rgb_to_colour_blindness(colour_3rd_band, blindness_type, blind_conversion_method); % [0...255]
    colour_1st_band = colour_1st_band/255; % [0...1]
    colour_2nd_band = colour_2nd_band/255; % [0...1]
    colour_3rd_band = colour_3rd_band/255; % [0...1]
    % Build colormap for colorbar
    modRGBcolorbar = zeros(size(RGBcolorbar), class(RGBcolorbar));
    for idx_band = 1:3
        % Modified RGB
        modRGBcolorbar(:,idx_band) = ...
            colour_1st_band(idx_band)*RGBcolorbar(:,1) + ...
            colour_2nd_band(idx_band)*RGBcolorbar(:,2) + ...
            colour_3rd_band(idx_band)*RGBcolorbar(:,3); % [0...255]
    end
    % Calculate images
    azFilteredRGBmod = zeros(size(azFilteredRGBunique), 'uint8');
    azFilteredRGBmod_unique = zeros(size(azFilteredRGBunique), 'uint8');
    azFilteredRGBmod_unique_max = zeros(size(azFilteredRGBunique), 'uint8');
    for idx_band = 1:3
        % Modified RGB
        azFilteredRGBmod(:,:,idx_band) = ...
            colour_1st_band(idx_band)*azFilteredRGB(:,:,1) + ...
            colour_2nd_band(idx_band)*azFilteredRGB(:,:,2) + ...
            colour_3rd_band(idx_band)*azFilteredRGB(:,:,3); % [0...255]
        % Intensity of most intense band
        azFilteredRGBmod_unique(:,:,idx_band) = ...
            colour_1st_band(idx_band)*azFilteredRGBunique(:,:,1) + ...
            colour_2nd_band(idx_band)*azFilteredRGBunique(:,:,2) + ...
            colour_3rd_band(idx_band)*azFilteredRGBunique(:,:,3); % [0...255]
        % Maximum intensity of most intense band
        azFilteredRGBmod_unique_max(:,:,idx_band) = ...
            colour_1st_band(idx_band)*(azFilteredRGBunique(:,:,1)*255) + ...
            colour_2nd_band(idx_band)*(azFilteredRGBunique(:,:,2)*255) + ...
            colour_3rd_band(idx_band)*(azFilteredRGBunique(:,:,3)*255); % [0...255]
    end
    azFilteredRGBmod_unique = uint8(azFilteredRGBmod_unique);
    % Modified RGB
    if RGB_full_universal_disp_flag
        current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
        figure_objects = [figure_objects current_figure];
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod),grid on,
            ylabel('Depth (km)')
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod),grid on,axis xy
            ylabel('Height a.s.l. (km)'),
        end
        set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
        xlabel(xLabelString),
        title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
        colormap(modRGBcolorbar),
        auxCBar=colorbar('FontSize',18,'Direction','Reverse');
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
    end
    % Display image with the intensity of the most intense band in each SAR sample
    if RGB_unique_universal_disp_flag
        current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
        figure_objects = [figure_objects current_figure];
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod_unique),grid on,
            ylabel('Depth (km)')
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod_unique),grid on,axis xy
            ylabel('Height a.s.l. (km)'),
        end
        set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
        xlabel(xLabelString),
        title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
        colormap(modRGBcolorbar),
        auxCBar=colorbar('FontSize',18,'Direction','Reverse');
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
    end
    % Display image with the the most intense band in each SAR sample
    if RGB_max_universal_disp_flag
        current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
        figure_objects = [figure_objects current_figure];
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod_unique_max),grid on,
            ylabel('Depth (km)')
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod_unique_max),grid on,axis xy
            ylabel('Height a.s.l. (km)'),
        end
        set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
        xlabel(xLabelString),
        title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
        colormap(modRGBcolorbar),
        auxCBar=colorbar('FontSize',18,'Direction','Reverse');
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
    end
end
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% RGB Doppler-Decompositions representations with overlapping bands.
%   Displays:
%   - image with RGB Doppler-Decomposition
%   - image with the intensity of the most intense band in each SAR sample
%
% The RGB Doppler-Decompositions are represented with other than the pure 
%   Red, Green and Blue colours, depending on the parameters chosen in the
%   function 'universally_readable_colourmap.m', according to the integer
%   parameter 'palette_code':
%
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
percentageShift_over = (approxApBandProcessed/(foundPRF/recPulseDiff))*100*[-0.3 0 0.3];
ratioAzDop_over = (foundPRF/recPulseDiff)./(approxApBandProcessed*[0.4 0.6 0.4]);
azFilteredRGB_over = rgb_doppler_decomposition_dB_min_max_v1(recSARimage/recMaxAbsoluteSAR, ratioAzDop_over, percentageShift_over, recPulseDouble, recRadialDepthGrid/1e3, dBLim, equalization_flag);
if aslFlag % thickness
    for azkk = 1:size(recSARimage,2)
        azFilteredRGB_over(:,azkk,:) = circshift(azFilteredRGB_over(:,azkk,:),[foundHeightSurface_shiftSamples(azkk) 0 0]);
    end
end
close(gcf); close(gcf); close(gcf) % 8,32; 25,15;
numberOfColors_over = 256;
RGBFreqLimits_over = repmat(percentageShift_over(:).', 2, 1) + [-0.5; 0.5]*(100./ratioAzDop_over(:).');
RGBFreqLimitsNorm = (RGBFreqLimits_over-min(RGBFreqLimits_over(:)))/(max(RGBFreqLimits_over(:))-min(RGBFreqLimits_over(:)))*(numberOfColors_over-1)+1;
RGBFreqLimitsNorm_int_over = round(RGBFreqLimitsNorm);
%------------------------------------------------
% Image with the intensity of the most intense band in each SAR sample
%------------------------------------------------
NN = size(azFilteredRGB_over,1);
MM = size(azFilteredRGB_over,2);
azFilteredRGBunique_over = zeros(NN, MM, 3, 'uint8');
% Get most intense band
[~,indColourMax] = max(azFilteredRGB_over,[],3); 
% Index in 3D matrix
indexMax = (1:NN*MM).'+(indColourMax(:)-1)*NN*MM;
% Assign intensity
azFilteredRGBunique_over(indexMax) = azFilteredRGB_over(indexMax);
% Loop all palette codes
for palette_code = palette_codes
    % Initialize colourmap
    RGBcolorbar_over = zeros(numberOfColors_over, 3);
    % Palettes
    [readable_colourmap,...
    first_band_template, second_band_template, third_band_template] = universally_readable_colourmap(...
        palette_code,...
        blind_conversion_method,...
        wanted_colours_per_band,...
        kovesi_modification_flag,...
        display_single_palette_flag,...
        display_multiple_palette_flag);
    % Colours representing each of the three Frequency sub-bands
    colour_1st_band = first_band_template(1,:); % [0...255]
    colour_2nd_band = second_band_template(1,:); % [0...255]
    colour_3rd_band = third_band_template(1,:); % [0...255]
    colour_1st_band = rgb_to_colour_blindness(colour_1st_band, blindness_type, blind_conversion_method); % [0...255]
    colour_2nd_band = rgb_to_colour_blindness(colour_2nd_band, blindness_type, blind_conversion_method); % [0...255]
    colour_3rd_band = rgb_to_colour_blindness(colour_3rd_band, blindness_type, blind_conversion_method); % [0...255]
    colour_1st_band = colour_1st_band/255; % [0...1]
    colour_2nd_band = colour_2nd_band/255; % [0...1]
    colour_3rd_band = colour_3rd_band/255; % [0...1]
    main_colours = [colour_1st_band;...
                    colour_2nd_band;...
                    colour_3rd_band];
    % Build colormap for colorbar
    for channel = 1:3
        low_lim = RGBFreqLimitsNorm_int_over(1,channel);
        top_lim = RGBFreqLimitsNorm_int_over(2,channel);
        RGBcolorbar_over(low_lim:top_lim,:) = RGBcolorbar_over(low_lim:top_lim,:) + ...
            repmat(main_colours(channel,:),[top_lim-low_lim+1 1]);
    end
    RGBcolorbar_over = min(RGBcolorbar_over,1);
    % Calculate images
    azFilteredRGBmod_over = zeros(size(azFilteredRGBunique_over), 'uint8');
    azFilteredRGBmod_unique_over = zeros(size(azFilteredRGBunique_over), 'uint8');
    for idx_band = 1:3
        % Modified RGB
        azFilteredRGBmod_over(:,:,idx_band) = ...
            colour_1st_band(idx_band)*azFilteredRGB_over(:,:,1) + ...
            colour_2nd_band(idx_band)*azFilteredRGB_over(:,:,2) + ...
            colour_3rd_band(idx_band)*azFilteredRGB_over(:,:,3); % [0...255]
        % Intensity of most intense band
        azFilteredRGBmod_unique_over(:,:,idx_band) = ...
            colour_1st_band(idx_band)*azFilteredRGBunique_over(:,:,1) + ...
            colour_2nd_band(idx_band)*azFilteredRGBunique_over(:,:,2) + ...
            colour_3rd_band(idx_band)*azFilteredRGBunique_over(:,:,3); % [0...255]
    end
    azFilteredRGBmod_unique_over = uint8(azFilteredRGBmod_unique_over);
    % Modified RGB
    if RGB_over_full_universal_disp_flag
        current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
        figure_objects = [figure_objects current_figure];
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod_over),grid on,
            ylabel('Depth (km)')
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod_over),grid on,axis xy
            ylabel('Height a.s.l. (km)'),
        end
        set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
        xlabel(xLabelString),
        title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
        colormap(RGBcolorbar_over),
        auxCBar=colorbar('FontSize',18,'Direction','Reverse');
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm_int_over)-1)/(numberOfColors_over-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits_over.'/100*foundPRF/recPulseDiff*10)/10)))))
        set(auxCBar,'TickLength',0)
    end
    % Display image with the intensity of the most intense band in each SAR sample
    if RGB_over_unique_universal_disp_flag
        current_figure = figure; set(gcf,'Color','w', 'Units', figure_units, 'Position', figure_position_norm),
        figure_objects = [figure_objects current_figure];
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod_unique_over),grid on,
            ylabel('Depth (km)')
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod_unique_over),grid on,axis xy
            ylabel('Height a.s.l. (km)'),
        end
        set(gca,'FontSize',18, 'Color', [0.8 0.8 0.8]),
        xlabel(xLabelString),
        title(['SAR image, RGB Doppler, ' recFlight ', ' TXWVstring '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
        colormap(RGBcolorbar_over),
        auxCBar=colorbar('FontSize',18,'Direction','Reverse');
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm_int_over)-1)/(numberOfColors_over-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits_over.'/100*foundPRF/recPulseDiff*10)/10)))))
        set(auxCBar,'TickLength',0)
    end
end
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% Mosaic of axes, enumerated as:
%   from left to right (if any)
%   from top to bottom (if any)
%
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% Mosaic type
if mosaic_type
    % Figure size
    mosaic_figure_units = 'centimeters';
    mosaic_figure_position_norm = [0    0.0565    1    0.8722]; % normalised units
    % Axes sizes
    mosaic_axes_units = 'centimeters';
    mosaic_colorbar_units = 'centimeters';
    % Assign labels
    mosaic_xlabel_text = xLabelString;
    if ~aslFlag % thickness
        mosaic_ylabel_text = 'Depth (km)';
    else
        mosaic_ylabel_text = 'Height a.s.l. (km)';
    end
    % Axes sizes
    switch mosaic_type
        case 1 % 2x2
            mosaic_axes = gobjects(4,1);
            axis_position_cm = cell(numel(mosaic_axes),1);
            % Labels
            xlabel_cell = {'', '', mosaic_xlabel_text, mosaic_xlabel_text};
            ylabel_cell = {mosaic_ylabel_text, '', mosaic_ylabel_text, ''};
            % Definitions
            top_offset = 15; % [axes_units]
            left_offset = 4; % [axes_units]
            bottom_offset = 6.5; % [axes_units]
            right_offset = 14.5; % [axes_units]
            axis_height = 7; % [axes_units]
            axis_width = 7; % [axes_units]
            colorbar_x_offset = 0.1; % [axes_units]
            colorbar_width = 0.3; % [axes_units]
            % Assignations
            axis_position_cm{1} = [ left_offset,    top_offset, axis_width, axis_height];
            axis_position_cm{2} = [right_offset,    top_offset, axis_width, axis_height];
            axis_position_cm{3} = [ left_offset, bottom_offset, axis_width, axis_height];
            axis_position_cm{4} = [right_offset, bottom_offset, axis_width, axis_height];
            colorbar_position_cm{1} = [sum(axis_position_cm{1}([1 3]))+colorbar_x_offset, axis_position_cm{1}(2), colorbar_width, axis_height];
            colorbar_position_cm{2} = [sum(axis_position_cm{2}([1 3]))+colorbar_x_offset, axis_position_cm{2}(2), colorbar_width, axis_height];
            colorbar_position_cm{3} = [sum(axis_position_cm{3}([1 3]))+colorbar_x_offset, axis_position_cm{3}(2), colorbar_width, axis_height];
            colorbar_position_cm{4} = [sum(axis_position_cm{4}([1 3]))+colorbar_x_offset, axis_position_cm{4}(2), colorbar_width, axis_height];
        case 2 % 3x1
            mosaic_axes = gobjects(3,1);
            axis_position_cm = cell(numel(mosaic_axes));
            % Figure size
            mosaic_figure_position_norm = [2    10   17.4   6.4]; % 'mosaic_figure_units'
            % Labels
            xlabel_cell = {mosaic_xlabel_text, mosaic_xlabel_text, mosaic_xlabel_text};
            ylabel_cell = {mosaic_ylabel_text, '', ''};
            % Definitions
            top_offset = 1.0; % [axes_units]
            left_offset = 1.6; % [axes_units]
            axis_height = 4.0; % [axes_units]
            axis_width = 4.0; % [axes_units]
            horizontal_offset = 1.3 + axis_width; % [axes_units]
            colorbar_width = 0.3; % [axes_units]
            colorbar_x_offset = 0.1; % [axes_units]
            % Assignations
            axis_position_cm{1} = [left_offset+0*horizontal_offset, top_offset, axis_width, axis_height];
            axis_position_cm{2} = [left_offset+1*horizontal_offset, top_offset, axis_width, axis_height];
            axis_position_cm{3} = [left_offset+2*horizontal_offset+0.7, top_offset, axis_width, axis_height];
            colorbar_position_cm{1} = [axis_position_cm{1}(1)-colorbar_x_offset-colorbar_width, axis_position_cm{1}(2), colorbar_width, axis_height];
            colorbar_position_cm{2} = [sum(axis_position_cm{2}([1 3]))+colorbar_x_offset, axis_position_cm{2}(2), colorbar_width, axis_height];
            colorbar_position_cm{3} = [axis_position_cm{3}(1)-colorbar_x_offset-colorbar_width, axis_position_cm{3}(2), colorbar_width, axis_height];
    end
    axis_font_size = 9;
    bar_font_size = 9;
    % Build mosaic
    RGBcolorbar_over = 'jet';
    current_figure = figure;
    set(gcf,'Color','w', 'Units', mosaic_figure_units, 'Position', mosaic_figure_position_norm),
    % Axis 1
    axis_index = 1;
    if length(mosaic_axes) >= axis_index
        mosaic_axes(axis_index) = axes;
        if ~aslFlag % thickness
            imagesc(alongTrackAxis, recRadialDepthGrid/1e3, 20*log10(abs(recSARimage)/recMaxAbsoluteSAR)),
        else % above sea level
            imagesc(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, 20*log10(abs(recSARimage_SL)/recMaxAbsoluteSAR)), axis xy,
        end
        grid on, caxis(-dBLim(end:-1:1)),
        set(gca, 'FontSize', axis_font_size, 'Color', [0.8 0.8 0.8], 'Units', mosaic_axes_units, 'Position', axis_position_cm{axis_index}),
        %set(gca, 'XAxisLocation', 'top'),
        set(gca, 'YAxisLocation', 'Right', 'GridColor', [1 1 1], 'GridLineStyle', '--')
        xlabel(xlabel_cell{axis_index}), ylabel(ylabel_cell{axis_index}),
        title({[recFlight ', ' TXWVstring(1:6) '-' RXstring],['Range-focussed only']})
        used_title = get(mosaic_axes(axis_index), 'Title');
        set(used_title,'Units','centimeters')
        set(used_title,'Position',get(used_title,'Position') - [0 0.1 0])
        colormap(gca, 'gray'),
        auxCBar = colorbar('FontSize', bar_font_size, 'Units', mosaic_colorbar_units, 'Position', colorbar_position_cm{axis_index});
        auxCBarLabel = auxCBar.Label;
        set(auxCBarLabel,'String','Amplitude (dB)','Rotation',-90,'VerticalAlignment','top');
        set(auxCBar,'Ticks',linspace(-max(dBLim),-min(dBLim),3),'TickLabels',linspace(-max(dBLim),-min(dBLim),3)), clear auxCBar,
    end
    % Axis 2
    axis_index = 2;
    if length(mosaic_axes) >= axis_index
        mosaic_axes(axis_index) = axes;
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGB), grid on,
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGB), grid on, axis xy
        end % above sea level
        set(gca, 'FontSize', axis_font_size, 'Color', [0.8 0.8 0.8], 'Units', mosaic_axes_units, 'Position', axis_position_cm{axis_index}),
        set(gca, 'GridColor', [1 1 1], 'YTickLabels',[], 'GridLineStyle', '--')
        xlabel(xlabel_cell{axis_index}), ylabel(ylabel_cell{axis_index}),
        title({[recFlight ', ' TXWVstring(1:6) '-' RXstring], ['RGB Doppler, aperture ' num2str(decomposition_processed_band) '\circ']})
        used_title = get(mosaic_axes(axis_index), 'Title');
        set(used_title,'Units','centimeters')
        set(used_title,'Position',get(used_title,'Position') - [0 0.1 0])
        colormap(gca, RGBcolorbar),
        auxCBar = colorbar('FontSize', bar_font_size, 'Direction', 'Reverse', 'Units', mosaic_colorbar_units, 'Position', colorbar_position_cm{axis_index});
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
    end
    % Axis 3
    axis_index = 3;
    if length(mosaic_axes) >= axis_index
        mosaic_axes(axis_index) = axes;
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod_unique), grid on,
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod_unique), grid on, axis xy
        end % above sea level
        set(gca, 'FontSize', axis_font_size, 'Color', [0.8 0.8 0.8], 'Units', mosaic_axes_units, 'Position', axis_position_cm{axis_index}),
        set(gca, 'YAxisLocation', 'Right', 'GridColor', [1 1 1], 'YTickLabels',[], 'GridLineStyle', '--')
        xlabel(xlabel_cell{axis_index}), ylabel(ylabel_cell{axis_index}),
        title({[recFlight ', ' TXWVstring(1:6) '-' RXstring], ['RGB Doppler, aperture ' num2str(decomposition_processed_band) '\circ']})
        used_title = get(mosaic_axes(axis_index), 'Title');
        set(used_title,'Units','centimeters')
        set(used_title,'Position',get(used_title,'Position') - [0 0.1 0])
        colormap(gca, modRGBcolorbar),
        auxCBar = colorbar('FontSize', bar_font_size, 'Direction', 'Reverse', 'Units', mosaic_colorbar_units, 'Position', colorbar_position_cm{axis_index});
        %auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'TickLabels','')
        %set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
    end
    % Axis 4
    axis_index = 4;
    if length(mosaic_axes) >= axis_index
        mosaic_axes(axis_index) = axes;
        if ~aslFlag % thickness
            image(alongTrackAxis, recRadialDepthGrid/1e3, azFilteredRGBmod), grid on,
        else
            image(alongTrackAxis, heightASLrecRadialDepthGrid/1e3, azFilteredRGBmod), grid on, axis xy
        end % above sea level
        set(gca, 'FontSize', axis_font_size, 'Color', [0.8 0.8 0.8], 'Units', mosaic_axes_units, 'Position', axis_position_cm{axis_index}),
        xlabel(xlabel_cell{axis_index}), ylabel(ylabel_cell{axis_index}),
        title(['RGB Doppler, ' recFlight ', ' TXWVstring(1:6) '-' RXstring ', beam center ' num2str(recAlongTrackCentralAngle) '\circ, aperture ' num2str(recAlongTrackApertureAngle) '\circ'])
        used_title = get(mosaic_axes(axis_index), 'Title');
        set(used_title,'Units','centimeters')
        set(used_title,'Position',get(used_title,'Position') - [0 0.1 0])
        colormap(gca, modRGBcolorbar),
        auxCBar = colorbar('FontSize', bar_font_size, 'Direction', 'Reverse', 'Units', mosaic_colorbar_units, 'Position', colorbar_position_cm{axis_index});
        auxCBarLabel = auxCBar.Label; set(auxCBarLabel,'String','Doppler spectrum (Hz)','Rotation',-90,'VerticalAlignment','bottom'); 
        set(auxCBar,'Ticks',(unique(RGBFreqLimitsNorm64)-1)/(numberOfColors-1),'TickLabels',strtrim(cellstr(num2str(unique(round(RGBFreqLimits.'/100*foundPRF/recPulseDiff*10)/10)))))
    end
end
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% Link axes properties
%
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
% Single axes
single_axes = findobj(figure_objects, 'Type', 'Axes');
% Single and mosaic axes
if mosaic_type
    all_axes_interest = [single_axes; mosaic_axes];
else
    all_axes_interest = single_axes;
end
% Link all axes and limits
linkaxes(all_axes_interest, 'xy');
ylim(all_axes_interest, vertical_limits_shown);
xlim(all_axes_interest, along_track_limits);
% Change how to display the along-track axes
if display_track_distance_flag
    shown_along_track_limits = xlim;
    % Indexes of first and last along-track indexes shown
    [~, displayed_track_first] = min(abs(alongTrackAxis - shown_along_track_limits(1))); % [samples]
    [~, displayed_track_last] = min(abs(alongTrackAxis - shown_along_track_limits(end))); % [samples]
    xlim(all_axes_interest, alongTrackAxis([displayed_track_first displayed_track_last]));
    % Track distance between the first and last along-track indexes shown
    displayed_track_distance = track_distance(displayed_track_last) - track_distance(displayed_track_first); % [m]
    % Distance corrected to the first shown along-track indexe
    track_distance_corrected = track_distance - track_distance(displayed_track_first); % [m]
    % Interval to be shown within the shown track distance
    track_distance_interval = floor(displayed_track_distance/number_display_intervals/precission_m)*precission_m; % [m]
    % Distances to be shown
    track_distance_ticks = (0:number_display_intervals - 1)*track_distance_interval; % [m]
    % Indexes of the distances to be shown
    index_closest_distance = zeros(1, length(track_distance_ticks)); % [samples]
    for index_closest = 1:length(index_closest_distance)
        [~, index_closest_distance(index_closest)] = min(abs(track_distance_corrected - track_distance_ticks(index_closest))); % [samples]
    end
    % Labels of x-axis
    xlabel(all_axes_interest, 'Along-track (km)');
    if mosaic_type
        for axis_index = 1:length(mosaic_axes)
            if isempty(xlabel_cell{axis_index})
                xlabel(mosaic_axes(axis_index), '');
            end
        end
    end
    % X-ticks of the along-track indexes corresponding to the distances to be shown
    xticks(all_axes_interest, alongTrackAxis(index_closest_distance)); % [kilotraces]
    % Labels of the x-ticks
    xticklabels(all_axes_interest, track_distance_ticks/1e3); % [km]
    % Special labels of the upper x-ticks, made with new axis
    upper_axes = gobjects(length(mosaic_axes),1);
    top_axis_normalized_xticks = (alongTrackAxis(index_closest_distance) - alongTrackAxis(displayed_track_first))/diff(alongTrackAxis([displayed_track_first displayed_track_last])); % [0-1]
    for upper_index = 1:length(upper_axes)
        % Create new axis
        upper_axes(upper_index) = axes('Units', mosaic_axes(upper_index).Units);
        % Properties of new axis
        set(upper_axes(upper_index),...
            'Position', [mosaic_axes(upper_index).Position(1) sum(mosaic_axes(upper_index).Position([2 4])) mosaic_axes(upper_index).Position(3) 0], ...
            'XAxisLocation', 'top', ...
            'YAxisLocation', 'left', ...
            'Color', 'none', ...
            'XTick', top_axis_normalized_xticks, ...
            'XTickLabels', round(alongTrackAxis(index_closest_distance),2), ...
            'YTick', []);
        xlabel(upper_axes(upper_index), xlabel_cell{upper_index});
        if upper_index == 2
            title(upper_axes(upper_index), [recFlight ', ' TXWVstring(1:6) '-' RXstring ...
                ', Range-focussed only, RGB Doppler aperture ' num2str(decomposition_processed_band) '\circ']);
        end
        title(mosaic_axes(upper_index), '');
    end
end
% set(gca, 'XTickLabelMode', 'auto')
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
%
%
% Text-box labels
%
%
%----------------------------------------------------------------
%----------------------------------------------------------------
%----------------------------------------------------------------
if add_box_labels_in_mosaic_flag
    % Row-wise format
    textbox_FontSize = 9;
    textbox_units = 'centimeters';
    textbox_width = 0.43;
    textbox_height = 0.41;
    textbox_h_offset = 0.3;
    textbox_v_offset = 0.1;
    horizontal_shift = 0:length(mosaic_axes)-1;
    for box_index = 1:length(horizontal_shift)
        % char(97) = 'a'
        t = annotation('textbox', 'String', ['(' char(97+box_index-1) ')'],...
            'FontSize', textbox_FontSize, 'Units', textbox_units,...
            'LineStyle', 'none', 'Background', 'w', 'FitBoxToText', 'off',...
            'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Middle');
        t.Position = [left_offset + horizontal_shift(box_index)*horizontal_offset + axis_width - (textbox_h_offset + textbox_width), ...
                      top_offset + axis_height - (textbox_v_offset + textbox_height), ...
                      textbox_width,...
                      textbox_height ];
        if box_index == length(horizontal_shift)
            t.Position = t.Position + [0.7 0 0 0];
        end
    end
end