%-------------------------------------------------------------------
% Parameters for SAR processing algorithmBack-Projection
%
%   '...'
%-------------------------------------------------------------------
seasonFlight = '1920';
% output directory of processed SAR image
out_dir = '...'; 
% output file format
netCDF_flag = 1; % 0- Matlab files ('.mat'); 1- netCDF file ('.nc')
% label for the short name of the output processed SAR image
out_label = '1920T04';
% extension for the short name of the output processed SAR image
out_extension = 'nc';
% sufix for the short name of the output processed SAR image
out_fileSufix = '';
% file with 'Timestamp in GPS', 'Lat', 'Long', 'Ellipsoidal Height(WGS84)', 'Roll', 'Pitch', 'Heading'...
attribute_file = '.../T04.ASC';
% file with the terrain clearance:
%   '.csv' file: obtained after initial processing with doppler processing without range migration correction:
%       traceNum,Lon,Lat,Time,dateNum,PriNum,Eht,terrainClear,surfElev,surfSource,bedElev,iceThick,surfPickEff,botPickLoc,sampOffset
%   '.mat' file: obtained after initial processing with radargram:
%       'Terrain_Clearance_est' and 'pulses_of_Terrain_Clearance_est'
bot_file = '.../1920T04_bot.csv';
%-------------------------------------------------------------
% environment parameters
%-------------------------------------------------------------
% relative refractive index (ice/free space)
nr = 1.78; % [-]
% speed of light
c0 = 299707760; % [m/s]
%-------------------------------------------------------------
% diffraction regions (incidence perpendicular to surface, followed by diffraction below surface)
%-------------------------------------------------------------
% time interval of 'bot_file'
tracesPRI = 0.4; % [s]
% vector of initial trace of diffraction region
initialTraceDiffractionRegion = -1; % number of trace
% vector of end trace of diffraction region
endTraceDiffractionRegion = 0; %[1e3*100]; % number of trace
%-------------------------------------------------------------
% radar parameters
%-------------------------------------------------------------
% carrier frequency
f0 = 150e6; % [Hz]
% pulse repetition frequency before range stacking
prf_before_rstacking = 15625; % [Hz]
% range stacking factor
rstacking = 25; % [pulses]
% number of transmitted sequences
nseqs = 5;
% sampling frequency (real)
Sfreal = 120e6; % [Hz]
% video frequency (real)
fvideoReal = 30e6; % [Hz]
% system delay
wd = 1.9e-6; % [s]
%-------------------------------------------------------------
% azimuth processing parameters
%-------------------------------------------------------------
% lenght (power of 2) of azimuth processing block
proc_block_size = 8*1024; % [samples]
% azimuth presumming factor
az_pres = 2; 
% center of angular aperture in doppler processing (from target point of view, i.e. positive when radar looks backwards)
center_angle = 0; % [deg]
% full aperture beamwidth in doppler processing
integr_angle = 30; % [deg]
%-------------------------------------------------------------
% Raw data files parameters
%-------------------------------------------------------------
dirIn = {'...','...','...'};
labelTransmitters = {'TxP1234','TxS9ABC'};
labelChirp = {'C13L4','J13L4','C13L1'};
labelReceivers = {'RxP1','RxP2','RxP3','RxP4','RxB5','RxB6','RxB7','RxB8','RxS9','RxSA','RxSB','RxSC'};
% indexes within 'labelReceivers'
rx = [1 ]; 
% label for wanted sets: 'TxP1234_C13L4','TxP1234_J13L4','TxS9ABC_C13L4','TxS9ABC_J13L4','TxP1234_C13L1'
labelSets = {'TxP1234_C13L4','TxP1234_J13L4','TxS9ABC_C13L4','TxS9ABC_J13L4','TxP1234_C13L1'}; 
% sets for beamforming, indexes within 'labelSets'
setsBeamf = [3 4 ]; 
% blocks of pulses (block index starts at 0) of raw data
blockIni = 151; % 17;
blockEnd = 166; % 33;
pulsesPerBlock = 2500; % pulses with prf
% blockIni = floor((azGrid(1)-1)/pulsesPerBlock); % (index starts at 0)
% blockEnd = floor((azGrid(end)-1)/pulsesPerBlock); % (index starts at 0)
% reference receiver antenna index (within 'rx') for calibration
refRX = 1;
% reference set (TX and pulse) index, within 'setsBeamf'
refSet = 1;
%-------------------------------------------------------------
% Signal parameters
%-------------------------------------------------------------
labelSignals = {'C', 'J'};
amplitudeSignals = [1, -1];
%-------------------------------------------------------------
% back-projected output grid
%-------------------------------------------------------------
% depth grid (sampling (c0/(2*nr*Brg)))
% azimuth grid of original raw data (azGrid=1 corresponds to blockIni=0), without presumming
azGrid = (154*pulsesPerBlock+1):az_pres*1:(164*pulsesPerBlock); % [samples], (index starts at 1)
%azGrid = (710*pulsesPerBlock):az_pres*1:(718*pulsesPerBlock); % [samples], (index starts at 1)
%azGrid = azGrid([31 32 33]);
%azGrid = 777500:777500; % [samples], (index starts at 1)
%azGrid = 686200:686200; % [samples], (index starts at 1)
%azGrid = 694000:694000; % [samples], (index starts at 1)
% cross-track angle (nadir = 0, portwards is positive)
crossTrackAngleGrid = 0; % [deg]
% length of output block processed azimuths
lengthAzOutBlock = 512;
% radial distance from platform center
%radialGrid = 1300:1:1700; % [m]
radialDepthGrid = -180:6.48/1:3000; % [m]
%radialDepthGrid = radialDepthGrid(893);
%---------------------------------------------------
% Antenna parameters
%---------------------------------------------------
% Antenna positions when aircraft is on ground:
% [x; y; z] coordinates relative to [Xcenter, Ycenter, 0] (aircraft on ground)
% 1st coordinate [m]: orhogonal to wing, positive to front (nose) of the aircraft
% 2nd coordinate [m]: along wing, positive to port tip of the aircraft
% 3rd coordinate [m]: height, positive towards sky
portPos = [0.01 -0.01 0.003 0.000;...
           8.37510 6.76510 5.1653 3.5400;...
           2.614 2.548 2.454 2.356].'; % [m], [x; y; z]
starPos = [0.000 -0.005 -0.003 0.005;...
           -3.5100 -5.1295 -6.7403 -8.3695;...
           2.356 2.454 2.548 2.614].'; % [m], [x; y; z]
bellyPos = [-2.915 -2.915 -2.915 -2.915;...
            1.469 0.500 -0.505 -1.469;...
            0.850 0.850 0.850 0.850].'; % [m], [x; y; z]
leverArms = [-2.323 0.012 1.435]; % [m], [x y z]
%-------------------------------------------------------------
% navigation data parameters
%-------------------------------------------------------------
firstRowIMUfile = '7052'; % row in 'IMUfile'
lastRowIMUfile = '10373'; % row in 'IMUfile'
columnUTCtime = 'A'; % column in 'IMUfile'
columnGPSCOG = 'P'; % column in 'IMUfile'
columnroll = 'E'; % column in 'IMUfile'
columnpitch = 'F'; % column in 'IMUfile'
columnheading = 'G'; % column in 'IMUfile'
columnHEllipsoid = 'D'; % column in 'IMUfile'
columnlongitude = 'C'; % column in 'IMUfile'
columnlatitude = 'B'; % column in 'IMUfile'
IMUGPSheight = 1.2; % [m], probably is in 'IMUFILE'
% reference ellipsoid
referenceEllipsoidLabel = 'WGS84';
% time format:
%   0: number of days from January 0, year 0000
%   1: seconds of week (week starts on Sunday!)
%   2: fractional part of day since current day at 00:00
timeFormatFlag = 0;
% sampling frequency of GPS
GPSfreq = 1; % [Hz]
% offset from IMU/GPS to raw data [pulses after range stacking]
GPSManualOffsetInPulses = 2545; % [pulses]
% interpolate from GPS, with GPS_time = 'GPS_delay' + UTC_time
GPS_delay = 18; % (s)
%-------------------------------------------------------------
% 'bot' (from 'bottom') data parameters, after previous fast doppler processing
%-------------------------------------------------------------
firstRowBotfile = 2; % row in 'BotFile'
lastRowBotfile = 8251; % row in 'BotFile'
numColumnDateNumber = 5; % column-th in 'BotFile'
numColumnTerrainClear = 8; % column-th in 'BotFile'
%-------------------------------------------------------------
% processing flags
%-------------------------------------------------------------
% flag for debugging: plots echoes and echoe locations
debug1Flag = 0; % 1:debug mode; 0:no debug
indAzDebug = 2800; %2;%2275; % sample within 'azGrid'
anglesDebug = [-10 0 +10]; % (deg)
% flag for debugging: stores azimuth reference function and samples
debug2Flag = 0; % 1:debug mode; 0:no debug
% subtract background noise in each receiver
filterBGNoise_flag = 0; % 1:filter; 0:no filter (use drift values)
% peaks filtering real (the filter generates broadband spectral noise and performs badly when the peak is within a high frequency signal!!!)
filterNoisyPeakReal_flag = 1; % 1:filter; 0:no filter
filterNoisyPeakRealManual_flag = 1; % 1:manual filtering; 0:automatic
% IQ conversion parameters
IQconversion_flag = 1; % 0: not convert to IQ; 1:convert to IQ;
flagDC = 1; % 0:filter DC; 1:no filter DC
IQmethode = 1; % 1:with FFT; 2:with Hilbert filter
% peaks filtering IQ (the filter generates broadband spectral noise !!!)
filterNoisyPeakIQ_flag = 0; % 1:filter; 0:no filter
% doppler filter (remove doppler frequencies higher than the +-2*v/lambda)
dopplerFilter_flag = 1; % 1:boxcar filter; 0:no filter
% Rotation correction
rotationCorrection_flag = 0; % 0:not perform correction; 1: perform
% Constant calibration correction
constantAPDCalibrationCorrection_flag = 1; % 0:not perform correction; 1: perform
% Constant calibration correction definition
RxCalDependOnTx_flag = 1; % 0: as a column vector for the sets and a row vector for the receivers; 1: as matrices (receivers x sets)
% Constant calibration correction source
constantAPDCalibrationCorrectionFile_flag = 0; % 0: manually (within the actual file); 1: from calibration file
if constantAPDCalibrationCorrectionFile_flag
constantAPDCalibrationMatFile = '...'; 
else % manual
    labelReceiversConstantCalibration = labelReceivers; % row cell
    labelSetsConstantCalibration = {'TxP1234_C13L4';'TxP1234_J13L4';'TxP1234_C13L1';'TxS9ABC_C13L4';'TxS9ABC_J13L4'}; % column cell
    if ~RxCalDependOnTx_flag % 0: as a column vector for the sets and a row vector for the receivers; 1: as matrices (receivers x sets)
        amplitudeErrorRX = 1*ones(size(labelReceiversConstantCalibration)); % [-] (Receivers x 1)
        amplitudeErrorSet = 1*ones(size(labelSetsConstantCalibration)); % [-] (1 x Sets)
        phaseErrorRX = 0*ones(size(labelReceiversConstantCalibration)); % (rad) (Receivers x 1)
        phaseErrorSet = 0*ones(size(labelSetsConstantCalibration)); % (rad) (1 x Sets)
        timeDelayErrorRX = 0*ones(size(labelReceiversConstantCalibration)); % (s) (Receivers x 1)
        timeDelayErrorSet = 0*ones(size(labelSetsConstantCalibration)); % (s) (1 x Sets)
    else % matrices for correction parameters
        load('...');
        %amplitudeErrorSetRX = repmat([amplitudeErrorSetRX(1:4,1); amplitudeErrorSetRX(5:8,4)], [1 length(labelSetsConstantCalibration)]); % [-] (Receivers x Sets)
        amplitudeErrorSetRX = repmat(10.^([-2.57;-1.92;0;-5.02;-1.40;0.24;0;0.49]/20), [1 length(labelSetsConstantCalibration)]); % [-] (Receivers x Sets)
        %phaseErrorSetRX = [[102.98;118.11;0;-67.21;0;0;0;0] 180+[102.98;118.11;0;-67.21;0;0;0;0] [102.98;118.11;0;-67.21;0;0;0;0] 180+[102.98;118.11;0;-67.21;0;0;0;0] [102.98;118.11;0;-67.21;0;0;0;0]]*pi/180; % (rad) (Receivers x Sets)
        phaseErrorSetRX = phaseErrorSetRX; % (rad) (Receivers x Sets)
        %timeDelayErrorSetRX = repmat([9.60;10.09;0;7.59;9.60;10.09;0;7.59]*1e-9, [1 length(labelSetsConstantCalibration)]); % (s) (Receivers x Sets)
        timeDelayErrorSetRX = timeDelayErrorSetRX; % (s) (Receivers x Sets)
        %amplitudeErrorSetRX = 1*ones(length(labelReceivers), length(labelSetsConstantCalibration)); % [-] (Receivers x Sets)
        %phaseErrorSetRX = 0*ones(length(labelReceivers), length(labelSetsConstantCalibration)); % (rad) (Receivers x Sets)
        %timeDelayErrorSetRX = 0*ones(length(labelReceivers), length(labelSetsConstantCalibration)); % (s) (Receivers x Sets)
    end
end
% angular back-scattering correction
angularBackScatteringCorrection_flag = 0; % 0:not perform correction; 1: perform
angularBackScatteringFile = '...';
% Beamforming
beamforming_flag = 0; % 0:not perform beamforming; 1: perform
% range compression
compression_flag = 1; % 0:not perform range compression; 1:perform
slopeTX = 1; % slope sign of transmitted pulse
nominalChirp_flag = 1; % 0:chirp on file; 1:nominal (expected) chirp
minGD_flag = 1; % 0- do not shift; 1- shift, calculating delay with instant frequencies
                % 2- shift, calculating delay with group delay
win_flag = 1;
% presumming filter
azPres_flag = 0; 
filterPres_flag = 3; % 3- boxcar in doppler; 2- triangle in time; 1- boxcar in time; 0- none
% range interpolation factor
intFactor = 1; 
% Back-projection flags
filterMigrationDC_flag = 0; % 0: no filter; 1: filter DC of BP samples
% beam-steering or null-steering
flagNS = 0;  % 0- beam-steering; 1- null-steering
flagMultipleSections = 0; % 0- 1array; 1- several subarrays
lengthLinearArray = 101;%51; % [samples]
indexCentralElement = 51;%26; % [sample]
alpha2Cancel = 0; % [deg]
cancelOrderDer = 0; % [deg] cancel until 0, 1st or 2nd order derivative of the pattern at 'alpha2Cancel'
flagMethodeNoMultipleSections = 0;
% refraction angle estimation
refractionAngleIterativeFlag = 0; % 0: small-angle aprox; 1: iterative methode
numOfIntervalsRefAngle = 1*1024; % number of intervals for searching the elevation angle (power of 2) in the iterative methode
% propagation delay
propagationDelayPolynomialInterpFlag = 0; % 0: no polynomial interpolation; 1: cubic interpolation
% fractional delay filter parameters
fractFilterLength = 20; % length of causal FIR filter
fractDelays = -0.5:0.01:0.5; % [samples], vector of fractional (decimal) delays, within the interval [-0.5, 0.5]
fractDelayPolyDegree = 0; % >0: degree of polynomial interpolation; 0: no interpolation
% cross-track elevation (above ellipsoid)
crossTrackElevationFlag = 0; %1-calculate; 0-not calcula
%---------------------------------------------
% Parameters for netCDF file
%---------------------------------------------
% --- References:
%     netCDF User's Guide (NUG): main information about netCDF
%           https://www.unidata.ucar.edu/software/netcdf/documentation/4.7.4-pre/index.html
%           https://www.unidata.ucar.edu/software/netcdf/documentation/4.7.4-pre/attribute_conventions.html
% --- Data types:
% 'NC_BYTE'     int8
% 'NC_UBYTE'    uint8
% 'NC_CHAR'     char
% 'NC_SHORT'	int16
% 'NC_USHORT'	uint16
% 'NC_INT'      int32
% 'NC_UINT'     uint32
% 'NC_INT64'    int64
% 'NC_UINT64'   uint64
% 'NC_FLOAT'	single
% 'NC_DOUBLE'	double
if netCDF_flag
    %---------------------------------------------
    % DIMENSIONS in netCDF File: cell of structs, with:
    %   - 'name' field: indicator for the name of the DIMENSION;
    %   - 'length' field: indicator for the length of the DIMENSION. Interger number or the standard 'NC_UNLIMITED'
    % The unlimited dimension must be the last dimension when defining the variables, because MATLAB
    %   uses FORTRAN-style ordering for netCDF (fastest-varying dimension comes first and the slowest comes last)
    %---------------------------------------------
    indexDimension = 0; % to grow, starting at 1 before the first DIMENSION
    cellDimNetCDF = {};
    % --- 'Limits' dimension, limited
    indexDimension = indexDimension + 1;
    cellDimNetCDF{indexDimension}.name = 'limits';
    cellDimNetCDF{indexDimension}.length = 2; % lower and greater limits
    indexLimitsDimension = indexDimension - 1;
    % --- 'Depth grid' dimension, limited
    indexDimension = indexDimension + 1;
    cellDimNetCDF{indexDimension}.name = 'radialDepthGrid';
    cellDimNetCDF{indexDimension}.length = length(radialDepthGrid);
    indexRadialDepthGridDimension = indexDimension - 1;
    % --- 'Profile' dimension (along-track), unlimited
    indexDimension = indexDimension + 1;
    cellDimNetCDF{indexDimension}.name = 'profile';
    cellDimNetCDF{indexDimension}.length = 'NC_UNLIMITED';
    indexProfileDimension = indexDimension - 1;
    %---------------------------------------------
    % Global ATTRIBUTES in netCDF File: cell of structs, with:
    %   - 'name' field: indicator for the name of the ATTRIBUTE;
    %   - 'value' field: indicator for the ATTRIBUTE value
    %---------------------------------------------
    indexGA = 0; % to grow, starting at 1 before the first Global ATTRIBUTES
    cellGlobalAttNetCDF = {};
    keySep = ', ';
    % --- Title (defined by NUG (NetCDF User's Guide))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'title';
    cellGlobalAttNetCDF{indexGA}.value = ['Ice-sounding depth profiles from airborne synthetic-aperture-radar (SAR) PASIN2 instrument, ' ...
                                          'over the Rutford Ice Stream in the Antarctic summer 2019/2020'];
    % ******** Summary (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'summary';
    cellGlobalAttNetCDF{indexGA}.value = [... % summary of PASIN2
              'PASIN2 (Polarimetric Airborne Scientific Instrument, mark 2) is a 150 MHz coherent pulsed radar with the purpose ' ...
              'of deep ice sounding for bedrock, subglacial channels and ice-water interface detection in Antarctica. PASIN2 aims ' ...
              'to map in 3D the bottom of glaciers, either bedrock or basal water. It is designed and operated by the British ' ...
              'Antarctic Survey (BAS). With multiple antennas for transmission and reception, PASIN2 enables polarimetric ' ...
              '3D estimation of the bottom profile with a single pass, reducing the gridding density of the survey paths. ' ...
              'Single 2D SAR images present ambiguities for distinguishing between scatterers from port and starboard directions. ' ...
              'This occurs because the antennas are pointing vertically downwards, to take advantage of the greater radar cross-section ' ...
              'under perpendicular incidence. The solution adopted for PASIN2 was an across-track physical array. After processing several ' ...
              '2D SAR images (range and along-track dimensions) with transmitter-receiver pairs, the directional ambiguities are ' ...
              'resolved, obtaining the across-track Direction of Arrival (DoA, elevation angle) estimation. Finally, from the 3D ' ...
              'geometry of range, along-track and across-track angles, the real depths and across-track distances are estimated, ' ...
              'regarding the case of wrongly assumed vertical DoA of a single SAR image.' newline...
          ... % summary of PASIN2 processing 
              'The data processing is performed off-line, consisting in: 1) channel calibration; '...
              '2) 2D synthetic aperture radar (SAR) imaging, based on pulse-compression for range dimension, and back-projection for along-track dimension; ' ...
              '3) Direction of Arrival estimation (DoA) of the remaining across-track angle, based on the non-linear MUSIC algorithm, ' ...
                  'by using several SAR images from different transmitters or receivers; ' ...
          'and 4) 3D-mapping from the three dimensions range, along-track and across-track angles, by correcting the real depths and ' ...
                  'across-track distances, regarding the case of wrongly-assumed vertical DoA of a single SAR image' ...
              'Calibration flights, during the Antarctic Summer campaigns in 16/17 and 19/20 seasons, assessed and validated the ' ...
              'instrument and processing performances.' newline...
          ... % summary of current netCDF 
              'This netCDF file corresponds to a SAR image of the data collected during the Antarctic Summer of 2019/2020, ' ...
              'in the Rutford Ice Stream, near the grounding line. The SAR image is 2D, with ' ...
              'vertical ice-sounding profiles along the aircraft trajectory, one profile per along-track location. Each profile ' ...
              'belongs to the netCDF dimension and variable ''profile'', and it contains the vertical grid defined by ' ...
              'the netCDF dimension and variable ''radialDepthGrid''. For each of these profiles, besides the real and imaginary ' ...
              'components of the complex SAR-processed profiles, other auxiliary data are included, such as time, pulse number, ' ...
              'positioning (latitude, longitude, height, speed, terrain clearance) and attitude data (roll, pitch, yaw) of the aircraft, ' ...
                          'the electromagnetic propagation model and others.'];
    % --- history (defined by NUG (NetCDF User's Guide))
    % Provides an audit trail for modifications to the original data.
    % Well-behaved generic netCDF filters will automatically append their name and the parameters with which they were invoked to the global history attribute of an input netCDF file.
    % We recommend that each line begin with a timestamp indicating the date and time of day that the program was executed.
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'history';
    cellGlobalAttNetCDF{indexGA}.value = 'To be defined!!!';
    % --- conventions (defined by NUG (NetCDF User's Guide))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'conventions';
    cellGlobalAttNetCDF{indexGA}.value = 'CF-1.8 ACDD-1.3 NCEI-2.0';
    % --- institution (defined by CF (Climate and Forecast) Metadata Conventions)
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'institution';
    cellGlobalAttNetCDF{indexGA}.value = 'British Antarctic Survey';
    % --- program (defined by ACDD (Attribute Convention for Dataset Discovery))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'program';
    cellGlobalAttNetCDF{indexGA}.value = 'Scientific Program';
    % --- project (defined by ACDD (Attribute Convention for Dataset Discovery))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'project';
    cellGlobalAttNetCDF{indexGA}.value = 'BEAMISH, Bed Access, Monitoring and Ice Sheet History';
    % --- source (defined by CF (Climate and Forecast) Metadata Conventions)
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'source';
    cellGlobalAttNetCDF{indexGA}.value = 'Airborne ice-penetrating synthetic aperture radar, image by Backprojection';
    % --- instrument/sensor (defined by ACDD (Attribute Convention for Dataset Discovery))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'instrument';
    cellGlobalAttNetCDF{indexGA}.value = ['Earth Remote Sensing Instruments' keySep 'Active Remote Sensing' keySep 'Profilers/Sounders' keySep 'Radar Sounders' keySep 'PASIN2' keySep newline ...
                                          'Polarimetric radar Airborne Science Instrument (PASIN2)'];
    % --- platform (defined by ACDD (Attribute Convention for Dataset Discovery))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'platform';
    cellGlobalAttNetCDF{indexGA}.value = 'De Havilland DHC-6 Twin Otter aircraft';
    % --- season
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'season';
    cellGlobalAttNetCDF{indexGA}.value = seasonFlight;
    % --- campaign (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'campaign';
    cellGlobalAttNetCDF{indexGA}.value = seasonFlight;
    % --- flight
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'flight';
    cellGlobalAttNetCDF{indexGA}.value = out_label;
    % --- reference ellipsoid
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'reference_ellipsoid';
    cellGlobalAttNetCDF{indexGA}.value = referenceEllipsoidLabel;
    % --- Radar parameters (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'radar_parameters';
    cellGlobalAttNetCDF{indexGA}.value = ['Bistatic PASIN2 radar system operating with a centre frequency of ' num2str(f0/1e6) ' MHz and using five interleaved pulses: ' ...
                                          'Chirp: 4 microseconds, 13 MHz bandwidth linear chirp, with 0deg and 180deg modulations, from port and starboard arrays; ' ...
                                          'Chirp: 1 microseconds, 13 MHz bandwidth linear chirp, with 0deg modulation, from port array; ' ...
                                          'System Pulse Repetition Frequency: ' num2str(prf_before_rstacking) ' Hz (pulse repetition interval: ' num2str(1e6/prf_before_rstacking) ' microseconds)'];
    % --- Antennas (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'antennas';
    cellGlobalAttNetCDF{indexGA}.value = ['8 folded dipole elements: 4 transmitters/receivers (port side), 4 transmitters/receivers (starboard side); ' ...
                                          '4 printed dipole elements: 4 receivers (belly); ' ...
                                          'Antenna gain: 11 dBi wing arrays, 8 dBi Belly array; Transmit power: 0.5 kW into each 4 antennas; Maximum transmit duty cycle: 10% at full power (4 x 0.5 kW)'...
                                          'Antenna labels, from port tip to starboard tip: P1, P2, P3, P4, B5, B6, B7, B8, S9, SA, SB, SC'];
    % --- Digitiser (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'digitiser';
    cellGlobalAttNetCDF{indexGA}.value = ['Receiver sampling frequency: ' num2str(Sfreal/1e6) ' MHz (resulting in sampling interval of ' num2str(1e9/Sfreal) ' nanoseconds); '...
                                          'Receiver coherent stacking: ' num2str(rstacking) ' pulses; ' ...
                                          'Receiver digital filtering: 13 MHz; Effective PRF per each of the five pulses: ' num2str(prf_before_rstacking/(nseqs*rstacking)) ' Hz (post-hardware stacking); '...
                                          'Sustained data rate per pulse and receiver: 150 MB/km'];
    % --- Processing (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'processing';
    cellGlobalAttNetCDF{indexGA}.value = ['Processing after conversion from real domain to in-phase and quadrature (IQ) domain; chirp compression; ' ...
                                          'transmitter/receiver correction of amplitude, phase and delay from calibration values; waveform summation; '...
                                          'synthetic aperture radar (SAR) processing applied to the scatterer grid locations. The grid is defined for depth and along-track locations'];
    % --- Resolutions (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'resolution';
    cellGlobalAttNetCDF{indexGA}.value = ['Vertical resolution in ice (velocity in ice of ' num2str(c0/nr/1e6) ' meters/microseconds): ' num2str(c0/(2*13e6*nr)) ' meters; '...
                                          'along-track resolution: ' num2str((c0/f0)/(4*sind(integr_angle/2))) ' meters; '...
                                          'vertical sampling: ' num2str(mean(diff(radialDepthGrid))) ' meters; '...
                                          'along-track sampling: ' num2str(mean(diff(azGrid))/(prf_before_rstacking/(nseqs*rstacking))*1e3) ' milisseconds'];
    % --- GPS sensor (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'GPS_IMU_LIDAR';
    %cellGlobalAttNetCDF{indexGA}.value = ['Dual-frequency GPS (Leica SR530 and Ashtec Z-12) Absolute GPS positional accuracy: ~0.1 m (relative accuracy is one order of magnitude better). ' ...
    %                                      'Banking angle was limited to 10deg during aircraft turns to avoid phase issues between GPS receiver and transmitter'];
    cellGlobalAttNetCDF{indexGA}.value = ['Global Positioning System (GPS) and Inertial Measurement Unit (IMU) devices are integrated, to obtain geographic ' ...
                                           'coordinates (latitude, longitude and elevation) and rotation angles (roll, pitch and yaw), sampled at ' num2str(GPSfreq) ' Hz. For ' ...
                                           'terrain clearance (distance from radar to surface) measurement, it is used a lidar (Riegl LMS-Q240i) in near-infrared, ' ...
                                           'with a maximum range of 650 m and 2 mm accuracy'];
    % --- keywords (defined by ACDD (Attribute Convention for Dataset Discovery))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'keywords';
    cellGlobalAttNetCDF{indexGA}.value = ['Earth Science' keySep 'Cryosphere' keySep 'Glaciers/Ice Sheets' keySep 'Glacier Topography/Ice Sheet Topography' keySep newline ... % GCMD location
                                          'Continent' keySep 'Antarctica' keySep 'Rutford Ice Stream' keySep newline ... % GCMD science
                                          'Antarctic' keySep 'Aerogeophysics' keySep 'Ice thickness' keySep 'Radar' keySep 'Surface elevation'];
    % --- Location (-----) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'location';
    cellGlobalAttNetCDF{indexGA}.value = 'Rutford Ice Stream';
    % --- Date of data take start (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'time_coverage_start';
    cellGlobalAttNetCDF{indexGA}.value = 'YYYY/MM/DD';
    % --- Date of data take end (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'time_coverage_end';
    cellGlobalAttNetCDF{indexGA}.value = 'YYYY/MM/DD';
    % --- featureType (defined by CF (Climate and Forecast) Metadata Conventions)
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'featureType';
    cellGlobalAttNetCDF{indexGA}.value = 'profile';
    % --- processing_level (Attribute Convention for Dataset Discovery)
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'processing_level';
    cellGlobalAttNetCDF{indexGA}.value = 'L1b';
    % --- keywords (defined by ACDD (Attribute Convention for Dataset Discovery))
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'keywords_vocabulary';
    cellGlobalAttNetCDF{indexGA}.value = 'GCMD Science Keywords Version 9.1.5';
    % --- List of standard names (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'standard_name_vocabulary';
    cellGlobalAttNetCDF{indexGA}.value = 'CF Standard Name Table v77';
    % --- comment (defined by CF (Climate and Forecast) Metadata Conventions)
    indexGA = indexGA + 1;
    commentString = ['Airborne SAR 2D-image with the incorrectly vertical-assumed depth profile (nadir direction, range) relative to the air/ice interface, along the aircraft trajectory (Doppler). ' newline ...
                     'Supporting data with the positioning and attitude of the platform are also included. ' newline ...
                     'The SAR 2D-image (complex domain) has been processed with Backprojection algorithm. ' newline ...
                     'Coregistred SAR 2D-images from different transmitters/receivers can be combined to estimate the direction of arrival (DoA) of each 2D-image pixel, and hence its true depth. ' newline ...
                     'The main processing script is included (full name and content) in the ''GLOBAL'' ATTRIBUTE ''processing_main''. ' newline ...
                     'The parameter files are included (full name and content) in the ''GLOBAL'' ATTRIBUTE ''parameter_file_ii''. '];
    cellGlobalAttNetCDF{indexGA}.name = 'comment';
    cellGlobalAttNetCDF{indexGA}.value = commentString;
    % --- references (defined by CF (Climate and Forecast) Metadata Conventions)
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'references';
    cellGlobalAttNetCDF{indexGA}.value = ['http://cfconventions.org/Data/cf-conventions/cf-conventions-1.8/cf-conventions.html' newline ...         % 'variables', 'attributes', etc.
                                    'http://cfconventions.org/Data/cf-standard-names/73/build/cf-standard-name-table.html' newline ...              % 'standard_name' conventions.
                                    'https://www.unidata.ucar.edu/software/netcdf/documentation/4.7.4-pre/index.html' newline ...                   % main information about netCDF
                                    'https://www.unidata.ucar.edu/software/netcdf/documentation/4.7.4-pre/attribute_conventions.html' newline ... 
                                    'https://wiki.esipfed.org/Category:Attribute_Conventions_Dataset_Discovery' newline ...                         % 'variables', 'attributes', etc.
                                    'https://wiki.esipfed.org/Attribute_Convention_for_Data_Discovery_1-3' newline...
                                    'https://www.nodc.noaa.gov/data/formats/netcdf/v2.0/' newline ...                                               % 'variables', 'attributes', etc.
                                    'https://earthdata.nasa.gov/earth-observation-data/find-data/gcmd/gcmd-keywords' newline ...                    % 'keywords'
                                    'https://gcmd.earthdata.nasa.gov/static/kms/'];
    % --- Acknowledgement (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'acknowledgement';
    cellGlobalAttNetCDF{indexGA}.value = ['This work was funded by the British Antarctic Survey in support of the Natural Environment Research Council (NERC). ' ...
                                          'The campaign was carried out by the British Antarctic Survey'];
    % --- License (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'license';
    cellGlobalAttNetCDF{indexGA}.value = 'http://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/';
    % --- Name of publisher (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'publisher_name';
    cellGlobalAttNetCDF{indexGA}.value = 'UK Polar Data Centre';
    % --- Type of publisher (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'publisher_type';
    cellGlobalAttNetCDF{indexGA}.value = 'institution';
    % --- Email of publisher (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'publisher_email';
    cellGlobalAttNetCDF{indexGA}.value = 'polardatacentre@bas.ac.uk';
    % --- Link of publisher (defined by ACDD (Attribute Convention for Dataset Discovery)) ********
    indexGA = indexGA + 1;
    cellGlobalAttNetCDF{indexGA}.name = 'publisher_link';
    cellGlobalAttNetCDF{indexGA}.value = 'https://www.bas.ac.uk/data/uk-pdc/';
    %---------------------------------------------
    % VARIABLES and ATTRIBUTES in netCDF File: cell of structs, with:
    %   - 'name' field: indicator for the name of the VARIABLE;
    %   - 'xtype' field: indicator for the type of the VARIABLE;
    %   - 'dimid' field: indicator for the dimension id of the VARIABLE, according to the index of 'cellDimNetCDF' - 1.
    %       In the standard, 'dimid' is indexed starting with zero (C-based), and here too, but the index of 'cellDimNetCDF' starts with 1 (Matlab-based)
    %       The unlimited dimension must be the last dimension, because MATLAB uses FORTRAN-style
    %           ordering for netCDF (fastest-varying dimension comes first and the slowest comes last)
    %   - 'FillValue' field: indicator for the standard ATTRIBUTE '_FillValue';
    %   - 'DefVarChunking' field: indicator for defining chunking with function 'netcdf.defVarChunking'
    %   - 'DefVarDeflate' field: indicator for defining the compression with function 'netcdf.defVarDeflate'
    %   - other fields: standard ATTRIBUTE names
    %---------------------------------------------
    indexVar = 0; % to grow, starting at 1 before the first VARIABLE
    cellVariableNetCDF = {};
    %---------------------
    % Fixed-size variables (non-record variables)
    %---------------------
    % --- Depth grid: nadir bottom depth below (positive) and above (negative) surface ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'radialDepthGrid';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexRadialDepthGridDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Depth grid in SAR image: nadir (vertical) direction below (positive) and above (negative) surface';
    cellVariableNetCDF{indexVar}.units = 'meters';
    cellVariableNetCDF{indexVar}.positive = 'down';
    cellVariableNetCDF{indexVar}.axis = 'Z';
    % --- Centre of angular aperture in Doppler processing (from target point of view, i.e. positive when radar looks backwards)
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'alongTrackCentralAngle';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = []; % scalar variable
    cellVariableNetCDF{indexVar}.long_name = 'Centre of angular beamwidth in Doppler processing (from target point of view, i.e. positive when radar looks backwards)';
    cellVariableNetCDF{indexVar}.units = 'degrees';
    cellVariableNetCDF{indexVar}.valid_min = -90;
    % --- Full aperture beamwidth in Doppler processing
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'alongTrackApertureAngle';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = []; % scalar variable
    cellVariableNetCDF{indexVar}.long_name = 'Full aperture beamwidth in Doppler processing';
    cellVariableNetCDF{indexVar}.units = 'degrees';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    % --- Carrier (central) frequency
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'centralFrequency';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = []; % scalar variable
    cellVariableNetCDF{indexVar}.long_name = 'Carrier (central) frequency of the system';
    cellVariableNetCDF{indexVar}.units = 'Hz';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    % --- Minimum and maximum of positive values of real SAR image
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'intervalAbsoluteSAR';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = indexLimitsDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Interval of values (minimum and maximum) of absolute value of SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    % --- Minimum and maximum of positive values of real SAR image
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'intervalPositiveRealSAR';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = indexLimitsDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Interval of values (minimum and maximum) of positive values of real SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    % --- Minimum and maximum of negative values of real SAR image
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'intervalNegativeRealSAR';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = indexLimitsDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Interval of values (minimum and maximum) of negative values of real SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.valid_max = 0;
    % --- Minimum and maximum of positive values of imaginery SAR image
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'intervalPositiveImagSAR';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = indexLimitsDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Interval of values (minimum and maximum) of positive values of imaginery SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    % --- Minimum and maximum of negative values of imaginery SAR image
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'intervalNegativeImagSAR';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = indexLimitsDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Interval of values (minimum and maximum) of negative values of imaginery SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.valid_max = 0;
    %---------------------
    % Unlimited-size variables (record variables)
    %---------------------
    % --- Profile --- (start at 0)
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'profile';
    cellVariableNetCDF{indexVar}.xtype = 'NC_UINT';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.cf_role = 'profile_id';
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Transmitted pulse --- (start at 1)
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'pulse';
    cellVariableNetCDF{indexVar}.xtype = 'NC_UINT';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Transmitted pulse within the effective PRF sequence (for the current waveform, after pulse stacking)';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    cellVariableNetCDF{indexVar}.FillValue = 2^32-1; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Pulse Time ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'pulseTime';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Time at transmitted pulse within the effective PRF sequence (for the current waveform, after pulse stacking)';
    cellVariableNetCDF{indexVar}.units = 'number of days from January 0, year 0000';
    cellVariableNetCDF{indexVar}.FillValue = -1; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Latitude ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'latitudeAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'latitude';
    cellVariableNetCDF{indexVar}.long_name = 'Latitude of aircraft at transmitted pulse';
    cellVariableNetCDF{indexVar}.units = 'degrees_north';
    cellVariableNetCDF{indexVar}.valid_min = -90;
    cellVariableNetCDF{indexVar}.valid_max = +90;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Longitude ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'longitudeAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'longitude';
    cellVariableNetCDF{indexVar}.long_name = 'Longitude of aircraft at transmitted pulse';
    cellVariableNetCDF{indexVar}.units = 'degrees_east';
    cellVariableNetCDF{indexVar}.valid_min = -360;
    cellVariableNetCDF{indexVar}.valid_max = +360;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Height of aircraft ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'heightAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'height_above_mean_sea_level';
    cellVariableNetCDF{indexVar}.long_name = 'Height of aircraft above mean sea level (asl)';
    cellVariableNetCDF{indexVar}.units = 'meters';
    cellVariableNetCDF{indexVar}.positive = 'up';
    cellVariableNetCDF{indexVar}.FillValue = -99999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Terrain clearance of aircraft ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'terrainClearanceAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Terrain clearance, distance from platform to air interface with ice, sea or ground';
    cellVariableNetCDF{indexVar}.units = 'meters';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Aircraft speed ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'speedAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'platform_speed_wrt_air';
    cellVariableNetCDF{indexVar}.long_name = 'Aircraft speed';
    cellVariableNetCDF{indexVar}.units = 'meters/second';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Aircraft roll ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'rollAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'platform_roll_starboard_down';
    cellVariableNetCDF{indexVar}.long_name = 'Aircraft roll, positive according to right-hand rule pointing towards aircraft direction (along-track)';
    cellVariableNetCDF{indexVar}.units = 'degrees';
    cellVariableNetCDF{indexVar}.valid_min = -360;
    cellVariableNetCDF{indexVar}.valid_max = +360;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Aircraft pitch ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'pitchAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'platform_pitch_fore_up';
    cellVariableNetCDF{indexVar}.long_name = 'Aircraft pitch, positive according to left-hand rule pointing towards port direction';
    cellVariableNetCDF{indexVar}.units = 'degrees';
    cellVariableNetCDF{indexVar}.valid_min = -360;
    cellVariableNetCDF{indexVar}.valid_max = +360;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Aircraft yaw ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'yawAircraft';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.standard_name = 'platform_yaw_fore_starboard';
    cellVariableNetCDF{indexVar}.long_name = 'Aircraft yaw, positive according to left-hand rule pointing upwards, perpendicular to flat sea level';
    cellVariableNetCDF{indexVar}.units = 'degrees';
    cellVariableNetCDF{indexVar}.valid_min = -360;
    cellVariableNetCDF{indexVar}.valid_max = +360;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- Ice thickness ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'iceThickness';
    cellVariableNetCDF{indexVar}.xtype = 'NC_DOUBLE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Ice thickness initially estimated, from ''bot'' file';
    cellVariableNetCDF{indexVar}.units = 'meters';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    cellVariableNetCDF{indexVar}.FillValue = -999; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- EM wave propagation model ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'EMpropagationModel';
    cellVariableNetCDF{indexVar}.xtype = 'NC_UBYTE';
    cellVariableNetCDF{indexVar}.dimid = indexProfileDimension; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'EM wave propagation model within ice: 0, by pure-refraction; 1, by refraction-diffraction';
    cellVariableNetCDF{indexVar}.valid_min = 0;
    cellVariableNetCDF{indexVar}.FillValue = 2^8-1; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = lengthAzOutBlock; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    % --- SAR image (real of complex) ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'SARimageReal';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = [indexRadialDepthGridDimension indexProfileDimension]; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Real part of complex processed SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.FillValue = 0; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = [cellDimNetCDF{indexRadialDepthGridDimension+1}.length lengthAzOutBlock]; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    cellVariableNetCDF{indexVar}.DefVarDeflate = {true 9}; % compression with {shuffleFlag, deflateLevel(0:noComp -> 9:maxComp)}, with function 'netcdf.defVarDeflate'
    % --- SAR image (imaginary of complex) ---
    indexVar = indexVar + 1;
    cellVariableNetCDF{indexVar}.name = 'SARimageImag';
    cellVariableNetCDF{indexVar}.xtype = 'NC_FLOAT';
    cellVariableNetCDF{indexVar}.dimid = [indexRadialDepthGridDimension indexProfileDimension]; % index of 'cellDimNetCDF' - 1
    cellVariableNetCDF{indexVar}.long_name = 'Imaginary part of complex processed SAR image';
    cellVariableNetCDF{indexVar}.units = 'amplitude';
    cellVariableNetCDF{indexVar}.FillValue = -1; % '_FillValue' is the standard name
    cellVariableNetCDF{indexVar}.DefVarChunking = [cellDimNetCDF{indexRadialDepthGridDimension+1}.length lengthAzOutBlock]; % chunk size according to '.dimid', with function 'netcdf.defVarChunking'
    cellVariableNetCDF{indexVar}.DefVarDeflate = {true 9}; % compression with {shuffleFlag, deflateLevel(0:noComp -> 9:maxComp)}, with function 'netcdf.defVarDeflate'
end