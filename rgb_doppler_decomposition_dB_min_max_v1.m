function azFilteredRGB = rgb_doppler_decomposition_dB_min_max_v1(SAR, ratioAz, percentageShift, xGrid, yGrid, dBLim, equalization_flag)
%
% Obtains a coloured RGB-SAR from a 'SAR' image in the along-track grid
%   'xGrid' and range grid 'yGrid', by splitting the full doppler bandwidth
%   in 3 sub-bands, asigning colours red, green and blue to the lower,
%   medium, and upper doppler frequency bands. Each band has a width
%   defined by 'ratioAz' over the full bandwidth, and the central frequency 
%   for each band is defined by 'percentageShift', a percentage over the
%   full bandwidth. Each band has a colour defined by the RGB colour model,
%   with pure red, green and blue colours, using 8 bits per colour band,
%   with integer numbers between 0 and 255 levels:
%               lower sub-band:  0 <= R <= 255
%               medium sub-band: 0 <= G <= 255
%               upper sub-band:  0 <= B <= 255
% The limits to set the darkest (0) and brightest (255) levels are given
%   by the parameter 'dBLim'>=0, in dB units below the maximum of the
%   SAR image (when the input 'equalization_flag' is false), or each
%   sub-band (when the input 'equalization_flag' is true).
%
% Steps:
%   1. The image is normalized to its maximum
%   2. For each of the 3 sub-bands:
%       a. Filtering: filter the image for the sub-band
%       b. Limit the filtered image to amplitudes within the upper limit
%          -dBLim(1) and lower limit -dBLim(2), so that
%                   -dBLim(2) <= filtered_image <= -dBLim(1)
%       c. Scaling: convert to integer numbers between 0 and 255
%       d. Assign to the sub-band
%   3. Calculate a new image with the intensity (0-255) of the most intense
%       sub-band. This image is displayed, but it is not part of the
%       output.
%   4. Calculate a new image with the maximum intensity (255) of the most 
%       intense sub-band. This image is displayed, but it is not part of
%       the output.
%
%   Input: SAR (float): [R x A] 2D-matrix with the processed SAR image in
%               complex domain, and format (range, azimuth).
%          ratioAz (float): [1 x 3] vector or number with the ratio (>1)
%               of SAR full processing band over each of the 
%               doppler sub-bandwidths.
%               Example: ratioAz = [6.25 6.25 6.25]; % ratioAz > 1
%          percentageShift (float): [1 x 3] vector or number with the 
%               shift of each central frequency, in percentage relative to
%               the full processing band. The sign of each percentage
%               depends on the sign of the central frequency of the filter.
%               Example: percentageShift = [-16, 0, 16]; % -50 < percentageShift < 50
%          xGrid (float): [1 x A] vector with the regular grid to be
%               displayed as x-coordinates.
%               Example: xGrid = 0:size(SAR,1)-1; % [samples]
%          yGrid (float): [1 x R] vector with the regular grid to be
%               displayed as y-coordinates.
%               Example: yGrid = 0:size(SAR,2)-1; % [samples]
%          dBLim (float): [1 x 2] vector with the upper and lower
%               limits in dB relative to the SAR image or sub-band maximum.
%               Example: dBLim = [0 100]; % (dB) limit below the maximum: upper=max-0; lower=max-100
%               Example: dBLim = [10 70]; % (dB) limit below the maximum: upper=max-10; lower=max-70
%          equalization_flag(logical): flag to equalize the sub-bands (true)
%               or not (false). The upper and lower limits defined by
%               'dBLim' will be applied relative to the maximum of the
%               image when the flag is not active (false), or relative to
%               the maximum of each individual sub-band when active (true).
%               Example: equalization_flag = true;
%   Output: azFilteredRGB (uint8): [R x A x 3] 3D-matrix RGB matrix with
%               the format (range, azimuth, 3=RGB).
%
%   Calling: 'azFilteredRGB = rgb_doppler_decomposition_dB_min_max_v1(SAR, ratioAz, percentageShift, xGrid, yGrid, dBLim, equalization_flag)'
%
%   Dependencies:
%           Based on: -
%           Called by: scriptBackProjectionImage3_MultipleBlocks.m
%           Calls: -
%
%   Date: 11.10.2016
%   Author: Alvaro Arenas-Pingarron, British Antarctic Survey
%
%------------------------------------------------
% Examples for 'xGrid', 'yGrid', 'dBLim' and 'equalization_flag'
%------------------------------------------------
% Examples
% xGrid = 0:size(SAR,1)-1;
% yGrid= 0:size(SAR,2)-1;
% dBLim = [0 100]; % (dB) limit below the maximum: upper=max-0; lower=max-100
% equalization_flag = true;
%------------------------------------------------
% Examples for 'ratioAz' and 'percentageShift'
%------------------------------------------------
% The sampling frequency 'Fs' of the domain of interest is always greater
%   than the band of interest 'Bw' of the data, so that Fs > Bw. If the
%   Doppler band of interest is
%           Bw = 30 Hz,
%   the band to cover with the 3 sub-bands of decomposition is the interval
%           [-Bw/2, +Bw/2] = [-15, +15]
%
% NON-OVERLAPPING BANDS:
%
%   .......\______Red______/\_____Green_____/\______Blue_____/..........
% -PRF/2  -Bw/2                     0                      +Bw/2   +PRF/2
%
% For non-overlapping bands seamlessly covering the interval [-Bw/2, +Bw/2],
%   the three sub-bandwidths would be
%           SBw1 = SBw2 = SBw3 = Bw/3 = 10 Hz,
%   and the central frequencies
%           [Fc1, Fc2, Fc3] = [-Bw/3, 0, +Bw/3] = [-10, 0, +10] Hz
% For a sampling frequency Fs = 62.5 Hz, the total available band would
%   be [-Fs/2, +Fs/2]. The input parameters 'ratioAz' and 'percentageShift'
%   would be:
%           ratioAz = Fs./[SBw1, SBw2, SBw3] = [6.25 6.25 6.25]
%           percentageShift = 100*[Fc1, Fc2, Fc3]/Fs = [-16, 0, 16]
%   so that ratioAz > 1 and -100 < percentageShift < 100
%
% For a sampling frequency Fs = 125 Hz, the input parameters would be:
%           ratioAz = Fs./[SBw1, SBw2, SBw3] = [12 12 12]
%           percentageShift = 100*[Fc1, Fc2, Fc3]/Fs = [-8, 0, 8]
%
% OVERLAPPING BANDS:
%
%   .......\_________Red_______/         \_______Blue________/..........
%   .......          \____________Green____________/         ........... 
%   .......          \__Yellow_/         \__Cyan___/          .......... 
% -PRF/2  -Bw/2                     0                      +Bw/2   +PRF/2
%
% For overlapping sub-bands seamlessly covering the interval [-Bw/2, +Bw/2],
%   in a way that the overlapped intervals have the same width as the
%   intervals without overlapping, the three sub-bandwidths would be
%           SBw1 = SBw3 = 0.4Bw = 12 Hz,
%           SBw2 = 0.6Bw = 18 Hz,
%   and the central frequencies
%           [Fc1, Fc2, Fc3] = [-0.3Bw, 0, +0.3Bw] = [-9, 0, +9] Hz
% For a sampling frequency Fs = 62.5 Hz, the total available band would
%   be [-Fs/2, +Fs/2]. The input parameters 'ratioAz' and 'percentageShift'
%   would be:
%           ratioAz = Fs./[SBw1, SBw2, SBw3] = [5.2083  3.4722  5.2083]
%           percentageShift = 100*[Fc1, Fc2, Fc3]/Fs = [-14.4  0  14.4]
%   so that ratioAz > 1 and -100 < percentageShift < 100
%
% For a sampling frequency Fs = 125 Hz, the input parameters would be:
%           ratioAz = Fs./[SBw1, SBw2, SBw3] = [10.4167  6.9444  10.4167]
%           percentageShift = 100*[Fc1, Fc2, Fc3]/Fs = [-7.2, 0, 7.2]
%
%------------------------------------------------
% Initialization
%------------------------------------------------
NN = size(SAR,1); % range samples
MM = size(SAR,2); % along-track samples
azFilteredRGB = zeros(NN,MM,3, 'uint8');
% Normalize input SAR image
SAR = SAR/max(max(abs(SAR)));
% If these input parameters are scalars, convert to [1x3] vector
if length(ratioAz) == 1
    ratioAz = ratioAz*ones(1,3);
end
if length(percentageShift) == 1
    percentageShift = percentageShift*[-1 0 1]; 
end
% First and last samples of band rejection when the template filter is
%   centered at 0. The first sample of the filter is the frequency 0 Hz.
minZero = round(MM./(2*ratioAz)+1); % [samples]
maxZero = round(MM*(2*ratioAz-1)./(2*ratioAz)); % [samples]
%------------------------------------------------
% Get each doppler band
%------------------------------------------------
for channel = 1:3 % shift for negative, zero and positive
    shiftAz = percentageShift(channel);
    % Template of doppler filter before shifting the filter
    doppFilter1 = ones(1,MM);
    doppFilter1(minZero(channel):maxZero(channel)) = 0;
    % Shift template filter
    doppFilter2 = circshift(doppFilter1,[0 round((shiftAz/100)*MM)]);
    % a. Filtering
    azFilteredRGBfloat = ifft(fft(SAR.').*repmat(doppFilter2.',[1 NN])).';
    % b. Limits for the filtered image (sub-band)
    % b1. Minimum and maximum limits
    sar_min_db = 20*min(min(log10(abs(azFilteredRGBfloat)))); % [dB]
    if equalization_flag
        sar_max_db = 20*max(max(log10(abs(azFilteredRGBfloat)))); % [dB]
    else
        sar_max_db = 0; % [dB]
    end
    % b2. Upper and lower limits for the scaling
    rgb_max_db = max(sar_min_db, sar_max_db - dBLim(1)); % [dB]
    rgb_min_db = max(sar_min_db, sar_max_db - dBLim(2)); % [dB]
    % c. Scaling
    azFilteredRGBfloat = 255*(20*log10(abs(azFilteredRGBfloat)) - rgb_min_db)/(rgb_max_db - rgb_min_db); % [0 - 255]
    % d. Assign doppler channel
    azFilteredRGB(:,:,channel) = uint8(azFilteredRGBfloat);
end
clear 'azFilteredRGBfloat',
stringTitle = ['RGB doppler bands (Bw=' num2str(100./ratioAz) '%, Fc=' num2str(percentageShift) '%)'];
figure,image(xGrid, yGrid, azFilteredRGB),grid on,
xlabel('Azimuth (samples)'), ylabel('Range (samples)'), title(stringTitle)
set(gca,'FontSize',18), set(gcf,'Color','w')
%------------------------------------------------
% Image with the intensity of the most intense band in each SAR sample
%------------------------------------------------
azFilteredRGBunique = zeros(NN,MM,3, 'uint8');
% Get most intense band
[~,indColourMax] = max(azFilteredRGB,[],3); 
% Index in 3D matrix
indexMax = (1:NN*MM).'+(indColourMax(:)-1)*NN*MM;
% Assign intensity
azFilteredRGBunique(indexMax) = azFilteredRGB(indexMax);
clear 'indexMax'
indColourMax = uint8(indColourMax);
% Represents
stringTitle = ['Maximum intensity of RGB doppler bands (Bw=' num2str(100./ratioAz) '%, Fc=' num2str(percentageShift) '%)'];
figure,image(xGrid, yGrid, azFilteredRGBunique),grid on
xlabel('Azimuth (samples)'), ylabel('Range (samples)'), title(stringTitle)
set(gca,'FontSize',18), set(gcf,'Color','w')
%------------------------------------------------
% Image with the the most intense band in each SAR sample
%------------------------------------------------
stringTitle = ['Most intense RGB doppler band (Bw=' num2str(100./ratioAz) '%, Fc=' num2str(percentageShift) '%)'];
figure,imagesc(xGrid, yGrid, indColourMax),grid on, colormap([1 0 0; 0 1 0; 0 0 1])
xlabel('Azimuth (samples)'), ylabel('Range (samples)'), title(stringTitle)
set(gca,'FontSize',18), set(gcf,'Color','w')
set(gca,'FontSize',18), set(gcf,'Color','w')