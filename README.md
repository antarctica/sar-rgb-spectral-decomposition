# sar-rgb-spectral-decomposition
Matlab software for generating the Doppler Decomposition of Synthetic Aperture Radar (SAR) images by splitting the full Doppler bandwidth in 3 sub-bands, assigning colours red, green and blue to the lower, medium, and upper doppler frequency bands. For a universal interpretation by colour-blind users, different colourmaps with a variation palette of triplets of colours are included, together with their estimated visualisations.

Cite the code: [![DOI](https://zenodo.org/badge/667776596.svg)](https://zenodo.org/badge/latestdoi/667776596)

<p align="center">
  <img src="/sar_vs_rgb.png" alt="comparison of beamwidths of SAR and RGB Decomposition">
</p>

## Scope
The repository contains Matlab functions for generating the RGB SAR Doppler Decomposition, different colourmaps with triplets of colours, and an estimator of the visualisation by colour blind interpreters with protanopia (insensitivity to red light), deuteranopia (to green) and tritanopia (to blue).

(To be submitted) Alvaro Arenas-Pingarrón, Alex M. Brisbourne, Carlos Martín, Hugh F.J. Corr, Carl Robinson, Tom A. Jordan, Paul V. Brennan: 'An Alternative Representation of Synthetic Aperture Radar Images for Englacial Observations', EGU The Cryosphere.

The dataset for testing the algorithm is composed of two SAR images from the [British Antarctic Survey](https://www.bas.ac.uk) (BAS) ice-sounding airborne Synthetic Aperture Radar (SAR) PASIN2 (Polarimetric Airborne Scientific INstrument, mark 2). PASIN2 is an active 150MHz pulsed airborne SAR (Synthetic Aperture Radar) for shallow, englacial and deep subglacial imagery. The SAR products are formatted as [NetCDF](https://www.unidata.ucar.edu/software/netcdf/) (Network Common Data Format) files (.nc) for data sharing. They contain standard data and its associated metadata, including units, conventions, intervals and human-readable descriptions, among others. NetCDF file content can be explored and plotted with software such as the free stand-alone Panoply, and the content can be read and created using available interfaces in C, C++, Java, Fortran, Python, IDL, MATLAB, and R, among others. The NetCDF files are archived by the [UK Polar Data Centre](https://www.bas.ac.uk/data/uk-pdc/) (UK PDC)

Arenas Pingarron, A., Brisbourne, A., Corr, H., Jordan, T., Robinson, C., Martin, C., Nicholls, K., & Smith, A. (2023). Ice-sounding airborne synthetic aperture radar depth profiles from Recovery Ice Stream 2016/17 and Rutford Ice Strem 2019/20 to test the RGB-Doppler-Decomposition method. (Version 1.0) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/40c2f86b-1a02-4106-934a-42769682df66

## Objectives
The users can:	
1.	Obtain the SAR image parameters indicated for the SAR Doppler Decomposition.
2.	Create different colourmaps.
3. Estimate the visualisation by colour blind interpreters.

## Outcomes
The users will:
1.	Recognise how to generate the SAR Doppler Decomposition.
2.	Identify the optimum colourmap and representation of the decomposition, depending on the interpreters.

## Content
- `rgb_doppler_decomposition_dB_min_max_v1.m`: Matlab function to generate the RGB SAR Doppler Decomposition.
- `sar_rgb_doppler_decomposition_dataset\`: folder with the two parameter files in Matlab language (.m) to process the dataset test with the SAR-product images.
- `sar_rgb_doppler_decomposition_scripts\`: folder with Matlab scripts to calculate the RGB SAR Doppler Decompositions from the SAR products (in NetCDF format) of the dataset, and represent them with different colourmaps.
- `colourmaps\`:  folder with the software related to the colourmap generation and visualisation.
- `colourmaps\universally_readable_colourmap.m`: Matlab function to generate colourmaps for the SAR Decomposition. Attempt to build an universally readable colourmap by using three colours whose combinations can be distinguished for colour-blindness due to protanopia (lack of sensitivity to red), deuteranopia (green) and tritanopia (blue). The function creates three palettes of shades for a triplet of colours, each colour to represent each of the three bands of the Decomposition. The aim of the colourmap is that not only each of the three basic colours can be distinguished by colour blind interpreters, but also the combinations of the palettes, which will occurr depending on the Decomposition.
- `colourmaps\rgb_to_colour_blindness.m`: Matlab function to calculate the RGB colours as seen by  colour-blind interpreters with protanopia, deuteranopia and tritanopia.

## RGB SAR Doppler Decomposition
 ```MATLAB
 azFilteredRGB = rgb_doppler_decomposition_dB_min_max_v1(SAR, ratioAz, percentageShift, xGrid, yGrid, dBLim, equalization_flag)
```

Obtains a coloured RGB-SAR from a `SAR` image in the along-track grid `xGrid` and range grid `yGrid`, by splitting the full doppler bandwidth in 3 sub-bands, asigning colours red, green and blue to the lower, medium, and upper doppler frequency bands. Each band has a width defined by `ratioAz` over the full bandwidth, and the central frequency for each band is defined by `percentageShift`, a percentage over the full bandwidth. Each band has a colour defined by the RGB colour model, with pure red, green and blue colours, using 8 bits per colour band, with integer numbers between 0 and 255 levels:

1. lower sub-band:  0 <= R <= 255
2. medium sub-band: 0 <= G <= 255
3. upper sub-band:  0 <= B <= 255

The limits to set the darkest (0) and brightest (255) levels are given by the parameter `dBLim`>=0, in dB units below the maximum of the SAR image (when the input `equalization_flag` is false), or each sub-band (when the input `equalization_flag` is true).

Steps:
  1. The image is normalized to its maximum
  2. For each of the 3 sub-bands:
  	a. Filtering: filter the image for the sub-band
  	b. Limit the filtered image to amplitudes within the upper limit  `-dBLim(1)` and lower limit `-dBLim(2)`, so that `-dBLim(2) <= filtered_image <= -dBLim(1)`
  	c. Scaling: convert to integer numbers between 0 and 255
  	d. Assign to the sub-band
  3. Calculate a new image with the intensity (0-255) of the most intense sub-band. This image is displayed, but it is not part of the output.
  4. Calculate a new image with the maximum intensity (255) of the most intense sub-band. This image is displayed, but it is not part of the output.

### Example parameters
```MATLAB
xGrid = 0:size(complex_SARimage,1)-1;
yGrid= 0:size(complex_SARimage,2)-1;
dBLim = [10 90]; % (dB) limit below the maximum: upper=max-10; lower=max-90
equalization_flag = true;
```
The sampling frequency `Fs` of the domain of interest is always greater than the band of interest `Bw` of the data, so that `Fs > Bw`. If the Doppler band of interest is
```MATLAB
% Bandwidth of interest
Bw = 30; % [Hz]
```
then the band to cover with the 3 sub-bands of decomposition is the interval `[-Bw/2, +Bw/2] = [-15, +15]`, with non-overlapping or overlapping bands.

#### Non-overlapping bands:
<p align="center">
  <img src="/non_overlapping_bands.svg" width="230" height="150" alt="spectrum without overlapping bands">
</p>

For non-overlapping bands seamlessly covering the interval [-Bw/2, +Bw/2], the three sub-bandwidths would be `SBw1 = SBw2 = SBw3 = Bw/3 = 10 Hz`, and the central frequencies `[Fc1, Fc2, Fc3] = [-Bw/3, 0, +Bw/3] = [-10, 0, +10] Hz`.
For a sampling frequency `Fs = 62.5 Hz`, the total available band would be `[-Fs/2, +Fs/2]`. The input parameters `ratioAz` and `percentageShift` would be:
```MATLAB
% Non-overlapping with Fs = 62.5 Hz and Bw = 30 Hz
% Sampling frequency
Fs = 62.5; % [Hz]
% Subband-widths
SBw = Bw/3; % [Hz]
% Central frequencies
Fc1 = -Bw/3; % [Hz]
Fc2 = 0; % [Hz]
Fc3 = Bw/3; % [Hz]
% Doppler decomposition parameters
ratioAz = Fs./[SBw, SBw, SBw]; % [6.25 6.25 6.25] [-]
percentageShift = 100*[Fc1, Fc2, Fc3]/Fs; % [-16, 0, 16] [%]
% Doppler decomposition
azFilteredRGB = rgb_doppler_decomposition_dB_min_max_v1(complex_SARimage, ratioAz, percentageShift, xGrid, yGrid, dBLim, equalization_flag);
```
so that `ratioAz > 1` and `-100 < percentageShift < 100`

For a sampling frequency `Fs = 125 Hz`, the input parameters would be:
```MATLAB
% Non-overlapping with Fs = 125 Hz and Bw = 30 Hz
% Sampling frequency
Fs = 125; % [Hz]
% Subband-widths
SBw = Bw/3; % [Hz]
% Central frequencies
Fc1 = -Bw/3; % [Hz]
Fc2 = 0; % [Hz]
Fc3 = Bw/3; % [Hz]
% Doppler decomposition parameters
ratioAz = Fs./[SBw, SBw, SBw]; % [12 12 12] [-]
percentageShift = 100*[Fc1, Fc2, Fc3]/Fs; % [-8, 0, 8] [%]
```

#### Overlapping bands
<p align="center">
  <img src="/overlapping_bands.svg" width="230" height="150" alt="spectrum with overlapping bands">
</p>

For overlapping sub-bands seamlessly covering the interval `[-Bw/2, +Bw/2]`, in a way that the overlapped intervals have the same width as the intervals without overlapping, the three sub-bandwidths would be  `SBw1 = SBw3 = 0.4Bw = 12 Hz` and `SBw2 = 0.6Bw = 18 Hz`, and the central frequencies `[Fc1, Fc2, Fc3] = [-0.3Bw, 0, +0.3Bs] = [-9, 0, +9] Hz`
For a sampling frequency `Fs = 62.5 Hz`, the total available band would be `[-Fs/2, +Fs/2]`. The input parameters `ratioAz` and `percentageShift` would be:
          ratioAz = Fs./[SBw1, SBw2, SBw3] = [5.2083  3.4722  5.2083]
          percentageShift = 100*[Fc1, Fc2, Fc3]/Fs = [-14.4  0  14.4]
```MATLAB
% Overlapping with Fs = 62.5 Hz and Bw = 30 Hz
% Sampling frequency
Fs = 62.5; % [Hz]
% Subband-widths
SBw1 = 0.4*Bw; % [Hz]
SBw2 = 0.6*Bw; % [Hz]
SBw3 = 0.4*Bw; % [Hz]
% Central frequencies
Fc1 = -0.3*Bw; % [Hz]
Fc2 = 0; % [Hz]
Fc3 = 0.3*Bw; % [Hz]
% Doppler decomposition parameters
ratioAz = Fs./[SBw1, SBw2, SBw3]; % [5.2083  3.4722  5.2083]  [-]
percentageShift = 100*[Fc1, Fc2, Fc3]/Fs; % [-14.4  0  14.4] [%]
```
so that `ratioAz > 1` and `-100 < percentageShift < 100`

For a sampling frequency `Fs = 125 Hz`, the input parameters would be:
```MATLAB
% Overlapping with Fs = 125 Hz and Bw = 30 Hz
% Sampling frequency
Fs =125; % [Hz]
% Subband-widths
SBw1 = 0.4*Bw; % [Hz]
SBw2 = 0.6*Bw; % [Hz]
SBw3 = 0.4*Bw; % [Hz]
% Central frequencies
Fc1 = -0.3*Bw; % [Hz]
Fc2 = 0; % [Hz]
Fc3 = 0.3*Bw; % [Hz]
% Doppler decomposition parameters
ratioAz = Fs./[SBw1, SBw2, SBw3]; % [10.4167  6.9444  10.4167] [-]
percentageShift = 100*[Fc1, Fc2, Fc3]/Fs; % [-7.2, 0, 7.2] [%]
```
## Universally readable colourmaps
```MATLAB
function [readable_colourmap,...
    first_band_template,...
    second_band_template,...
    third_band_template] = universally_readable_colourmap(...
    palette_code,...
    blind_conversion_method,...
    wanted_colours_per_band,...
    kovesi_modification_flag,...
    display_single_palette_flag,...
    display_multiple_palette_flag)
```
<p align="center">
  <img src="/colourmaps/colourmaps.png" alt="colourmaps and interpretation for colour-blind interpreters">
</p>

Attempt to build an universally readable colourmap for RGB SAR Doppler Decomposition images, by using three colours whose combinations can be distinguished for colour-blindness due to protanopia (lack of sensitivity to red), deuteranopia (green) and tritanopia (blue).

The function creates templates with the variation palette of a triplet of colours according to the input code `palette_code`, each colour to represent each of the three bands of a Frequency Decomposition analysis. The palette of each template has 8 colour shades related to the intensity of the frequency band: the maximum intensity of each band is represented by the primary colour of that band, and the minimum intensity by a dark version of this primary colour.

The output `readable_colourmap` is obtained from interpolating `wanted_colours_per_band` from each of the three shade templates. The aim of the colourmap is that not only each of the three primary colours can be distinguished by colour-blind interpreters, but also the combinations of the palettes, which will occurr depending on the Frequency Decomposition. The colours of `readable_colourmap` may not coincide with the colours from the output templates `first_band_template`, etc.

1.   Frequency band 1 -> primary colour 1 ->  output `first_band_template` (primary colour to dark) -> from the template, interpolate `wanted_colours_per_band` colours to assign to the entries `2` to `wanted_colours_per_band+1` in output `readable_colourmap`
2.   Frequency band 2 -> primary colour 2 ->  output `second_band_template` (primary colour to dark) -> from the template, interpolate `wanted_colours_per_band` colours to assign to the entries `wanted_colours_per_band+2` to `2*wanted_colours_per_band+1` in output `readable_colourmap`
3.   Frequency band 3 -> primary colour 3 ->  output `third_band_template` (primary colour to dark) -> from the template, interpolate `wanted_colours_per_band` colours to assign to the entries `2*wanted_colours_per_band+2` to `3*wanted_colours_per_band+1` in output `readable_colourmap`

The output `readable_colourmap` is useful for comparing the colours when scaled depending on the Frequency Decomposition. The outputs from this function can be used in two ways:
1. the three primary colours (the first entry from the template outputs `first_band_template`, `second_band_template` and `third_band_template`) can be used directly to scale the outputs of each band of the Frequency Decomposition,
2. or the `wanted_colours_per_band` interpolated colours can be assigned to `wanted_colours_per_band` output bins of each band of the Frequency Decomposition.

Possible triplets from `palette_code`, with RGB integer values between 0 and 255:
+ `palette_code=0`: Pure Red, Green and Blue (R-G-B):
	+ R = RGB(1,0,0) x 255
	+ G = RGB(0,1,0) x 255
	+ B = RGB(0,0,1) x 255
+ `palette_code=-1`: Modified Red, Green and Blue (Rm-Gm-Bm, `kovesi_modification_flag=0`):
	+ Rm = RGB(0.9, 0.0, 0.0) x 255
	+ Gm = RGB(0.0, 0.8, 0.0) x 255
	+ Bm = RGB(0.1, 0.2, 1.0) x 255
+ `palette_code=-2`: Olive, Teal and Purple (O-T-P):
	+ O = RGB(0.5, 0.5, 0) x 255
	+ T = RGB(0.0, 0.5, 0.5) x 255
	+ P= RGB(0.5, 0.0, 0.5) x 255
+ `palette_code=-3`: Interpolations between Rm-Gm-Bm and O-T-P:
+ `palette_code=+1`: Universally readable v1, whose total summation of colours is NOT equal to pure white.
+ `palette_code=+2`: Universally readable v2, whose total summation of colours is NOT equal to pure white.
+ `palette_code=+3`: Universally readable v3 with dark yellow, dark grey and violet blue (Yd-Gd-Bv), whose total summation of colours is equal to pure white, with :
	+ Yd = RGB(0.55, 0.55, 0) x 255
	+ Gd = RGB(0.25, 0.25, 0.25) x 255
	+ Bv = RGB(0.20, 0.20, 0.75) x 255

To estimate the visualisation by color-blind interpreters, the function `colourmaps\rgb_to_colour_blindness.m` is called. The algorithm for this estimation is chosen from the input parameter `blind_conversion_method`. 

If The flag `kovesi_modification_flag` is active it modifies the pure red, green and blue to make them have uniform lightness (see [Kovesi, P. Good colour maps: how to design them. CoRR abs/1509.03700 (2015)](https://arxiv.org/abs/1509.03700v1 "Kovesi, P. Good colour maps: how to design them. CoRR abs/1509.03700 (2015)") or the [Peter Kovesi repository](https://github.com/peterkovesi/PerceptualColourMaps.jl)), and deactivate the flag to adopt other modification with less uniform lightness.

To display colourmaps with the indivudual interpolated palettes as expected to be seen according to each blindness type, activate the input boolean flag `display_single_palette_flag`. To display in one image the colourmaps of the single interpolated palette for each blindness type, and in another image the mixing-circles with the colour combinations of the first colour from the output templates, activate the input boolean flag `display_multiple_palette_flag`.

### Examples
##### Triplet with modified R-G-B primary colours
```MATLAB
% Triplet with modified Red-Green-Blue as primary colours
% Palette code
palette_code = -1;
% Conversion algorithm
blind_conversion_method = 4; % the recommended is the fourth
% Number of colours from primary to darkest
wanted_colours_per_band = 8;
% Uniform lightness of modified Red-Green-Blue, or less uniform
kovesi_modification_flag = false; % less uniform
% Display options
display_single_palette_flag = true;
display_multiple_palette_flag = true;
% Run
[readable_colourmap,...
	first_band_template,...
	second_band_template,...
	third_band_template] = universally_readable_colourmap(...
	palette_code,...
	blind_conversion_method,...
	wanted_colours_per_band,...
	kovesi_modification_flag,...
	display_single_palette_flag,...
	display_multiple_palette_flag)
```
##### Get the three primary colours
From the templates of each of the three sub-bands, take the first colour:
```MATLAB
% Primary colours, the first of each template
primary_colours = zeros(3, 3);
primary_colours(1,:) = first_band_template(1,:);  % [0...255]
primary_colours(2,:) = second_band_template(1,:);  % [0...255]
primary_colours(3,:) = third_band_template(1,:);  % [0...255]
```
##### Apply primary colours to Doppler Decomposition
From the Doppler Decomposition, the output `azFilteredRGB` is a variable of three dimensions: range (SAR dimension), along-track (SAR dimension) and sub-band (three Doppler sub-bands). This output is normalised with values between 0 and 255. If the primary colours were pure red, green and blue, the coloured output would be the same as for the Doppler Decomposition, as it would exactly correspond with the pure colours. If the primary colours are not pure red, green and blue, a colour conversion must be done:
```MATLAB
% Convert primary colours into values from 0 to 1
colour_1st_band = primary_colours(1, : ) / 255; % [0...1]
colour_2nd_band = primary_colours(2, :) / 255; % [0...1]
colour_3rd_band = primary_colours(3, :) / 255; % [0...1]
% Initialise re-coloured image
azFilteredRGBmod = zeros(size(azFilteredRGB), 'uint8');
for idx_band = 1:3
	% Modified RGB
	azFilteredRGBmod(:, :, idx_band) = ...
	colour_1st_band(idx_band)*azFilteredRGB(:, :, 1) + ...
	colour_2nd_band(idx_band)*azFilteredRGB(:, :, 2) + ...
	colour_3rd_band(idx_band)*azFilteredRGB(:, :, 3); % [0...255]
end
```

### References

## Estimation of visualisation for colour-blind interpreters
```MATLAB
function rgb_256_blind = rgb_to_colour_blindness(rgb_256_matrix, blindness_type, conversion_method)
```
Estimates the RGB colours of the the input `rgb_256_matrix` as seen by the colour-blind interpreters with total either protanopia (lack of sensitivity to red), deuteranopia (to green) or tritanopia (to blue) depending on the input `blindness_type`. There are four estimation methods, one chosen from the parameter `conversion_method`.

[Nicolas Burrus](http://nicolas.burrus.name/) provides a comprehensive [description of algorithms](https://daltonlens.org/understanding-cvd-simulation/) for simulating colour vision anomalies, their [comparison and scope](https://daltonlens.org/opensource-cvd-simulation/), an interactive [simulator](https://daltonlens.org/colorblindness-simulator), and Python code in their [repository](https://github.com/DaltonLens/DaltonLens-Python).

The four selectable methods are:
1. **Method 1** `conversion_method = 1`,  explained by Martin Krzywinski in his [website](http://mkweb.bcgsc.ca/), in [Designing For Color Blindness](http://mkweb.bcgsc.ca/colorblind/math.mhtml#projecthome).
2. **Method 2** `conversion_method = 2`, explained in the [authors repository](https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html):
G. M. Machado, M. M. Oliveira and L. A. F. Fernandes, "[A Physiologically-based Model for Simulation of Color Vision Deficiency](https://doi.org/10.1109/TVCG.2009.113)," in IEEE Transactions on Visualization and Computer Graphics, vol. 15, no. 6, pp. 1291-1298, Nov.-Dec. 2009, https://doi.org/10.1109/TVCG.2009.113.
3. **Method 3**, `conversion_method = 3`:
Brettel, H., Vienot, F., & Mollon, J. D. (1997). [Computerized simulation of color appearance for dichromats](http://vision.psychol.cam.ac.uk/jdmollon/papers/Dichromatsimulation.pdf). Journal of the Optical Society of America. A, Optics, Image Science, and  Vision, 14(10), 2647–2655. https://doi.org/10.1364/josaa.14.002647
3. **Method 4**, `conversion_method = 4`, **our preferred method** because the article is open access, explained [here](https://daltonlens.org/understanding-cvd-simulation/#The-simplification-of-Vi%C3%A9not-1999-for-protanopes-and-deuteranopes "here") in the website https://daltonlens.org/:
Viénot, F., Brettel, H. and Mollon, J.D. (1999), [Digital video colourmaps for checking the legibility of displays by dichromats](https://doi.org/10.1002/(SICI)1520-6378(199908)24:4%3C243::AID-COL5%3E3.0.CO;2-3). Color Res. Appl., 24: 243-252. http<area>s://doi.org/10.1002/(SICI)1520-6378(199908)24:4<243::AID-COL5>3.0.CO;2-3.

### Example

##### Get the three primary colours for interpreters with protanopia
From the outputs of calling the function `colourmaps\universally_readable_colourmap.m` take the first colour from the templates of each of the three sub-bands, and run:

```MATLAB
% Blindness type
blindness_type = 'protanopia'; % 'P', 'p', 'pro', etc
% Conversion algorithm
blind_conversion_method = 4; % the recommended is the fourth
% Colour interpretation
new_colours = uint8(rgb_to_colour_blindness(primary_colours, blindness_type, blind_conversion_method));  % [0...255]
```
### References
##### General on colour blindness:
1. https://daltonlens.org/understanding-cvd-simulation/
2. http://www.vischeck.com/
3. http://www.daltonize.org/
4. https://daltonlens.org/opensource-cvd-simulation/
5. https://daltonlens.org/colorblindness-simulator

##### Source code
1. https://github.com/DaltonLens/DaltonLens-Python
2. https://github.com/joergdietrich/daltonize/tree/main
3. https://github.com/MaPePeR/jsColorblindSimulator

##### Visual references:
1. [DaltonLens](https://daltonlens.org/opensource-cvd-simulation/): a great number of displays according to the several algorithms.
2. [ConvertingColors](https://convertingcolors.com/rgb-color-255_0_0.html): Coblis Version2 (see https://daltonlens.org/opensource-cvd-simulation/)
3. [Coblis v2](https://www.color-blindness.com/coblis-color-blindness-simulator/): Coblis Version2


## LICENSE
MIT license (https://github.com/antarctica/sar-rgb-spectral-decomposition/blob/main/LICENSE)
