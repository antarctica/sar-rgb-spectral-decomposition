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
%
% Attempt to build an universally readable colourmap for RGB SAR images,
%   by using three colours whose combinations can be distinguished for
%   colour-blindness due to protanopia (lack of sensitivity to red),
%   deuteranopia (green) and tritanopia (blue).
%
% The function creates templates 'first_band_template',
%   'second_band_template' and 'third_band_template' with the variation
%   palette of a triplet of RGB colours according to 'palette_code', each 
%   colour to represent each of the three bands of a Frequency 
%   Decomposition analysis. The palette of each band has 8 colour shades 
%   related to the  intensity of the frequency band: the maximum intensity 
%   of each band is represented by the primary colour of that band, and the
%    minimum intensity by a dark version of this primary colour.
%
% The output 'readable_colourmap' is obtained from interpolating 
%   'wanted_colours_per_band' from each of the three shade templates. The 
%   aim of the colourmap is that not only each of the three primary colours 
%   can be distinguished by colour-blind interpreters, but also the 
%   combinations of the palettes, which will occurr depending on the 
%   Frequency Decomposition. The colours of 'readable_colourmap' may not 
%   coincide with the colours from the output templates
%   'first_band_template', etc.
%
%   Frequency band 1 -> primary colour 1 -> output 'first_band_template'
%       (primary colour to dark) -> from the template, interpolate
%       'wanted_colours_per_band' colours to assign to the entries '2' to
%       'wanted_colours_per_band+1' in the output 'readable_colourmap'.
%   Frequency band 2 -> primary colour 2 -> output 'second_band_template'
%       (primary colour to dark) -> from the template, interpolate
%       'wanted_colours_per_band' colours to assign to the entries
%       'wanted_colours_per_band+2' to '2*wanted_colours_per_band+1' in
%       the output 'readable_colourmap'.
%   Frequency band 3 -> primary colour 3 -> output 'third_band_template'
%       (primary colour to dark) -> from the template, interpolate
%       'wanted_colours_per_band' colours to assign to the entries
%       '2*wanted_colours_per_band+2' to '3*wanted_colours_per_band+1' in
%       the output 'readable_colourmap'.
%
% To estimate the visualisation by color-blind interpreters, the function
%   'rgb_to_colour_blindness.m is called'. The algorithm for this
%   estimation is chosen from the input parameter 'blind_conversion_method'.
%
% As reference for the colours, the next websites can be used:
%   (not preferred, based on Coblis Version2) www.convertingcolors.com
%   (preferred) https://www.color-blindness.com/coblis-color-blindness-simulator/
%   (preferred) https://daltonlens.org/colorblindness-simulator
%
% The primary colours described by their Red-Green-Blue components, each
%   defined with an integer between 0 and 255, are combinations and
%   variatons of the next:
%   RGB Primary colours
%       red = RGB(255,0,0)
%       green = RGB(0,255,0)
%       blue = RGB(0,0,255)
%   OTP primary colours
%       olive = RGB(128,128,0)
%       teal = RGB(0,128,128)
%       purple = RGB(128,0,128)
%
% The output 'readable_colourmap' is useful for comparing the colours when
%   scaled depending on the Frequency Decomposition. The outputs from this
%   function can be used in two ways:
%       1. the three primary colours (the first entry from the template
%           outputs 'first_band_template', 'second_band_template' and
%           'third_band_template') can be used directly to scale the
%           outputs of each band of the Frequency Decomposition,
%       2. or the 'wanted_colours_per_band' interpolated colours can be
%           assigned to 'wanted_colours_per_band' output bins of each band
%           of the Frequency Decomposition.
%
% Input:
%     palette_code (integer): number to select each of the three primary
%           (representing) colours of the triplet and the palette extended 
%           from each of these colours. There is a palette (range of 
%           colours) for each of the three primary colours of the triplet
%           (of the three bands of a Frequency Decomposition analysis). 
%           Each palette is a range of colours interpolated from the 
%           primary colour to its dark version.
%          -3: Interpolation of secondary (olive, teal, purple) and
%              modified-primary (red, green, blue) colours, with the total
%              summation of colours equals to pure white ([1,1,1]). The
%              weight of the secondary colours is 75% and the weight of the
%              primary 25%. These weights are given by the fixed parameter:
%              secondary_weight = 0.75;
%          -2: Secondary colours (olive, teal, purple), with the total
%              summation of colours equals to pure white ([1,1,1]).
%          -1: Modified Red, Green, Blue for uniform perception, still with
%              the total summation of colours equals to pure white ([1,1,1]).
%              Reference: Kovesi, P. Good colour maps: how to design them. CoRR abs/1509.03700 (2015).
%           0: Pure Red, Green, Blue
%           1: Universally readable, v1:
%              [goldenish, light brownish, violetish-bluish], the total
%              summation of colours is NOT equal to pure white.
%           2: Universally readable, v2:
%              [goldenish, light brownish, violetish-bluish], the total
%              summation of colours is NOT equal to pure white.
%           3: Universally readable, v3:
%              [yellowish, dark_grey, bluish], the total
%              summation of colours equals to pure white ([1,1,1]).
%           Example: palette_code = -3;
%     blind_conversion_method (positive integer): algorithm number to
%           convert colours to colour-blindess equivalent.
%           1: Uses LMSD65 to convert from XYZ to LMS (normalized to D65),
%              but there are other options such as CIECAM97 and CIECAM02.
%              http://mkweb.bcgsc.ca/colorblind/math.mhtml#projecthome
%              http://mkweb.bcgsc.ca/colorblind/code.mhtml
%           2: Machado 2009
%              https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html
%              https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/Machado_Oliveira_Fernandes_CVD_Vis2009_final.pdf
%              Reference:
%               G. M. Machado, M. M. Oliveira and L. A. F. Fernandes, "A
%               Physiologically-based Model for Simulation of Color Vision
%               Deficiency," in IEEE Transactions on Visualization and Computer
%               Graphics, vol. 15, no. 6, pp. 1291-1298, Nov.-Dec. 2009,
%               doi: 10.1109/TVCG.2009.113.
%           3: Brettel 1997. The Vischeck display filter for GIMP.
%               Brettel, H., Vienot, F., & Mollon, J. D. (1997). Computerized
%                   simulation of color appearance for dichromats. Journal of the
%                   Optical Society of America. A, Optics, Image Science, and
%                   Vision, 14(10), 2647–2655.
%                   https://doi.org/10.1364/josaa.14.002647
%              https://github.com/GNOME/gimp/blob/master/modules/display-filter-color-blind.c
%              https://daltonlens.org/understanding-cvd-simulation/#Comparing-Vi%C3%A9not-1999-and-Brettel-1997
%              http://www.vischeck.com/
%           4: Based on Vie'not 1999, but roughly
%              Reference:
%               Viénot, F., Brettel, H. and Mollon, J.D. (1999), Digital video
%                   colourmaps for checking the legibility of displays by
%                   dichromats. Color Res. Appl., 24: 243-252.
%                   https://doi.org/10.1002/(SICI)1520-6378(199908)24:4<243::AID-COL5>3.0.CO;2-3
%               https://daltonlens.org/understanding-cvd-simulation/#The-simplification-of-Vi%C3%A9not-1999-for-protanopes-and-deuteranopes
%               https://github.com/joergdietrich/daltonize/tree/main/doc
%             Example: blind_conversion_method = 4;
%     wanted_colours_per_band (postive integer): numberof colours per band.
%             Example: wanted_colours_per_band = 8;
%     kovesi_modification_flag(logical): flag for modifying the pure red,
%           green and blue to make them have uniform lightness, according
%           to the reference article
%               Kovesi, P. Good colour maps: how to design them. CoRR abs/1509.03700 (2015).
%               https://arxiv.org/abs/1509.03700
%           To follow Kovesi recommendation, activate the flag (true); to
%           adopt other modification with less uniform lightness,
%           deactivate the flag (false).
%             Example: kovesi_modification_flag = false;
%     display_single_palette_flag(logical): flag for displaying (true) or
%           not (false) the colourmaps with the indivudual palettes as
%           expected to be seen according to the colour blindness.
%             Example: display_single_palette_flag = false;
%     display_single_palette_flag(logical): flag to display (true) or
%           not (false) the colourmaps with the indivudual interpolated 
%           palettes as expected to be seen according to each blindness type.
%             Example: display_single_palette_flag = false;
%     display_multiple_palette_flag(logical): flag to display (true) or
%           not (false) in one image the colormaps of the single
%           interpolated palette for each blindness type, and in another 
%           image the mixing-circles with the colour combinations of the
%           first colour from the output templates.
%             Example: display_multiple_palette_flag = false;
%
% Output:
%     readable_colourmap (uint8): [N x 3] matrix with the RGB colours of
%           the colourmap, with N = 1 + 3*wanted_colours_per_band. The
%           first colour is always black RGB(0,0,0), followed by
%           'wanted_colours_per_band' of the palette interpolated from the
%           first band template, 'wanted_colours_per_band' of the palette
%           interpolated from the second band template, and 
%           'wanted_colours_per_band' of the palette interpolated from the
%           third band template. Values between 0 and 255.
%             Example: readable_colourmap = [0 0 0; 148 93 0; ...]
%     first_band_template (uint8): [F x 3] matrix with the 'F' RGB colours 
%           of the template to be interpolated for obtaining the palette of
%           'wanted_colours_per_band' of the 1st band. Values between 0 and 255.
%             Example: first_band_template = [255 0 0; 200 0 0; ...]
%     second_band_template (uint8): [S x 3] matrix with the 'S' RGB colours
%           of the template to be interpolated for obtaining the palette of
%           'wanted_colours_per_band' of the 2nd band. Values between 0 and 255.
%             Example: second_band_template = [0 255 0; 0 200 0; ...]
%     third_band_template (uint8): [T x 3] matrix with the 'T' RGB colours 
%           of the template to be interpolated for obtaining the palette of
%           'wanted_colours_per_band' of the 3rd band. Values between 0 and 255.
%             Example: third_band_template = [0 0 255; 0 0 200; ...]

%
% Dependencies:
%         Called by: -
%         Calls: rgb_to_colour_blindness.m
%
% Date: 16.02.2023
%
% Author: Alvaro Arenas-Pingarron, British Antarctic Survey

% Beside the colour, to support the differentiation, in the second colour,
%   which is lighter in our version, will have gridded pattern with
%   horizontal and vertical lines. This is because the first (goldenish)
%   and second colour (light brownish) are less different than the third
%   (violetish-bluish).

% Colour definitions
if kovesi_modification_flag
    % Modified R-G-B for uniform lightness
    modified_red = [0.90 0.17 0.00];
    modified_green = [0.00 0.50 0.00];
    modified_blue = [0.10 0.33 1.00];
else
    % Not so uniform lightness
    modified_red = [0.90 0.00 0.00];
    modified_green = [0.00 0.80 0.00];
    modified_blue = [0.10 0.20 1.00];
    % modified_red = [0.90 0.10 0.00];
    % modified_green = [0.00 0.70 0.00];
    % modified_blue = [0.10 0.20 1.00];
end
olive_colour = [0.50 0.50 0.00];
teal_colour = [0.00 0.50 0.50];
purple_colour = [0.50 0.00 0.50];
% yellowish = [0.75 0.75 0];
% dark_grey = [0.25 0.25 0.35];
% bluish = [0 0 0.65];
yellowish = [0.75 0.60 0];
dark_grey = [0.25 0.30 0.35];
bluish = [0 0.10 0.65];
% yellowish = [1 1 0];
% dark_grey = [0 0 0];
% bluish = [0 0 1];
yellowish = [0.55 0.55 0];
dark_grey = [0.25 0.25 0.25];
bluish = [0.2 0.20 0.75];
%------------------------------------------------------------
% Templates for colour maps
%------------------------------------------------------------
% First band of colour triplet
switch palette_code
    case -3 % Of triplet with secondary colours: mod[Red, Green, Blue]
        % Averaged olive(128,128,0) and modified red
        secondary_weight = 0.75;
        interpolated_olive = secondary_weight*olive_colour + (1-secondary_weight)*modified_red;
        first_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(interpolated_olive);
    case -2 % Of triplet with secondary colours: mod[Red, Green, Blue]
        % Olive(128,128,0)
        first_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(olive_colour);
    case -1 % Of triplet with modified R-G-B for uniform perception: mod[Red, Green, Blue]
        % Modified red
        first_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(modified_red);
    case 0 % Of triplet with pure R-G-B: [Red, Green, Blue]
        % Pure red
        first_band_template = zeros(8,3);
        first_band_template(:,1) = 255*linspace(1, 0.1, size(first_band_template,1));
    case 1 % Of universally readable triplet, v1: [goldenish, light brownish, violetish-bluish]
        % Goldenish
        first_band_template = [229 205 84;... % +20%
                               199 178 56;... % +10%
                               170 151 25;... % base colour
                               141 125 0;...  % -10%
                               114 101 0;...  % -20%
                                86  77 0;...  % -30%
                                60  55 0];    % -40%
    case 2 % Of universally readable triplet, v2: [goldenish, light brownish, violetish-bluish]
        % Goldenish
        first_band_template = [240 201 84;... % +20%
                               210 173 57;... % +10%
                               180 147 27;... % base colour
                               151 121 0;...  % -10%
                               123  97 0;...  % -20%
                                95  73 0;...  % -30%
                                68  51 0];    % -40%
    case 3 % Of universally readable triplet, v3: [yellowish, dark_grey, bluish]
        first_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(yellowish);

end
% Second band of colour triplet
switch palette_code
    case -3 % Of triplet with secondary colours: mod[Red, Green, Blue]
        % Average teal(0,128,128) and modified green
        interpolated_teal = secondary_weight*teal_colour + (1-secondary_weight)*modified_green;
        second_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(interpolated_teal);
    case -2 % Of triplet with secondary colours: mod[Red, Green, Blue]
        % Teal(0,128,128)
        second_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(teal_colour);
    case -1 % Of triplet with modified R-G-B for uniform perception: mod[Red, Green, Blue]
        % Modified green
        second_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(modified_green);
    case 0 % Of triplet with pure R-G-B: [Red, Green, Blue]
        % Pure green
        second_band_template = zeros(8,3);
        second_band_template(:,2) = 255*linspace(1, 0.1, size(second_band_template,1));
    case 1 % Of universally readable triplet, v1: [goldenish, light brownish, violetish-bluish]
        % Light brownish
        second_band_template = [236 211 165;... % base colour
                               207 183 139;...  % -10%
                               179 157 113;...  % -20%
                               152 131  88;...  % -30%
                               125 105  64;...  % -40%
                                99  81  42];    % -50%
        %second_band_template(:,3) = second_band_template(:,3).*linspace(1,1.4,6).';
    case 2 % Of universally readable triplet, v2: [goldenish, light brownish, violetish-bluish]
        % Light brownish
        second_band_template = [240 211 165;... % base colour
                                211 183 138;... % -10%
                                183 157 113;... % -20%
                                155 131  88;... % -30%
                                129 105  64;... % -40%
                                103  81  41];   % -50%
    case 3 % Of universally readable triplet, v3: [yellowish, dark_grey, bluish]
        second_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(dark_grey);
end
% Third band of colour triplet
switch palette_code
    case -3 % Of triplet with secondary colours: mod[Red, Green, Blue]
        % Average purple (128,0,128) and modified blue
        interpolated_purple = secondary_weight*purple_colour + (1-secondary_weight)*modified_blue;
        third_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(interpolated_purple);
    case -2 % Of triplet with secondary colours: mod[Red, Green, Blue]
        % Purple colour (128,0,128)
        third_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(purple_colour);
    case -1 % Of triplet with modified R-G-B for uniform perception: mod[Red, Green, Blue]
        % Modified blue
        third_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(modified_blue);
    case 0 % Of triplet with pure R-G-B: [Red, Green, Blue]
        % Pure blue
        third_band_template = zeros(8,3);
        third_band_template(:,3) = 255*linspace(1, 0.1, size(third_band_template,1));
    case 1 % Of universally readable triplet, v1: [goldenish, light brownish, violetish-bluish]
        % Violetish-bluish
        third_band_template = [163 172 217;... % base colour
                              136 146 189;...  % -10%
                              110 120 162;...  % -20%
                               85  95 136;...  % -30%
                               60  72 110;...  % -40%
                               36  49  86];    % -50%
    case 2 % Of universally readable triplet, v2: [goldenish, light brownish, violetish-bluish]
        % Violetish-bluish
        third_band_template = [163 172 217;... % base colour
                              136 146 189;...  % -10%
                              110 120 162;...  % -20%
                               85  95 136;...  % -30%
                               60  72 110;...  % -40%
                               36  49  86];    % -50%
    case 3 % Of universally readable triplet, v3: [yellowish, dark_grey, bluish]
        third_band_template = 255*repmat(linspace(1, 0.1, 8).', 1, 3) * diag(bluish);
end
% Out of band colour
out_colour = [0 0 0];
%------------------------------------------------------------
% Interpolate band colours to the number of wanted colours
%------------------------------------------------------------
% First band
first_band_colours = zeros(wanted_colours_per_band, 3);
interpolation_interval = size(first_band_template,1) - 1;
if true
    first_start = 1.25; % start of interpolation
    first_end = 1 + interpolation_interval - 0.0*interpolation_interval; % end of interpolation
else
    first_start = 1 + 0.0*interpolation_interval; % start of interpolation
    first_end = 1 + interpolation_interval - 0.3*interpolation_interval; % end of interpolation
end
first_red_interpolator = griddedInterpolant(first_band_template(:,1));
first_green_interpolator = griddedInterpolant(first_band_template(:,2));
first_blue_interpolator = griddedInterpolant(first_band_template(:,3));
first_band_colours(:,1) = first_red_interpolator(linspace(first_start, first_end, wanted_colours_per_band));
first_band_colours(:,2) = first_green_interpolator(linspace(first_start, first_end, wanted_colours_per_band));
first_band_colours(:,3) = first_blue_interpolator(linspace(first_start, first_end, wanted_colours_per_band));
% Second band
second_band_colours = zeros(wanted_colours_per_band, 3);
interpolation_interval = size(second_band_template,1) - 1;
second_start = 1 + 0.0*interpolation_interval; % start of interpolation
second_end = 1 + interpolation_interval - 0.0*interpolation_interval; % end of interpolation
second_red_interpolator = griddedInterpolant(second_band_template(:,1));
second_green_interpolator = griddedInterpolant(second_band_template(:,2));
second_blue_interpolator = griddedInterpolant(second_band_template(:,3));
second_band_colours(:,1) = second_red_interpolator(linspace(second_start, second_end, wanted_colours_per_band));
second_band_colours(:,2) = second_green_interpolator(linspace(second_start, second_end, wanted_colours_per_band));
second_band_colours(:,3) = second_blue_interpolator(linspace(second_start, second_end, wanted_colours_per_band));
% Third band
third_band_colours = zeros(wanted_colours_per_band, 3);
interpolation_interval = size(third_band_template,1) - 1;
third_start = 1 + 0.0*interpolation_interval; % start of interpolation
third_end = 1 + interpolation_interval - 0.0*interpolation_interval; % end of interpolation
third_red_interpolator = griddedInterpolant(third_band_template(:,1));
third_green_interpolator = griddedInterpolant(third_band_template(:,2));
third_blue_interpolator = griddedInterpolant(third_band_template(:,3));
third_band_colours(:,1) = third_red_interpolator(linspace(third_start, third_end, wanted_colours_per_band));
third_band_colours(:,2) = third_green_interpolator(linspace(third_start, third_end, wanted_colours_per_band));
third_band_colours(:,3) = third_blue_interpolator(linspace(third_start, third_end, wanted_colours_per_band));
%------------------------------------------------------------
% Total colour map
%------------------------------------------------------------
% Add partial colour maps
readable_colourmap = [out_colour;...
                      first_band_colours;...
                      second_band_colours;...
                      third_band_colours];
readable_colourmap = round(readable_colourmap);
% Values for indexing the colour map
out_index = 0;
first_band_index  = 1 + (1/wanted_colours_per_band)*(0:wanted_colours_per_band-1);
second_band_index = 2 + (1/wanted_colours_per_band)*(0:wanted_colours_per_band-1);
third_band_index  = 3 + (1/wanted_colours_per_band)*(0:wanted_colours_per_band-1);
colour_indices = [out_index, first_band_index, second_band_index, third_band_index];
colour_indices = floor(colour_indices*wanted_colours_per_band) - (wanted_colours_per_band - 2);
colour_indices(1) = 1;
%------------------------------------------------------------
% Colour bar, plain
%------------------------------------------------------------
% Heights of bands
out_height = 200; % matrix row samples
first_band_height = 600; % matrix row samples
second_band_height = 600; % matrix row samples
third_band_height = 600; % matrix row samples
% Limits of bands
first_band_start = out_height + 1; % matrix row samples
first_band_end = first_band_start + first_band_height - 1; % matrix row samples
second_band_start = first_band_end + 1; % matrix row samples
second_band_end = second_band_start + second_band_height - 1; % matrix row samples
third_band_start = second_band_end + 1; % matrix row samples
third_band_end = third_band_start + third_band_height - 1; % matrix row samples
% Image for the colour bar
image_colorbar = zeros(out_height + first_band_height + second_band_height + third_band_height, ...
                       wanted_colours_per_band);
image_colorbar(1:out_height, :) = out_index;
image_colorbar(first_band_start:first_band_end, :) = repmat(first_band_index, first_band_height, 1);
image_colorbar(second_band_start:second_band_end, :) = repmat(second_band_index, second_band_height, 1);
image_colorbar(third_band_start:third_band_end, :) = repmat(third_band_index, third_band_height, 1);
% Convert image values to index within the colour map
image_index_colorbar = floor(image_colorbar*wanted_colours_per_band) - (wanted_colours_per_band - 2);
image_index_colorbar(image_index_colorbar < 0) = 1;
% Display typical RGB
if display_single_palette_flag
    figure,set(gcf,'Color','w'),
    imagesc(image_index_colorbar), colorbar
    title('Typical'), set(gca,'FontSize',16),
    colorbar_axis = gca;
    colormap(colorbar_axis, uint8(readable_colourmap))
    colorbar_axis_position = colorbar_axis.Position;
    colorbar_axis.Position = [colorbar_axis_position(1) colorbar_axis_position(2) 0.375 colorbar_axis_position(4)];
end
%------------------------------------------------------------
% Colourmap bar for colour-blindness
%------------------------------------------------------------
if display_single_palette_flag
    visuals = {'Protanopia','Deuteranopia','Tritanopia'};
    % Loop visions
    for visual_index = 1:length(visuals)
        figure,set(gcf,'Color','w'),
        imagesc(image_index_colorbar), colorbar
        title(visuals{visual_index}), set(gca,'FontSize',16),
        colorbar_axis = gca;
        colormap(colorbar_axis, uint8(rgb_to_colour_blindness(readable_colourmap, visuals{visual_index}, blind_conversion_method)))
        colorbar_axis_position = colorbar_axis.Position;
        colorbar_axis.Position = [colorbar_axis_position(1) colorbar_axis_position(2) 0.375 colorbar_axis_position(4)];
    end
end
%------------------------------------------------------------
% Colormap display: typical, protanopia, deuteranopia, tritanopia
%------------------------------------------------------------
if display_multiple_palette_flag
    % Colour-blindness
    visuals = {'None','Protanopia','Deuteranopia','Tritanopia'};
    % Image arrangement
    subplot_size = [2,2];
    % subplot_size = [1, 4];
    % Vertical ticks for band labelling
    ytick_label = out_height + [first_band_height/2 ...
        second_band_height+first_band_height/2 ...
        first_band_height+second_band_height+third_band_height/2];
    figure,set(gcf,'Color','w'),
    % Loop visions
    for visual_index = 1:length(visuals)
        subplot(subplot_size(1), subplot_size(2), visual_index);
        imagesc(image_index_colorbar),
        blindness_flag = ~contains(lower(visuals{visual_index}(1)), 'n');
        if blindness_flag
            title_visual = visuals{visual_index};
        else
            title_visual = 'Typical';
        end   
        title(title_visual), set(gca,'FontSize',16),
        colorbar_axis = gca;
        colormap_aux = uint8(rgb_to_colour_blindness(readable_colourmap, visuals{visual_index}, blind_conversion_method));
        colormap(colorbar_axis, colormap_aux)
        xlabel('Level'),ylabel('Band')
        set(gca,'XTick',1:wanted_colours_per_band,'YTick',ytick_label,'YTickLabel',1:3),
        colorbar_aux = colorbar;
        auxCBarLabel = colorbar_aux.Label;
        set(auxCBarLabel,'String','Colormap index','Rotation',-90,'VerticalAlignment','bottom'); clear auxCBar,
    end
end
%------------------------------------------------------------
% Mixing colours circles: typical, protanopia, deuteranopia, tritanopia
%------------------------------------------------------------
if display_multiple_palette_flag
    % Colour-blindness
    visuals = {'None','Protanopia','Deuteranopia','Tritanopia'};
    % Image arrangement
    subplot_size = [1, 4];
    % Length of square images
    canvas_lenght = 200; % [samples]
    % Colours to represent
    primary_colours = zeros(3, 3);
    primary_colours(1,:) = first_band_template(1,:); % [0...255]
    primary_colours(2,:) = second_band_template(1,:); % [0...255]
    primary_colours(3,:) = third_band_template(1,:); % [0...255]
    % Coordinates of colour circles
    distance_origin = 0.5/(1 + sqrt(2*(1 + sind(30))));
    circle_radius = distance_origin*sqrt(2*(1 + sind(30)));
    first_colour_centre =  [-cosd(30), -sind(30)]*distance_origin; % [x, y]
    second_colour_centre = [ cosd(90),  sind(90)]*distance_origin; % [x, y]
    third_colour_centre =  [ cosd(30), -sind(30)]*distance_origin; % [x, y]
    colour_x_centre = [first_colour_centre(1) second_colour_centre(1) third_colour_centre(1)];
    colour_y_centre = [first_colour_centre(2) second_colour_centre(2) third_colour_centre(2)];
    % Canvas with circles
    [canvas_coordinates_x, canvas_coordinates_y] = meshgrid(linspace(-0.5, 0.5, canvas_lenght));
    canvas_size = numel(canvas_coordinates_x);
    figure,set(gcf,'Color','w'),
    % Loop visuals
    for visual_index = 1:length(visuals)
        % Initialize canvas
        canvas_sheet = zeros(size(canvas_coordinates_x,1), size(canvas_coordinates_y, 2), 3, 'uint8');
        % Colour bands
        new_colours = uint8(rgb_to_colour_blindness(primary_colours, visuals{visual_index}, blind_conversion_method));
        % Loop colour bands
        for colour_index = 1:3
            % Circle
            in_circle = sqrt((canvas_coordinates_x - colour_x_centre(colour_index)).^2 + ...
                             (canvas_coordinates_y - colour_y_centre(colour_index)).^2) < circle_radius;
            in_circle = find(in_circle);
            % Loop RGB index
            for band_index = 1:3
                canvas_sheet(in_circle + (band_index-1)*canvas_size) = ...
                    canvas_sheet(in_circle + (band_index-1)*canvas_size) + ...
                    new_colours(colour_index, band_index);
            end
        end
        % Title
        blindness_flag = ~contains(lower(visuals{visual_index}(1)), 'n');
        if blindness_flag
            title_visual = visuals{visual_index};
        else
            title_visual = 'Typical';
        end
        % Draw circles
        subplot(subplot_size(1), subplot_size(2), visual_index);
        imagesc(canvas_coordinates_x(1,:), canvas_coordinates_y(:,1), canvas_sheet), axis xy, axis square
        title(title_visual), set(gca,'FontSize',16,'XTick',[],'YTick',[]),
    end
end
