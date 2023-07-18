function rgb_256_blind = rgb_to_colour_blindness(rgb_256_matrix, blindness_type, conversion_method)
%
% Calculates the RGB colours as seen by the colour-blindness due to
%   protanopia, deuteranopia, tritanopia and achromatopsia.
%
% Input:
%     rgb_256_matrix (numeric): [N x 3] matrix with the RGB colours to 
%             be converted, with values between 0 and 255.
%             Example: rgb_256_matrix = [10 80 95; ...]
%     blindness_type (string): case-insensitive colour-blindness kind:  
%             'none', 'protanopia', 'deuteranopia', 'tritanopia' and
%             'achromatopsia'.
%             Example: blindness_type = 'protanopia','P','p','Pro',etc.
%             Example: blindness_type = 'none','N','n','None',etc.
%     conversion_method (positive integer): number indexing the algorithm.
%             Method 4 is preferred, as it is properly referenced and
%             accesed.
%             Example: conversion_method = 4
%
% Output:
%     rgb_256_blind (numeric): [N x 3] matrix with the RGB colours after 
%             the blindness, with values between 0 and 255.
%             Example: rgb_256_blind = [10 80 95; ...]
%
% Dependencies:
%         Called by: universally_readable_colourmap.m
%         Calls: -
%
% Date: 16.02.2023
%
% Author: Alvaro Arenas-Pingarron, British Antarctic Survey
%
% Based on:
%   General:
%       https://daltonlens.org/understanding-cvd-simulation/
%       http://www.vischeck.com/
%       https://daltonlens.org/opensource-cvd-simulation/
%       https://daltonlens.org/colorblindness-simulator
%       https://www.color-blindness.com/coblis-color-blindness-simulator/
%       https://github.com/MaPePeR/jsColorblindSimulator
%
%   Visual references:
%       DaltonLens:
%           https://daltonlens.org/opensource-cvd-simulation/
%           A great number of displays according to the several algorithms.
%       ConvertingColors:
%           https://convertingcolors.com/rgb-color-255_0_0.html?search=RGB(255,0,0)
%           Uses Coblis Version2 (see https://daltonlens.org/opensource-cvd-simulation/)
%
%   Conversion algorithm 1:
%       Uses LMSD65 to convert from XYZ to LMS (normalized to D65), but
%       there are other options such as CIECAM97 and CIECAM02
%       http://mkweb.bcgsc.ca/colorblind/math.mhtml#projecthome
%       http://mkweb.bcgsc.ca/colorblind/code.mhtml
%       
%   Conversion algorithm 2: Machado 2009
%       https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/CVD_Simulation.html
%       https://www.inf.ufrgs.br/~oliveira/pubs_files/CVD_Simulation/Machado_Oliveira_Fernandes_CVD_Vis2009_final.pdf
%       G. M. Machado, M. M. Oliveira and L. A. F. Fernandes, "A
%           Physiologically-based Model for Simulation of Color Vision
%           Deficiency," in IEEE Transactions on Visualization and Computer
%           Graphics, vol. 15, no. 6, pp. 1291-1298, Nov.-Dec. 2009,
%           doi: 10.1109/TVCG.2009.113.
%
%   Conversion algorithm 3: Brettel 1997. The Vischeck display filter for GIMP.
%       Brettel, H., Vienot, F., & Mollon, J. D. (1997). Computerized
%           simulation of color appearance for dichromats. Journal of the
%           Optical Society of America. A, Optics, Image Science, and
%           Vision, 14(10), 2647–2655.
%           https://doi.org/10.1364/josaa.14.002647
%       https://github.com/GNOME/gimp/blob/master/modules/display-filter-color-blind.c
%       https://daltonlens.org/understanding-cvd-simulation/#Comparing-Vi%C3%A9not-1999-and-Brettel-1997
%       http://www.vischeck.com/
%
%   Conversion algorithm 4: Based on Vie'not 1999, but roughly
%       Viénot, F., Brettel, H. and Mollon, J.D. (1999), Digital video
%           colourmaps for checking the legibility of displays by
%           dichromats. Color Res. Appl., 24: 243-252.
%           https://doi.org/10.1002/(SICI)1520-6378(199908)24:4<243::AID-COL5>3.0.CO;2-3
%       https://daltonlens.org/understanding-cvd-simulation/#The-simplification-of-Vi%C3%A9not-1999-for-protanopes-and-deuteranopes
%       https://github.com/joergdietrich/daltonize/tree/main/doc
%

% Matrix to convert from typical linear-rgb to blind linear-rgb
%   output_columns = conversion_matrix * input_columns
%   output_rows = input_rows * transpose(conversion_matrix)

% If 'none'
if not(isempty(strfind(lower(blindness_type(1)), 'n')))
    rgb_256_blind = rgb_256_matrix;
    return
end

directly_in_rgb255_flag = false;
% Conversion from RGB to linear RGB, by Gamma expansion
% https://daltonlens.org/understanding-cvd-simulation/#From-sRGB-to-linearRGB
gamma_factor = 2.4;
ratio_factor = 12.92;
conv_limit = 0.04045;
% Conversion from linear RGB to LMS
if false
    % https://github.com/GNOME/gimp/blob/master/modules/display-filter-color-blind.c
    % https://github.com/DaltonLens/DaltonLens-Python/blob/master/daltonlens/convert.py
    %   class LMSModel_Vischeck_GIMP, parameter vischeck_LMS_from_linearRGB
    rgb2lms = [0.05059983, 0.08585369, 0.00952420; ...
               0.01893033, 0.08925308, 0.01370054; ...
               0.00292202, 0.00975732, 0.07145979];
else
    % https://daltonlens.org/understanding-cvd-simulation/#From-linearRGB-to-LMS
    % https://github.com/DaltonLens/DaltonLens-Python/blob/master/daltonlens/convert.py
    %   class LMSModel_sRGB_SmithPokorny75
    rgb2lms = [17.88240413, 43.51609057,  4.11934969; ...
               3.45564232, 27.15538246,  3.86713084; ...
               0.02995656,  0.18430896,  1.46708614];
    rgb2lms = rgb2lms/100;
end
% Conversion from LMS to linear RGB
lms2rgb = inv(rgb2lms);

% According to chosen method
switch conversion_method
    case 1
        %--------------------------------------------------------
        % not preferred
        %--------------------------------------------------------
        protanopia = [0.170556992 0.829443014 0;...
                      0.170556991 0.829443008 0;...
                     -0.004517144 0.004517144 1];
        deuteranopia = [0.33066007 0.66933993 0;...
                        0.33066007 0.66933993 0;...
                       -0.02785538 0.02785538 1];
        tritanopia = [1 0.1273989 -0.1273989;...
                      0 0.8739093  0.1260907;...
                      0 0.8739093  0.1260907];
        achromatopsia = [0.2126 0.7152 0.0722;...
                         0.2126 0.7152 0.0722;...
                         0.2126 0.7152 0.0722];
        % Chosen matrix
        switch lower(blindness_type(1))
            case 'p'
                blindness_matrix = protanopia;
            case 'd'
                blindness_matrix = deuteranopia;
            case 't'
                blindness_matrix = tritanopia;
            case 'a'
                blindness_matrix = achromatopsia;
        end
        % Conversion
        if directly_in_rgb255_flag
            % Blind conversion directly in rgb(0-255)
            rgb_256_blind = rgb_256_matrix * (blindness_matrix.');
        else
            % Step 1: From sRGB(0-255) to linear-rgb(0-1)
            % Gamma expansion
            linear_rgb = rgb_256_to_linear_01(rgb_256_matrix, gamma_factor, ratio_factor, conv_limit);
            % Step 2: From typical linear-rgb(0-1) to blind linear-rgb(0-1)
            linear_rgb_blind = linear_rgb * (blindness_matrix.');
            % Step 3: From linear-rgb(0-1) to sRGB(0-255)
            % Gamma compression
            rgb_256_blind = linear_01_to_rgb_256(linear_rgb_blind, gamma_factor, ratio_factor, conv_limit);
        end
    case 2
        %--------------------------------------------------------
        % Machado 2009
        %--------------------------------------------------------
        % References:
        % G. M. Machado, M. M. Oliveira and L. A. F. Fernandes, "A
        %   Physiologically-based Model for Simulation of Color Vision
        %   Deficiency," in IEEE Transactions on Visualization and Computer
        %   Graphics, vol. 15, no. 6, pp. 1291-1298, Nov.-Dec. 2009,
        %   doi: 10.1109/TVCG.2009.113.
        %
        protanopia = [0.152286	 1.052583	-0.204868;...
                      0.114503	 0.786281	 0.099216;...
                     -0.003882	-0.048116	 1.051998];
        deuteranopia = [0.367322	0.860646   -0.227968;...
                        0.280085	0.672501    0.047413;...
                       -0.011820	0.042940	0.968881];
        tritanopia = [1.255528	-0.076749	-0.178779;...
                     -0.078411	 0.930809	 0.147602;...
                      0.004733	 0.691367	 0.303900];
        achromatopsia = [0.2126 0.7152 0.0722;...
                         0.2126 0.7152 0.0722;...
                         0.2126 0.7152 0.0722];
        % Chosen matrix
        switch lower(blindness_type(1))
            case 'p'
                blindness_matrix = protanopia;
            case 'd'
                blindness_matrix = deuteranopia;
            case 't'
                blindness_matrix = tritanopia;
            case 'a'
                blindness_matrix = achromatopsia;
        end
        % Conversion
        if directly_in_rgb255_flag
            % Blind conversion directly in rgb(0-255)
            rgb_256_blind = rgb_256_matrix * (blindness_matrix.');
        else
            % Step 1: From sRGB(0-255) to linear-rgb(0-1)
            % Gamma expansion
            linear_rgb = rgb_256_to_linear_01(rgb_256_matrix, gamma_factor, ratio_factor, conv_limit);
            % Step 2: From typical linear-rgb(0-1) to blind linear-rgb(0-1)
            linear_rgb_blind = linear_rgb * (blindness_matrix.');
            % Step 3: From linear-rgb(0-1) to sRGB(0-255)
            % Gamma compression
            rgb_256_blind = linear_01_to_rgb_256(linear_rgb_blind, gamma_factor, ratio_factor, conv_limit);
        end
    case 3
        %--------------------------------------------------------
        % Based on Brettel 1997
        %--------------------------------------------------------
        % References:
        %   Brettel, H., Vienot, F., & Mollon, J. D. (1997). Computerized
        %       simulation of color appearance for dichromats. Journal of the
        %       Optical Society of America. A, Optics, Image Science, and
        %       Vision, 14(10), 2647–2655.
        %       https://doi.org/10.1364/josaa.14.002647
        %   https://github.com/GNOME/gimp/blob/master/modules/display-filter-color-blind.c
        %   https://daltonlens.org/understanding-cvd-simulation/#Comparing-Vi%C3%A9not-1999-and-Brettel-1997
        %
        % Convert from sRGB(0-255) to linear RGB(0-1), with gamma expansion
        if directly_in_rgb255_flag
            % Values of RGB constrained to 0 <= RGB <= 1
            linear_rgb = rgb_256_matrix/255; % 0...1
        else
            % Step 1: From RGB(0-255) to linear-rgb(0-1)
            % Gamma-expansion
            linear_rgb = rgb_256_to_linear_01(rgb_256_matrix, gamma_factor, ratio_factor, conv_limit);
        end
        % Convert from linear RGB(0-1) to LMS
        lms_matrix = linear_rgb * (rgb2lms.');
        red = lms_matrix(:,1);
        green = lms_matrix(:,2);
        blue = lms_matrix(:,3);
        % % * Load the LMS anchor-point values for lambda = 475 & 485 nm (for
        % % * protans & deutans) and the LMS values for lambda = 575 & 660 nm
        % % * (for tritans)
        % % */
        % anchor_0 = 0.08008;  anchor_1  = 0.1579;    anchor_2  = 0.5897;
        % anchor_3 = 0.1284;   anchor_4  = 0.2237;    anchor_5  = 0.3636;
        % anchor_6 = 0.9856;   anchor_7  = 0.7325;    anchor_8  = 0.001079;
        % anchor_9 = 0.0914;   anchor_10 = 0.007009;  anchor_11 = 0.0;
        % /* We also need LMS for RGB=(1,1,1)- the equal-energy point (one of
        % * our anchors) (we can just peel this out of the rgb2lms transform
        % * matrix)
        % */
        anchor_e_0 = sum(rgb2lms(1,:));
        anchor_e_1 = sum(rgb2lms(2,:));
        anchor_e_2 = sum(rgb2lms(3,:));
        switch lower(blindness_type(1))
            case 'p'
                factor_rgb_1 = [0, 2.27376148, -5.92721645];
                factor_rgb_2 = [0, 2.18595384, -4.10029338];
                inflection = anchor_e_2 / anchor_e_1;
                tmp = blue ./ green;
                mod_flag = tmp < inflection;
                lms_matrix(:, 1) = ...
                    factor_rgb_1(1) * red + factor_rgb_1(2) * green + factor_rgb_1(3) * blue;
                lms_matrix(mod_flag, 1) = ...
                    factor_rgb_2(1) * red + factor_rgb_2(2) * green + factor_rgb_2(3) * blue;
            case 'd'
                factor_rgb_1 = [0.43979987, 0, 2.60678902];
                factor_rgb_2 = [0.45746620, 0, 1.87574564];
                inflection = anchor_e_2 / anchor_e_0;
                tmp = blue ./ red;
                mod_flag = tmp < inflection;
                lms_matrix(:, 2) = ...
                    factor_rgb_1(1) * red + factor_rgb_1(2) * green + factor_rgb_1(3) * blue;
                lms_matrix(mod_flag, 2) = ...
                    factor_rgb_2(1) * red + factor_rgb_2(2) * green + factor_rgb_2(3) * blue;
            case 't'
                factor_rgb_1 = [-0.05574292, 0.15892917, 0];
                factor_rgb_2 = [-0.00254865, 0.05313210, 0];
                inflection = anchor_e_1 / anchor_e_0;
                tmp = green / red;
                mod_flag = tmp < inflection;
                lms_matrix(:, 3) = ...
                    factor_rgb_1(1) * red + factor_rgb_1(2) * green + factor_rgb_1(3) * blue;
                lms_matrix(mod_flag, 3) = ...
                    factor_rgb_2(1) * red + factor_rgb_2(2) * green + factor_rgb_2(3) * blue;
        end % End of switch lower(blindness_type(1))
        % Convert from LMS to linear RGB with 0<= RGB <=1
        rgb_01_blind = lms_matrix * (lms2rgb.');
        % Convert Convert from linear RGB(0-1) to sRGB(0-255), with gamma compression
        if directly_in_rgb255_flag
            % Values of RGB constrained to 0 <= RGB <= 255
            rgb_256_blind = 255*rgb_01_blind; % 0...255
        else
            % From linear-rgb(0-1) to RGB(0-255)
            % Gamma compression
            rgb_256_blind = linear_01_to_rgb_256(rgb_01_blind, gamma_factor, ratio_factor, conv_limit);
        end
    case 4
        %--------------------------------------------------------
        % Based on Vie'not 1999
        %--------------------------------------------------------
        % References:
        %       Viénot, F., Brettel, H. and Mollon, J.D. (1999), Digital video
        %           colourmaps for checking the legibility of displays by
        %           dichromats. Color Res. Appl., 24: 243-252.
        %           https://doi.org/10.1002/(SICI)1520-6378(199908)24:4<243::AID-COL5>3.0.CO;2-3
        %   https://daltonlens.org/understanding-cvd-simulation/#The-simplification-of-Vi%C3%A9not-1999-for-protanopes-and-deuteranopes
        %   https://github.com/joergdietrich/daltonize/tree/main/doc

        % Convert from sRGB(0-255) to linear RGB(0-1), with gamma expansion
        if directly_in_rgb255_flag
            % Values of RGB constrained to 0 <= RGB <= 1
            linear_rgb = rgb_256_matrix/255; % 0...1
        else
            % Step 1: From RGB(0-255) to linear-rgb(0-1)
            % Gamma-expansion
            linear_rgb = rgb_256_to_linear_01(rgb_256_matrix, gamma_factor, ratio_factor, conv_limit);
        end
        % Convert from linear RGB(0-1) to LMS
        lms_matrix = linear_rgb * (rgb2lms.');
        % Transform according to blindness type
        switch lower(blindness_type(1))
            case 'p'
                blindness_matrix = [0, 2.02344377, -2.52580405;...
                                    0, 1         ,  0;...
                                    0, 0         ,  1];
            case 'd'
                blindness_matrix = [1,          0,  0;...
                                    0.49420696, 0,  1.24826995;...
                                    0,          0,  1];
            case 't'
                blindness_matrix = [ 1,           0,           0;...
                                     0,           1,           0;...
                                    -0.01224491,  0.07203435,  0];
        end
        lms_matrix = lms_matrix * (blindness_matrix.');
        % Convert from LMS to linear RGB with 0<= RGB <=1
        rgb_01_blind = lms_matrix * (lms2rgb.');
        % Convert Convert from linear RGB(0-1) to sRGB(0-255), with gamma compression
        if directly_in_rgb255_flag
            % Values of RGB constrained to 0 <= RGB <= 255
            rgb_256_blind = 255*rgb_01_blind; % 0...255
        else
            % From linear-rgb(0-1) to RGB(0-255)
            % Gamma compression
            rgb_256_blind = linear_01_to_rgb_256(rgb_01_blind, gamma_factor, ratio_factor, conv_limit);
        end        
end % End of switch conversion_method
rgb_256_blind = max(0, rgb_256_blind);
rgb_256_blind = min(rgb_256_blind, 255);
end

function linear_rgb_01 = rgb_256_to_linear_01(rgb_256_matrix, gamma_factor, ratio_factor, conv_limit)
    % Converts from rgb(0-255) to linear rgb(0-1) with gamma factor, and
    %   applies there the blind conversion. Later goes back to rgb(0-255)
    % Step 1: From RGB(0-255) to linear-rgb(0-1)
    linear_rgb_01 = rgb_256_matrix/(255 * ratio_factor); % [0-1]
    linear_rgb_01(rgb_256_matrix > 255*conv_limit) = ...
        ((rgb_256_matrix(rgb_256_matrix > 255*conv_limit)/255 + 0.055)/1.055).^gamma_factor; % [0-1]
end

function rgb_256 = linear_01_to_rgb_256(linear_rgb_01, gamma_factor, ratio_factor, conv_limit)
    % From linear-rgb(0-1) to RGB(0-255) 
    rgb_256 = (255 * ratio_factor)*linear_rgb_01; % [0-255]
    rgb_256(linear_rgb_01 > (conv_limit/ratio_factor)) = ...
        255 * (1.055*linear_rgb_01(linear_rgb_01 > (conv_limit/ratio_factor)).^(1/gamma_factor) - 0.055); % [0-1]
    rgb_256(rgb_256 < 0) = 0;
end