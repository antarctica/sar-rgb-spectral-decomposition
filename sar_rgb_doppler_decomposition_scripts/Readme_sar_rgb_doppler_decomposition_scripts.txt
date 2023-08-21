Folder with Matlab scripts to calculate the RGB Synthetic Aperture Radar (SAR) Doppler Decompositions from SAR products
in NetCDF format, and represent them with different colourmaps. The SAR images correspond to two SAR products collected 
by the British Antarctic Survey (BAS) deep-ice sounding radar PASIN2 (Polarimetric Airborne Scientific INstrument, mark 2)
in Antarctica.

The SAR products generated for this datasets are in netCDF files ('.nc') archived by the UK Polar Data Centre (UK PDC)
(TBC)
Arenas Pingarron, A., Brisbourne, A., Corr, H., Jordan, T., Robinson, C., Martin, C., Nicholls, K., & Smith, A. (2023). Title to be defined (Version 1.0) [Data set]. NERC EDS UK Polar Data Centre. https://doi.org/10.5285/40c2f86b-1a02-4106-934a-42769682df66

and the processing-parameter files ('.m') in the folder 'sar_rgb_doppler_decomposition_dataset\' of the current repository.

Content:
- 'script_article_2016_2017_Flight11_Recovery_v1.m': The data set corresponds to a section of the
	flight F11 in Recovery Ice Stream, in 2016/17 season, flying across-flow in the ice Filchner Ice Shelf
	above a sub-glacial water channel near the grounding line. The data set is only range-compressed (RC),
	without SAR focussing, so it is not focussed in the along-track domain, but just in range domain.
- 'script_article_2019_2020_FlightT04_Rutford_v1.m': The data set corresponds
	to a section of the flight T04 in Rutford Ice Stream, in 2019/20 season, flying across-flow in the
	grounded ice from the center of the ice stream towards the ice margin next to the Fletcher Promontory. 
	The data set is a SAR image, thus focussed in the joint range and along-track domains.
