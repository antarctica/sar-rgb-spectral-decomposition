*******************
** Dataset scope **
*******************
Parameter files in Matlab language files ('.m') to process Synthetic Aperture Radar (SAR) products from the British
Antarctic Survey (BAS) airborne radar PASIN2 (Polarimetric Airborne Scientific INstrument, mark 2) for deep-ice sounding.

The SAR products, each matched with a processing-parameter file ('.m'), form a dataset to illustrate the
RGB Doppler Decomposition method on SAR images. This method is a generalised framework to interpret the SAR images,
by which the Doppler or range spectral domains are first split into three sub-bandwidths, later to each of the three
a colour of a triplet of colours is assigned, and finally the three are superposed into one single image by the addition
of the three colours. This is explained in the article (yet to be submitted)

(yet to be submitted)
Alvaro Arenas-Pingarrón, Alex M. Brisbourne, Carlos Martín, Hugh F.J. Corr, Carl Robinson, Tom A. Jordan, Paul V. Brennan: ‘A Generalised Representation Framework of Synthetic Aperture Radar Images for Improved Englacial Observations’, EGU The Cryosphere.

The SAR products generated from these processing-parameter files ('.m') are in netCDF files ('.nc') archived by the UK Polar Data Centre (UK PDC)
-------------- Data DOI to be added --------------

If the decomposition is applied on the Doppler spectrum, the new image contains the directional information
related to the Doppler frequencies: positive frequencies when the radar approaches the target, near-zero
frequencies when the relative distance from radar to target is near stationary, and negative when the radar
leaves it behind. If the backscattering is characterised by a very broad beamwidth the target will be 
gray/white, and if by a very narrow beamwidth then the target will be represented by one of the colours of
the triplet.

In any domain, the decomposition identifies the targets resulting from along-track (azimuth) and/or range
ambiguities. These ambiguous targets will appear as thick rainbows with the colours of the triplet.

With this dataset, the user will be able to apply the RGB Doppler Decomposition method, interpret the
results, modify the different parameters and colours, contrast the results, and conduct new decompositions
according to other datasets and needs.

************************
** PASIN2 description **
************************
PASIN2 (Polarimetric Airborne Scientific INstrument, mark 2) is designed for deep-ice sounding. It is 
a 150MHz pulsed airborne SAR (Synthetic Aperture Radar), with 12 antennas: 4 below port wing, 4 below the
fuselage (belly), and 4 below starboard wing. Wing antennas can transmit and receive, and belly antennas
are receiver-only. In transmission, the 4 antennas below each wing are driven simultaneously, effectively
acting as an array, and hence two array transmitters are considering: one from port, other from starboard.
A SAR image from any combination of transmitter/receivers will inform about the depth profile along the
trajectory of the aircraft, from the shallow and deep internal layers to the glaciar bottom, either
bedrock (ice sheets and ice streams) or ice/water interface (ice shelves). However, the depth profile
can be incorrectly interpreted, as the antenna radiation pattern is wide in the across-track plane
(perpendicular to the aircraft trajectory), and thus the echoes from the englacial or bottom reflections
might imping not only from the vertical direction of arrival (DoA), but from a wide DoA interval. To 
properly determine the depth of the scatterers or bottom, an array of antennas in reception must be used,
which will estimate the across-track DoA. After the DoA, and with the range and along-track dimensions, a
3D space is fully determined, and hence the 3D mapping is possible.

This processing is performed in 3 steps, being the inputs the outputs of the preceeding step:
	1) - SAR processing, from the radar raw data;
	        \___|___/_________________
					  \
				          _\/
	2) - DoA estimation, from several SAR images;
	        \___|___/_____________
				      \
				      _\/
	3) - 3D-mapping, from several DoA estimations.

SAR processing and DoA estimations need of antenna calibration, in amplitude, phase and delay. This is
to correct deviations due to the electronics, multipath reflections on aircraft hull and receiver drifts.
The calibration requires dedicated calibration flights.

Polarimetric, DoA and 3D-mapping products are not covered in this dataset, but they are based on the SAR images.

**************************
** Dataset description **
**************************
The SAR products here included are selected as data sets for testing the RGB Decomposition method in Doppler
domain, which is an alternative for representing the SAR images by decomposing the Doppler spectrum into
three sub-bands that will be overlapped into an RGB image.

The SAR products generated for this datasets are in netCDF files ('.nc') archived by the UK Polar Data Centre (UK PDC)
-------------- Data DOI to be added --------------

The data here included belong to the next two projects:
	1) FISS (Filchner Ice Shelf System), in 2016/17 season. The data set corresponds to a section of the
		flight F11 in Recovery Ice Stream, flying across-flow in the ice Filchner Ice Shelf above a
		sub-glacial water channel near the grounding line. The data set is only range-compressed (RC),
		without SAR focussing, so it is not focussed in the along-track domain, but just in range domain.

	2) BEAMISH (Bed Access, Monitoring and Ice Sheet History) in 2019/20 season. The data set corresponds
		to a section of the flight T04 in Rutford Ice Stream, flying across-flow in the grounded ice 
		from the center of the ice stream towards the ice margin next to the Fletcher Promontory. 
		The data set is a SAR image, thus focussed in the joint range and along-track domains.

The range-compressed (RC) product is here also defined as a SAR product, because althoug it does not include
the SAR processing, it results from a SAR data take.

The summary below aims to explain the scope of each product, but without detailing the processing procedures.

Polarimetric, DoA and 3D-mapping products are not covered here, but they are based on the SAR images.

Each SAR product is matched with a file with the processing parameters written in Matlab language ('.m'),
accessed from GitHub by following the external link
----------------- GITHUB LINK TO BE ADDED -----------------

******************************
** Organisation of products **
******************************
PASIN2 surveys are organised by flights. In season 16/17 there were 31 flights, labelled from F01 to F31. 
F01 to F29 were imaging flights, and F30 and F31 calibration flights. The calibration flights were
above the sea surface and while rolling the aircraft (rotation around the along-track axis).

Here the products are limited to one section of a single flight per season (flight F07 in 2016/17 season, 
and flight T11 in season 2019/20). For a wider list of products see
https://doi.org/10.5285/FAAC4156-047D-47BA-9E31-1A4F766BFDF8

The SAR products generated for this datasets are in netCDF files ('.nc') archived by the UK Polar Data Centre (UK PDC)
-------------- Data DOI to be added --------------

sar_rgb_doppler_decomposition_dataset
|
\--- 2016_2017_Flight11_SAR_ktrace33to35_TxPortC13L4PortJ13L4_RxSC_RangeCompression.nc
     2016_2017_Flight11_SAR_ktrace33to35_TxPortC13L4PortJ13L4_RxSC_RangeCompression_diary.txt
     2019_2020_FlightT04_SAR_ktrace154to163_TxStarC13L4StarJ13L4_RxP1_Backprojection.nc
     2019_2020_FlightT04_SAR_ktrace154to163_TxStarC13L4StarJ13L4_RxP1_Backprojection_diary.txt

This folder contains two kind of SAR products:
1) only with range-compression (RC, only range focussing), with suffix 'RangeCompression';
2) full SAR processing (range and along-track focussing), with suffix 'Backprojection'.

Each product (one RC image from season 2016/17, and one SAR image from season 2019/20) is paired with the diary file
('.txt') trimmed to the content of interest, and a parameter file ('.m', in Matlab language) with the processing
parameters. The '.m' file is accessed from GitHub by following the external link
----------------- GITHUB LINK TO BE ADDED -----------------

The file names contain two numbers related to the first and last trace of the along-track locations. A trace is
the time interval for measuring the along-track time during the flight. In season 16/17, a trace was 200 milisecons;
and in 19/20, 400 miliseconds. The numbers in the file names are kilotraces (ktrace, a thousand of traces), that
is 200 seconds in 16/17, and 400 seconds in 19/20. These numbers are not the same as the block number in the raw
data: each raw data block, as stored, is 20 seconds of flight, thus 100 traces (0.1 kilotraces) in 16/17; and
50 traces (0.05 kilotraces) in 19/20. The raw data blocks start by zero, and in 16/17 and 19/20 each block contains
2500 pulses. For example, if the product file name includes "ktrace33_35", the first kilotrace is 33, which corresponds
to 33000 traces; in 16/17 this is the first pulse of raw data block 330, being the pulse 2500*330+1 within the whole
flight sequence. Similarly, the last kilotrace is 35 (35000 traces), until the last pulse of raw data block 350,
which is the pulse 2500*(350+1). Then, for "ktraceA_B":
 - in 16/17: from first pulse of raw data block 10*A to last pulse of raw data block 10*B
 - in 19/20: from first pulse of raw data block 20*A to last pulse of raw data block 20*B

The SAR product images are in netCDF files, with extension ‘.nc’, and keyword ‘SAR’ in the file name. The complex-domain
images are stored with their real and imaginary parts, respectively in the two-dimensional (joint depth and along-track)
variables “SARimageReal” and “SARimageImag”. The depth is in the variable "radialDepthGrid", and the along-track index
in the variable "profile".

Examples of the content related to SAR products:
	- "2016_2017_Flight11_SAR_ktrace33to35_TxPortC13L4PortJ13L4_RxSC_RangeCompression.nc", with the range-compressed image
		of Flight 07, from kilotrace 73 to 76, transmitting from port side the pulses C13L4 and J13L4, receiving
		from starboard side in antenna SC, without SAR (along-track) focussing.
	- "2016_2017_Flight11_SAR_ktrace33to35_TxPortC13L4PortJ13L4_RxSC_RangeCompression_diary.txt", text file with the
		trimmed record-event information generated during the processing.
	- "2019_2020_FlightT04_SAR_ktrace154to163_TxStarC13L4StarJ13L4_RxP1_Backprojection.nc", with the processed SAR image
		of Flight T04, from kilotrace 154 to 164, transmitting from starboard side the pulses C13L4 and J13L4, receiving
		from port side in antenna P1, processing by back-projection method.
	- "2019_2020_FlightT04_SAR_ktrace154to163_TxStarC13L4StarJ13L4_RxP1_Backprojection_diary.txt", text file with the
		trimmed record-event information generated during the processing.

************************************************
** "SeasonCalibration" folder  (not included) **
************************************************
This folder contains the calibration data for the season, provided the electronics and the antenna locations
do not change throughout the season, for a given antenna configuration (like polarimetric).

For the season 2016/17 see
https://doi.org/10.5285/FAAC4156-047D-47BA-9E31-1A4F766BFDF8

For the season 2019/20 see
https://doi.org/10.5285/FAAC4156-047D-47BA-9E31-1A4F766BFDF8

************************************************
** "FlightCalibration" folder  (not included) **
************************************************
A specific calibration is needed, because the receivers suffered from a drift of delay and advance of 1 sample, 
ocurring at different times for each receiver. This 1-sample drift means a phase error of 90deg, hence needing
to calibrate each receiver independently. This drift is also different for each flight.

For the season 2016/17, flight F07, see
https://doi.org/10.5285/FAAC4156-047D-47BA-9E31-1A4F766BFDF8