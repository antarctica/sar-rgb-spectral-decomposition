Folder with the Matlab software related to the colourmap generation and visualisation:

- 'universally_readable_colourmap.m': Matlab function to generate colourmaps for the SAR Decomposition. Attempt to build
	an universally readable colourmap by using three colours whose combinations can be distinguished for colour-blindness
	due to protanopia (lack of sensitivity to red), deuteranopia (green) and tritanopia (blue). The function creates
	three palettes of shades for a triplet of colours, each colour to represent each of the three bands of the Decomposition.
	The aim of the colourmap is that not only each of the three basic colours can be distinguished by colour blind interpreters,
	but also the combinations of the palettes, which will occurr depending on the Decomposition.
- 'rgb_to_colour_blindness.m': Matlab function to calculate the RGB colours as seen by  colour-blind interpreters with protanopia,
	deuteranopia and tritanopia.