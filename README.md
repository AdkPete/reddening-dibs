# reddening-dibs
This repository contains codes for measuring reddening from equivalent
width measurements of Diffuse Interstellar Bands (DIBs) and the Na D 
doublet. Used for measuring E(B-V) towards classical novae, although in
principle it can apply to other sources as well.

Basic usage is to run the script measure_ebv.py. This will run through
a set of eight DIBs and the two Na lines for the spectra that are
provided. For each line a plot will be shown, and options are provided
for making equivalent width measurements. The blue line shows the
actual spectrum, the orange line is a continuum model, and the red line
will be integrated under to compute the EW. For each line, enter y
to accept the line measurement. Enter n to skip this line (useful for
non-detections, line blending etc.). Enter m to set the plot into manual
mode, at which point new line edges may be selected by clicking on the
plot wondow.

Once all lines have been selected, this will produce a table showing
the EW and E(B-V) for each line, and the combined results for 
E(B-V).

To run measure_ebv.py, enter three command line arguments. The first
is a boolean determining if DIB lines are to be used. The second is also
a boolean, and determines if Na lines are to be used. The third points
to the desired spectra, either provide a filename (or list of filenames
in the form of additional arguments), or a directory containing the 
spectra to be analyzed.

If you have questions, comments, feature requests or suggestions, please
reach out to Dr. Peter Craig at craigpe1@msu.edu.
