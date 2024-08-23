# Automatic photometry
 Semi-automatic photometry of astronomical images.
 
 **Better fit** asks for a directory with images that have relatively same center coordinates (+/- 1 arcminute), a catalog with coordinates and magnitudes of stars (Gaia in this case), and an export directory. In the process it automatically finds bright stars, makes a coordinate net and calculates flux of stars that have chosen as reference and check stars. In the end, it outputs photometry files with a zero-point correction of a reference star that has been found using least squares and deviations of check star magnitudes from catalogue ones.
 
 **Plot makers** make plots of photometry files that has been obtained by better fit, excluding **Plot maker** for a list which makes plots for ATLAS photometry files
