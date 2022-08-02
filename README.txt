Lightcurve_pipeline (Benjamin Amblard, adapted from Rachel Brown and Mike Mazur)

A pipeline to create lightcurves from a set of .fit (or .rcd) images.
This version of the code is meant to work on longer exposure images in an effort to detect variability in stars.

--------

main.py is the main file and will use the others to come up with a final target lightcurve from a target star and the exposure folder.
All variables are to be set there.


lightcurve_maker.py is responsible for coming up with the raw lightcurves of all the stars present in the image and save them as .txt files in a specified folder.


lightcurve_looker.py helps visualising the lightcurves with matplotlib


relative_photometry.py is where all the lightcurve correction happens, using the lightcurves of reference stars nearby the target.


systematics.py provides with additional function to deal with systematics that may exist after relative photometry has been performed.
Identifies periodicity in lightcurves. 
Was used to caracterise periodic systematics observed when plotting the lightcurves with respect to the star's pixel position on the detector.
Not used in main.py yet

--------

To do:

Tighten the selection of reference stars in relative_photometry.py (findReferenceStars)
Include detection and correction of systematics in main.py
Include full .rcd images support in main.py
Save corrected and average lightcurves as .txt files in lightcurve_maker.py
