''' Computes lightcurves on an imageset, and performs relative photometry on a target star '''

import lightcurve_maker as maker
import lightcurve_looker as looker
import relative_photometry as photo
import pathlib
import matplotlib.pyplot as plt
import os
import glob


###Specify data and save paths

#data_path = pathlib.Path(input('Path to images: ')) #folder name needs to be YYYY-MM-DD
#save_path = pathlib.Path(input('Path to save results: '))
data_path = pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/2022-06-18')
save_path = pathlib.Path('/Users/benji/Documents/UGA/M1 PHY/StageUWO/Lightcurve_Pipeline/Test')

###Create the lightcurves from the data
gain = 'high'
telescope = 'Red'
aperture_radius = 4 #pixels
detection_threshold = 3 #number of sigmas above background level
RCDfiles = False

#comment the next line if the lightcurves have already been computed
maker.getLightcurves(data_path, save_path, aperture_radius, gain, telescope, detection_threshold, RCDfiles)  #saved file format: 'star#_date_telescope_xpos-ypos.txt'


lightcurve_directory = save_path.joinpath('{}_{}sig_lightcurves'.format(gain, detection_threshold))

###Specify the target star for relative photometry. 
# Needs to be identified by user by finding the coordinates in the image and finding the corresponding lightcurve file
star_number = int(input('Index of target star: '))

refStars = photo.findReferenceStars(star_number, str(lightcurve_directory), findradius = 200) #finds reference stars

refStarfiles = []
for lightcurve in os.listdir(lightcurve_directory): #fetches the corresponding lightcurve files
    for i in refStars:
        if 'star{}'.format(i) in lightcurve:
            refStarfiles.append((i,lightcurve))

#Save lightcurves of all reference stars
reflightcurves = save_path.joinpath('refLightcurves')
if not reflightcurves.exists():
        reflightcurves.mkdir() 

for i in refStarfiles:
    star, file = i
    lightcurve = looker.get_Lightcurve(lightcurve_directory.joinpath(file))
    looker.lookLightcurve('{}'.format(star), lightcurve, reflightcurves, '0')

###Create median lightcurve
medLightcurve = photo.get_medianLightcurve(lightcurve_directory, refStars)

star_analysed = [ name for name in os.listdir(lightcurve_directory) if str(star_number) in name ][0]

###Correct target lightcurve
correctedLC = photo.correctLightcurve(lightcurve_directory.joinpath(star_analysed), medLightcurve)

###Save post correction lightcurve
looker.lookLightcurve(star_number+'_corrected', correctedLC, save_path)

