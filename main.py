import lightcurve_maker as make
import lightcurve_looker as look
import relative_photometry as photo
import pathlib
import matplotlib.pyplot as plt
import os
import glob

###Specify data and save paths
data_path = pathlib.Path(input('Path to images: '))
save_path = pathlib.Path(input('Path to save results: '))

###Create the lightcurves from the data
make.getLightcurves(data_path, save_path, 4, 'high', 'Red', 3, False)

lightcurve_directory = save_path.joinpath('/high_3sig_lightcurves')

###Specify the reference stars for relative photometry
stars = [2278, 2373, 2430, 2517, 2321, 2303, 2702, 2183]

#Save lightcurves of all reference stars for later checks


###Create median lightcurve
medLightcurve = photo.get_medianLightcurve(lightcurve_directory, stars)

star_number = input('Index of target star: ')

star_analysed = [ name for name in os.listdir(lightcurve_directory)[0] if star_number in name ]

###Correct target lightcurve
correctedLC = photo.correctLightcurve(lightcurve_directory.joinpath(star_analysed), medLightcurve)

###Save pre and post correction lightcurves
look.lookLightcurve(star_analysed+'_corrected', correctedLC, save_path)



fig, ax1 = plt.subplots()
plt.scatter([i/650*5.7 for i in range(len(medLightcurve))], medLightcurve)

ax1.set_xlabel('time (hours)')
ax1.set_ylabel('Counts/circular aperture (normalized)')
ax1.set_title('Normalised average lightcurve (weighted with brightness)')
plt.savefig(save_path.joinpath('averageLightcurve.png'))
plt.close()