from statistics import median
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pathlib
import linecache
from astropy.io import fits

def get_medianLightcurve(directory, stars):
    
    lightcurves = []    
    coordinates = []
    #get a list of files and sort properly
    files = os.listdir(directory)
    
    files = [x for x in directory.iterdir() if x.is_file()]
    files = [x for x in directory.iterdir() if 'stars_snr.txt' not in x.name]
    
    #loop through each star
    for filename in files:
        if filename.suffix == ".txt" and filename.name.split('_')[0]!='.' and int(filename.name.split('star')[1].split('_')[0]) in stars : 
            #make dataframe containing image name, time, and star flux value for the star
            flux = pd.read_csv(directory.joinpath(filename), delim_whitespace = True, names = ['filename', 'time', 'flux','x_drift', 'y_drift'], comment = '#')

            #star X, Y coordinates
                 #star number
            starnum = int(filename.name.split('star')[1].split('_')[0])
            coords = linecache.getline(str(directory.joinpath(filename)),6).split(': ')[1].split(' ')
            coordinates.append((starnum, coords))

            med = np.median(flux['flux'])  
            flux = flux['flux']                                 #median flux

            #add star 
            lightcurves.append((flux, med))
        else:
            continue

    medLightcurve = np.sum([i[0]*i[1] for i in lightcurves], axis=0)/np.sum([i[1] for i in lightcurves])
    medLightcurve = medLightcurve/np.median(medLightcurve)
    return medLightcurve

def correctLightcurve(star, medianLightcurve):
    flux = pd.read_csv(star, delim_whitespace = True, names = ['filename', 'time', 'flux','x_drift', 'y_drift'], comment = '#')
    flux['flux'] = flux['flux']/medianLightcurve
    coords = linecache.getline(str(star),6).split(': ')[1].split(' ')
    starnum = star.name.split('star')[1].split('_')
    med = np.median(flux['flux'])                           #median flux
    std = np.std(flux['flux'])                              #standard deviation of flux
    SNR = med/std                   

    return {'flux': flux, 'coords': coords, 'median': med, 'std': std, 'SNR': SNR}

def lookLightcurve(star, lightcurve):
    #get star info
    med = lightcurve['median']
    std = lightcurve['std']
    flux = lightcurve['flux']
    coords = lightcurve['coords']
    SNR = lightcurve['SNR']
    
    #fix time to account for minute rollover
    seconds = []    #list of times since first frame
    t0 = flux['time'][0]    #time of first frame
    passed0 = False
    for t in flux['time']:

        if t < t0:          #check if time has gone back to 0
            t = t + 60.
        
        if passed0:         #check if minute has rolled over
            while t - t0 < seconds[-1]:
                t = t + 60.
        
        passed0 = True
        
        seconds.append((t - t0)/60)
        
#make plot

    fig, ax1 = plt.subplots()

    ax1.scatter(seconds, flux['flux'])
    plt.ylim(bottom = 180000, top = 210000)
    ax1.hlines(med, min(seconds), max(seconds), color = 'black', label = 'median: %i' % med)
    ax1.hlines(med + std, min(seconds), max(seconds), linestyle = '--', color = 'black', label = 'stddev: %.3f' % std)
    ax1.hlines(med - std, min(seconds), max(seconds), linestyle = '--', color = 'black')
    ax1.set_xlabel('time (hours)')
    ax1.set_ylabel('Counts/circular aperture')
    ax1.set_title('Star #%s [%.1f, %.1f], SNR = %.2f' %(star, float(coords[0]), float(coords[1]), SNR))
    ax1.legend()
    directory = pathlib.Path(os.getcwd())
    plt.savefig(directory.joinpath('star' + star + '.png'), bbox_inches = 'tight')
    #plt.show()
    plt.close()


if __name__ == '__main__':
    directory = pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/lightcurves_2022-06-12')
    stars = [2278, 2373, 2430, 2517, 2321, 2303, 2702, 2183] #20220612
    #stars = [6688, 6952, 6769, 7123, 7990, 7213, 8039] #20220618
    medLightcurve = get_medianLightcurve(directory, stars)
    correctedLC = correctLightcurve(directory.joinpath(pathlib.Path('star2373_2022-06-12_Red.txt')), medLightcurve) #20220612
    #correctedLC = correctLightcurve(directory.joinpath(pathlib.Path('star7000_2022-06-18_Red.txt')), medLightcurve) #20220618
    lookLightcurve('2373_corrected', correctedLC)

    fig, ax1 = plt.subplots()
    plt.scatter([i/650*5.7 for i in range(len(medLightcurve))], medLightcurve)
    
    ax1.set_xlabel('time (hours)')
    ax1.set_ylabel('Counts/circular aperture (normalized)')
    ax1.set_title('Normalised average lightcurve (weighted with brightness)')
    plt.savefig(pathlib.Path(os.getcwd()).joinpath('averageLightcurve.png'))
    plt.close()

