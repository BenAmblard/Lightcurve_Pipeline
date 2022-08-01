from statistics import median
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import pathlib
import linecache
from astropy.io import fits
import lightcurve_looker as looker

def get_medianLightcurve(directory, stars):
    ''' Get median lightcurve of reference stars for relative photometry
    input: location of lightcurves (pathlib.Path object), reference stars (int list)
    output: median lightcurve (dict) '''
    
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
    ''' Correct the target lightcurve with a median lightcurve of reference stars
    input: target star number (int), median lightcurve (np.array)
    output: corrected lightcurve (python dict)'''

    flux = pd.read_csv(star, delim_whitespace = True, names = ['filename', 'time', 'flux','x_drift', 'y_drift'], comment = '#')
    x_drift = flux['x_drift']
    y_drift = flux['y_drift']
    #Field correctinon
    flux['flux'] = flux['flux']/medianLightcurve
    coords = linecache.getline(str(star),6).split(': ')[1].split(' ')
    starnum = star.name.split('star')[1].split('_')
    med = np.median(flux['flux'])                           #median flux
    
    std = np.std(flux['flux'])                              #standard deviation of flux
    SNR = med/std                   

    return {'flux': flux, 'coords': coords, 'median': med, 'std': std, 'SNR': SNR, 'x_drift': x_drift, 'y_drift': y_drift}

def lookLightcurve(star, lightcurve, drift = '0'):
    ''' used to plot lightcurve against time or x/y position 
    input: star name (string), star lightcurve (dict), 
    variable for abscissa axis ('0' for time , 'x' for x position, 'y' for y position) 
    output: None'''

    #get star info
    med = lightcurve['median']
    std = lightcurve['std']
    flux = lightcurve['flux']
    coords = lightcurve['coords']
    SNR = lightcurve['SNR']
    x_drift = [sum(flux['x_drift'][:i]) for i in range(len(flux['x_drift']))]
    #x_drift_mod = [i%(7.5) for i in x_drift]

    y_drift = [sum(flux['y_drift'][:i]) for i in range(len(flux['y_drift']))]
    #fix time to account for minute rollover
    seconds = []    #list of times since first frame
    t0 = flux['time'][0]    #time of first frame
    passed0 = False

    ''' Time correction (in .txt files, time is between 0 and 59)'''
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
    if drift == '0': #vs time
        seconds, new_flux = clipper(seconds, np.array(flux['flux']))
        ax1.scatter(seconds, new_flux)#/med -1)
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

    elif drift == 'x': #vs x position
        x_drift, new_flux = clipper(x_drift, flux['flux'])

        ax1.scatter(x_drift, new_flux)

        # ax1.hlines(med, min(seconds), max(seconds), color = 'black', label = 'median: %i' % med)
        # ax1.hlines(med + std, min(seconds), max(seconds), linestyle = '--', color = 'black', label = 'stddev: %.3f' % std)
        # ax1.hlines(med - std, min(seconds), max(seconds), linestyle = '--', color = 'black')
        
        ax1.set_xlabel('x position drift (px)')
        ax1.set_ylabel('Counts/circular aperture')
        ax1.set_title('Star #%s [%.1f, %.1f], SNR = %.2f' %(star, float(coords[0]), float(coords[1]), SNR))
        directory = pathlib.Path(os.getcwd())

        ax1.legend()
        plt.savefig(directory.joinpath('star' + star + '_xdrift.png'), bbox_inches = 'tight')
        #plt.show()
        plt.close()

    elif drift == 'y': #vs y position
        y_drift, new_flux = clipper(y_drift, flux['flux'])
        ax1.scatter(y_drift, new_flux)

        # ax1.hlines(med, min(seconds), max(seconds), color = 'black', label = 'median: %i' % med)
        # ax1.hlines(med + std, min(seconds), max(seconds), linestyle = '--', color = 'black', label = 'stddev: %.3f' % std)
        # ax1.hlines(med - std, min(seconds), max(seconds), linestyle = '--', color = 'black')

        ax1.set_xlabel('y position drift (px)')
        ax1.set_ylabel('Counts/circular aperture')
        ax1.set_title('Star #%s [%.1f, %.1f], SNR = %.2f' %(star, float(coords[0]), float(coords[1]), SNR))
        directory = pathlib.Path(os.getcwd())

        ax1.legend()
        plt.savefig(directory.joinpath('star' + star + '_ydrift.png'), bbox_inches = 'tight')
        #plt.show()
        plt.close()

def clipper(x, y):
    ''' clip outlying data'''

    med = np.median(y)
    std = np.std(y)
    threshold = 5
    xClipped, yClipped = [],[]
    for i in range(len(x)):
        if y[i] <= med + threshold*std and y[i] >= med - threshold*std :
            xClipped.append(x[i])
            yClipped.append(y[i])

    return np.array(xClipped), np.array(yClipped)

def saveLightcurve(savefile, lightcurve):
    ''' export corrected lightcurve to .txt '''

    
    flux = lightcurve['flux']
    x_drift = lightcurve['x_drift']
    y_drift = lightcurve['y_drift']
    coords = lightcurve['coords']

    with open(savefile, 'w') as filehandle:
            
            #file header
            filehandle.write('#\n#\n#\n#\n')
            filehandle.write('#    First Image File: %s\n' %(savefile))
            filehandle.write('#    Star Coords: %f %f\n' %(coords[0], coords[1]))
            filehandle.write('#\n#\n#\n')
            filehandle.write('#filename     time      flux      x_drift     y_drift\n')
        
      
            #loop through each frame to be saved
            for i in range(0, len(flux['flux'])):  
                filehandle.write('%s %f  %f  %f  %f\n' % (savefile.name(), float(flux['time'][i]), float(flux['flux'][i]), float(x_drift[i]), float(y_drift[i])))


def findReferenceStars(star, directory, findradius):
    ''' Automatically finds the reference stars in a radius around the target
    input: star number (int), lightcurve directory (str), radius (int)
    output: list of reference stars (list) '''

    ''' Create a dictionary of all stars in the frame'''
    allstars = {}

    for filename in os.listdir(directory):
        if ('star' in filename.split('_')[0] and filename.split('.')[1] == "txt"):
            starnum = filename.split('_')[0].split('star')[1]
            allstars['{}'.format(starnum)] = filename

    starfile = allstars['%i'%star]

    starcoords = (int(starfile.split('_')[3].split('.')[0].split('-')[0]), int(starfile.split('_')[3].split('.')[0].split('-')[1]))

    refStars = []

    ''' Filter stars based on distance to target star, maximum brightness and SNR '''
    for key, filename in allstars.items():
        coords = (int(filename.split('_')[3].split('.')[0].split('-')[0]), int(filename.split('_')[3].split('.')[0].split('-')[1])) 

        if coords != starcoords and (coords[0]-starcoords[0])**2+(coords[1]-starcoords[1])**2 <= findradius**2:
            starPath = pathlib.Path(directory+'/'+filename)
            flux = pd.read_csv(starPath, delim_whitespace = True, names = ['filename', 'time', 'flux','x_drift', 'y_drift'], comment = '#')['flux']
            std = np.std(flux)
            med = np.median(flux)
            SNR = med/std
            
            if max(flux) < 50000 and SNR>20:
                refStars.append(int(key))


    return refStars


if __name__ == '__main__':
    directory = pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/high_3sig_lightcurves')
    dirname = str(directory)

    star = 6577
    stardir = 'star6577_2022-06-18_Red_934-912.txt'

    refStars = findReferenceStars(star, dirname, 200)
    # #stars = [2278, 2373, 2430, 2321, 2303, 2183, 2441, 2423] #20220612 , removed 2702, 2517
    # stars_flat = [2528, 2634, 2694, 2784, 2579, 2556, 2995, 2424]
    # stars_bias = [1907, 1988, 2010, 2027, 2070, 2116, 2120, 2129, 2198]
    # stars_flatbias = [2518, 2624, 2684, 2774, 2569, 2546, 2985, 2413]
    # stars = [6387, 6682, 6779, 6846, 6876, 6878] #20220618
    # star = 2083
    # #star_bias = 2278
    # #star_flat = 2424

    lightcurve = looker.get_Lightcurve(pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/high_3sig_lightcurves/{}'.format(stardir)))
    looker.lookLightcurve('{}'.format(star), lightcurve, pathlib.Path(os.getcwd()), '0')

    medLightcurve = get_medianLightcurve(directory, refStars)
    correctedLC = correctLightcurve(directory.joinpath(pathlib.Path('{}'.format(stardir))), medLightcurve) #20220612

    lookLightcurve('{}_Bias_corrected'.format(star), correctedLC, '0')
    lookLightcurve('{}_Bias_corrected'.format(star), correctedLC, 'x')
    lookLightcurve('{}_Bias_corrected'.format(star), correctedLC, 'y')

    fig, ax1 = plt.subplots()
    x,y=clipper([i for i in range(len(medLightcurve))],medLightcurve)
    plt.scatter(x,y) #[i/650*5.7 for i in range(len(medLightcurve))], 
    
    ax1.set_xlabel('time (hours)')
    ax1.set_ylabel('Counts/circular aperture (normalized)')
    ax1.set_title('Normalised average lightcurve (weighted with brightness)')
    plt.savefig(pathlib.Path(os.getcwd()).joinpath('averageLightcurve.png'))
    #plt.show()
    plt.close()
