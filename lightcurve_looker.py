import matplotlib.pyplot as plt
import pathlib
import os
import numpy as np
import pandas as pd
import linecache

def get_Lightcurve(file):
    ''' retrieve lightcurve data from .txt files
    input: file name (pathlib.Path object)
    output: lightcurve (dict) '''

    lightcurve = {}

    flux = pd.read_csv(file, delim_whitespace = True, names = ['filename', 'time', 'flux', 'x_drift', 'y_drift'], comment = '#')

    starnum = int(file.name.split('star')[1].split('_')[0])
    coords = linecache.getline(str(file),6).split(': ')[1].split(' ')
    med = np.median(flux['flux'])                                   #median flux
    std = np.std(flux['flux'])
    SNR = med/std

    lightcurve = {'flux': flux, 'coords': coords, 'median': med, 'std': std, 'SNR': SNR}
    return lightcurve

def lookLightcurve(star, lightcurve, save_path, drift = '0'):
    ''' used to plot lightcurve against time or x/y position 
    input: star number (int), star lightcurve (dict), 
    variable for abscissa axis ('0' for time , 'x' for x position, 'y' for y position) 
    output: None'''

    #get star info
    med = lightcurve['median']
    std = lightcurve['std']
    flux = lightcurve['flux']
    coords = lightcurve['coords']
    SNR = lightcurve['SNR']
    x_drift = [sum(flux['x_drift'][:i]) for i in range(len(flux['x_drift']))]

    y_drift = [sum(flux['y_drift'][:i]) for i in range(len(flux['y_drift']))]
    
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
    if drift == '0':
        ax1.scatter(seconds, flux['flux'])
        ax1.hlines(med, min(seconds), max(seconds), color = 'black', label = 'median: %i' % med)
        ax1.hlines(med + std, min(seconds), max(seconds), linestyle = '--', color = 'black', label = 'stddev: %.3f' % std)
        ax1.hlines(med - std, min(seconds), max(seconds), linestyle = '--', color = 'black')
        ax1.set_xlabel('time (hours)')
        ax1.set_ylabel('Counts/circular aperture')
        ax1.set_title('Star #%s [%.1f, %.1f], SNR = %.2f' %(star, float(coords[0]), float(coords[1]), SNR))
        
        ax1.legend()
        directory = save_path
        plt.savefig(directory.joinpath('star' + star + '.png'), bbox_inches = 'tight')
        #plt.show()
        plt.close()

    elif drift == 'x':
        ax1.scatter(x_drift, flux['flux'])
        ax1.set_xlabel('x position drift (px)')
        ax1.set_ylabel('Counts/circular aperture')
        ax1.set_title('Star #%s [%.1f, %.1f], SNR = %.2f' %(star, float(coords[0]), float(coords[1]), SNR))
        ax1.legend()
        directory = save_path
        plt.savefig(directory.joinpath('star' + star + '_xdrift.png'), bbox_inches = 'tight')
        #plt.show()
        plt.close()

    elif drift == 'y':
        ax1.scatter(y_drift, flux['flux'])
        ax1.set_xlabel('y position drift (px)')
        ax1.set_ylabel('Counts/circular aperture')
        ax1.set_title('Star #%s [%.1f, %.1f], SNR = %.2f' %(star, float(coords[0]), float(coords[1]), SNR))
        ax1.legend()
        directory = save_path
        plt.savefig(directory.joinpath('star' + star + '_ydrift.png'), bbox_inches = 'tight')
        #plt.show()
        plt.close()

if __name__ == '__main__':
    star = 6620
    lightcurve = get_Lightcurve(pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/high_3sig_lightcurves/star{}_2022-06-18_Red.txt'.format(star)))
    lookLightcurve('{}'.format(star), lightcurve, pathlib.Path(os.getcwd()), '0')
    #lookLightcurve('2441', lightcurve, pathlib.Path(os.getcwd()), '0')