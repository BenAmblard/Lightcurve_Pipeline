from copy import deepcopy
import relative_photometry as photometry
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import lombscargle
import pandas as pd
import os
import pathlib
import linecache
from astropy.io import fits
from astropy.timeseries import LombScargle
import random 
import time
import matplotlib
from scipy.optimize import leastsq

def getLightcurve(star):
    ''' import lightcurve '''

    flux = pd.read_csv(star, delim_whitespace = True, names = ['filename', 'time', 'flux','x_drift', 'y_drift'], comment = '#')
    x_drift = flux['x_drift']
    y_drift = flux['y_drift']
    flux['flux'] = flux['flux']
    coords = linecache.getline(str(star),6).split(': ')[1].split(' ')
    starnum = star.name.split('star')[1].split('_')[0]
    med = np.median(flux['flux'])                           #median flux
    std = np.std(flux['flux'])                              #standard deviation of flux
    SNR = med/std                   

    return {'flux': flux, 'coords': coords, 'median': med, 'std': std, 'SNR': SNR, 'x_drift': x_drift, 'y_drift': y_drift}, starnum

def window(data,prec):
   
    """ Constructs the window function for a given set of x_coordinates

    Parameters
    ----------
    x : array
        x_coordinates of the sample datapoints.
    prec : int
        The length of the window function.

    Returns
    -------
    x_window : array
        .
    y_window : array
        Standard deviation of the period distribution.
   
    """
    x = deepcopy(data)
    #x.sort()
   
    x_w = np.linspace(min(x)-10,max(x)+10,prec)
    x_window = []
    y_window = []
   
    length = len(x)
    i = 0
    j = 0
    while i < prec:
        if(j >= length or x_w[i] < x[j]):
            x_window.append(x_w[i])
            y_window.append(0)
            i += 1
        elif(x_w[i] > x[j]):
            x_window.append(x[j])
            y_window.append(1)
            j += 1
        else:
            x_window.append(x[j])
            y_window.append(1)
            j += 1
            i += 1
    
    return x_window, y_window


def periodogram(lightcurve, savedir, save = True):
    ''' Compute periodogram of lightcurve '''

    starnum = lightcurve[1]
    lightcurve = lightcurve[0]
    med = lightcurve['median']
    std = lightcurve['std']
    flux = lightcurve['flux']
    coords = lightcurve['coords']
    SNR = lightcurve['SNR']
    x_drift = np.array([sum(lightcurve['x_drift'][:i+1]) for i in range(len(lightcurve['x_drift']))])
    y_drift = np.array([sum(lightcurve['y_drift'][:i+1]) for i in range(len(lightcurve['y_drift']))])
    data = flux['flux']

    y_drift, data = clipper(x_drift, data)

    data =  np.array(flux['flux'])/med - 1

    y_drift, data = clipper(y_drift, data)

    yShuffle = deepcopy(y_drift)

    n = 50
    duration = y_drift.ptp()

    freqs = np.linspace(1/duration, n/duration, 100*n)

    #Data shuffling for confidence level
    sigarray = [0 for i in data]
    for i in range(1000):
        #print(i)
        random.shuffle(yShuffle)
        periodShuffle = LombScargle(yShuffle, data).power(freqs, method = 'fast', normalization = 'psd')
        sigarray = [max(sigarray[i], periodShuffle[i]) for i in range(len(sigarray))]


    x_window, y_window = window(y_drift, 10000)
    periodWindow = LombScargle(x_window, y_window).power(freqs, method = 'fast', normalization = 'psd')
    periods = LombScargle(y_drift, data).power(freqs, method = 'fast', normalization = 'psd')
    

    perds = [1/i for i in freqs]
    plt.plot(perds , periods, label = 'Periodogram')
    plt.plot(perds, periodShuffle, label = 'Confidence level ')
    plt.plot(perds, periodWindow, label = 'Window function')
    plt.xlabel('Period (px)')
    plt.ylabel('Power (A.U.)')
    plt.title('Periodogram of star #{}'.format(starnum))
    plt.legend()
    
    plt.yscale('log')
    plt.xscale('log')
    plt.show()
    #plt.savefig('periodogram_star{}.png'.format(starnum))
    #plt.close()

def clipper(x, y):
    ''' clip outlying data '''

    med = np.median(y)
    std = np.std(y)
    threshold = 5
    xClipped, yClipped = [],[]
    for i in range(len(x)):
        if y[i] <= med + threshold*std and y[i] >= med - threshold*std :
            xClipped.append(x[i])
            yClipped.append(y[i])

    return np.array(xClipped), np.array(yClipped)


def sinefit(t,y, guess_f):
    guesses = [3*np.std(y)/np.sqrt(2), guess_f, 0, 0, np.mean(y)]
    print('----- Guesses: -----')
    print('Amplitude ', guesses[0])
    print('Frequency ', guesses[1])
    print('Phase ', guesses[2])
    print('Trend ', guesses[3])
    print('Offset ', guesses[4], '\n')
    to_minimize = lambda x: x[0]*np.sin(x[1]*t+x[2]) + x[3]*t + x[4] - y
    amp, freq, phase, gradient, mean = leastsq(to_minimize, guesses)[0]
    return amp, freq, phase, gradient, mean, guesses

if __name__ == '__main__':   
    directory = pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/lightcurves_2022-06-12')
    lightcurve = directory.joinpath(pathlib.Path('star2441_2022-06-12_Red.txt'))
    stars = [2278, 2373, 2430, 2321, 2303, 2183]
    stars_flat = [2528, 2634, 2694, 2784, 2579, 2556, 2995, 2424]
    medianLightcurve = photometry.get_medianLightcurve(directory, stars)
    lightcurve = photometry.correctLightcurve(lightcurve, medianLightcurve), '2441'
    periodogram(lightcurve, None)

    frequency = 2*np.pi/float(input('Period estimate? ') )#2*np.pi/8
    lightcurve = lightcurve[0]
    flux = lightcurve['flux']
    y_drift = np.array([sum(lightcurve['y_drift'][:i+1]) for i in range(len(lightcurve['y_drift']))])
    data = np.array(flux['flux'])/np.mean(flux['flux'])
    y_drift, data = clipper(y_drift, data)
    amp, freq, phase, gradient, offset, guesses = sinefit(y_drift, data, frequency)
    
    print('----- Best fit: -----')
    print('Amplitude ', amp)
    print('Frequency ', freq)
    print('Phase ', phase)
    print('Trend ', gradient)
    print('Offset ', offset, '\n')
    
    x = np.linspace(min(y_drift), max(y_drift), 1000)
    fit = lambda x: amp*np.sin(freq*x+phase) + gradient*x + offset

    guess = lambda x: guesses[0]*np.sin(guesses[1]*x+guesses[2]) + guesses[3]*x + guesses[4]
    plt.scatter(y_drift, data, marker = '.')
    plt.plot(x, fit(x), color = 'black', label = 'Best fit')
    plt.plot(x, guess(x), 'r--', label = 'Initial guess')
    plt.show()

    new_data = [data[i]/fit(y_drift[i]) for i in range(len(data))]
    plt.scatter(y_drift, new_data, marker = '.')
    plt.show()