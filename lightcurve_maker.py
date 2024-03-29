#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  8 13:32:08 2022

@author: Ben Amblard
"""

import sep
import numpy as np
import numba as nb
from glob import glob
from astropy.io import fits
from astropy.time import Time
from copy import deepcopy
import pathlib
import datetime
import matplotlib.pyplot as plt
import os
from tqdm import trange


def stackImages(folder, save_path, startIndex, numImages, bias, dark, flat, gain, isRCD):
    """make median combined image of first numImages in a directory
    input: current directory of images (path object), directory to save stacked image in (path object), starting index (int), 
    number of images to combine (int), bias image (2d numpy array), gain level ('low' or 'high')
    return: median combined bias subtracted image for star detection"""
    #for .fits files:
    if not isRCD:
        if folder.joinpath('converted.txt').is_file == False:
             with open(folder + 'converted.txt', 'a'):
                os.utime(folder + 'converted.txt')
                if gain == 'high':
                    os.system("python .\\RCDtoFTS.py " + str(folder) + ' ' + gain)
                else:
                    os.system("python .\\RCDtoFTS.py " + str(folder))
        
        fitsimageFileList = sorted(folder.glob('*.fits'))
        fitsimages = []   #list to hold bias data
        
        '''append data from each image to list of images'''
        for i in range(startIndex, numImages):
             fitsimages.append(fits.getdata(fitsimageFileList[i]))
            
        fitsimages = importFramesFITS(folder, fitsimageFileList, startIndex, numImages, bias, dark, flat)[0]
        
        '''take median of images and subtract bias'''
        imageMed = np.median(fitsimages, axis=0)
        imageMed = imageMed - bias
        
    
    else:    #for rcd files:

        '''get list of images to combine'''
        rcdimageFileList = sorted(folder.glob('*.rcd'))         #list of .rcd images
        rcdimages = importFramesRCD(folder, rcdimageFileList, startIndex, numImages, bias, gain)[0]     #import images & subtract bias

        imageMed = np.median(rcdimages, axis=0)          #get median value
        
        
    '''save median combined bias subtracted image as .fits'''
    hdu = fits.PrimaryHDU(imageMed)
    medFilepath = save_path.joinpath(gain + '_medstacked.fits')     #save stacked image
    hdu.writeto(medFilepath, overwrite = True)

    return imageMed

def initialFindFITS(data, detect_thresh):
    """ Locates the stars in the initial time slice 
    input: flux data in 2D array for a fits image, star detection threshold (float)
    returns: [x, y, half light radius] of all stars in pixels"""

    ''' Background extraction for initial time slice'''
    data_new = deepcopy(data)           #make copy of data
    bkg = sep.Background(data_new)      #get background array
    bkg.subfrom(data_new)               #subtract background from data
    thresh = detect_thresh * bkg.globalrms      # set detection threshold to mean + 3 sigma
    print(bkg.globalback)
    print(thresh)

    ''' Identify stars in initial time slice '''
    objects = sep.extract(data_new, thresh)#, deblend_nthresh = 1)


    ''' Characterize light profile of each star '''
    halfLightRad = np.sqrt(objects['npix'] / np.pi) / 2.  # approximate half light radius as half of radius

    
    ''' Generate tuple of (x,y,r) positions for each star'''
    positions = zip(objects['x'], objects['y'], halfLightRad)
    

    return positions


def refineCentroid(data, time, coords, sigma):
    """ Refines the centroid for each star for an image based on previous coords, used for tracking
    input: flux data in 2D array for single fits image, header time of image, 
    coord of stars in previous image, weighting (Gauss sigma)
    returns: new [x, y] positions, header time of image """

    '''initial x, y positions'''
    x_initial = [pos[0] for pos in coords]
    y_initial = [pos[1] for pos in coords]
    
    '''use an iterative 'windowed' method from sep to get new position'''
    new_pos = np.array(sep.winpos(data, x_initial, y_initial, sigma, subpix=5))[0:2, :]
    x = new_pos[:][0].tolist()
    y = new_pos[:][1].tolist()
    
    '''returns tuple x, y (python 3: zip(x, y) -> tuple(zip(x,y))) and time'''
    return tuple(zip(x, y)), time

def timeEvolveFITS(data, t, coords, r, stars, x_length, y_length):
    """ Adjusts aperture based on star drift and calculates flux in aperture 
    input: image data (flux in 2d array), image header times, star coords, 
    x per frame drift rate, y per frame drift rate, aperture length to sum flux in, 
    number of stars, x image length, y image length
    returns: new star coords [x,y], image flux, times as tuple"""

    '''get proper frame times to apply drift'''
    frame_time = Time(t, precision=9, format = 'jd').unix   #current frame time from file header (unix)
    drift_time = frame_time# - coords[1,3]    #time since previous frame [s]

    '''add drift to each star's coordinates based on time since last frame'''
    x = [coords[ind][0] for ind in range(stars)]
    y = [coords[ind][1] for ind in range(stars)]
    
    '''get list of indices near edge of frame'''
    EdgeInds = clipCutStars(x, y, x_length, y_length)
    EdgeInds = list(set(EdgeInds))
    EdgeInds.sort()
    
    '''remove stars near edge of frame'''
    xClip = np.delete(np.array(x), EdgeInds)
    yClip = np.delete(np.array(y), EdgeInds)
    
    '''add up all flux within aperture'''
    photometry = sep.sum_circle(data, xClip, yClip, r, bkgann = (r + 6., r + 11.))
    sepfluxes, sepsigma = (photometry[0]).tolist(), (photometry[1]).tolist()
    
    
    '''set fluxes at edge to 0'''
    for i in EdgeInds:
        sepfluxes.insert(i,0)
        sepsigma.insert(i,0)
        
    '''returns x, y star positions, fluxes at those positions, times'''
    star_data = tuple(zip(x, y, sepfluxes, sepsigma, np.full(len(sepfluxes), frame_time)))
    return star_data


def clipCutStars(x, y, x_length, y_length):
    """ When the aperture is near the edge of the field of view sets flux to zero to prevent 
    fadeout
    input: x coords of stars, y coords of stars, length of image in x-direction, 
    length of image in y-direction
    returns: indices of stars to remove"""

    edgeThresh = 20.          #number of pixels near edge of image to ignore
    
    '''make arrays of x, y coords'''
    xeff = np.array(x)
    yeff = np.array(y) 
    
    '''get list of indices where stars too near to edge'''
    ind = np.where(edgeThresh > xeff)
    ind = np.append(ind, np.where(xeff >= (x_length - edgeThresh)))
    ind = np.append(ind, np.where(edgeThresh > yeff))
    ind = np.append(ind, np.where(yeff >= (y_length - edgeThresh)))
    
    return ind


def fluxCheck(fluxProfile, num):
    """ Checks for problems with the star's light curve (tracking failures etc)
    input: light curve of star (array of fluxes in each image), current star number
    returns: -1 for empty profile,  -2 if too short, -3 for tracking failure, -4 for SNR too low, 0 if good
    if event detected returns frame number and star's light curve"""

    '''' Prunes profiles'''
    #light_curve = np.trim_zeros(fluxProfile)
    light_curve = fluxProfile
    
    if len(light_curve) == 0:
        return -1, light_curve, num
      
    FramesperMin = 2400      #ideal number of frames in a directory (1 minute)
    minSNR = 5             #median/stddev limit

    '''perform checks on data before proceeding'''
    
    #check for too short light curve
    if np.median(light_curve) == 0:
        return -2, light_curve,  num  # reject stars that go out of frame to rapidly
    
    #check for tracking failure
    if abs(np.mean(light_curve[:FramesperMin]) - np.mean(light_curve[-FramesperMin:])) > np.std(light_curve[:FramesperMin]):
        return -3, light_curve, num 
    
    #check for SNR level
    if np.median(light_curve)/np.std(light_curve) < minSNR:
        return -4, light_curve, num # reject stars that are very dim, as SNR is too poor

    #if all good
    return 0, light_curve, num


def getSizeFITS(filenames):
    """ gets dimensions of fits 'video' 
    input: list of filenames in directory
    returns: width, height of fits image, number of images in directory, 
    list of header times for each image"""
    
    '''get names of first and last image in directory'''
    filename_first = filenames[0]
    frames = len(filenames)     #number of images in directory

    '''get width/height of images from first image'''
    file = fits.open(filename_first)
    header = file[0].header
    width = header['NAXIS1']
    height = header['NAXIS2']

    return width, height, frames


def importFramesFITS(parentdir, filenames, start_frame, num_frames, bias, dark, flat):
    """ reads in frames from fits files starting at frame_num
    input: parent directory (minute), list of filenames to read in, starting frame number, how many frames to read in, 
    bias image (2D array of fluxes)
    returns: array of image data arrays, array of header times of these images"""

    imagesData = []    #array to hold image data
    imagesTimes = []   #array to hold image times
    
    '''list of filenames to read between starting and ending points'''
    files_to_read = [filename for i, filename in enumerate(filenames) if i >= start_frame and i < start_frame + num_frames]

    '''get data from each file in list of files to read, subtract bias frame'''
    for filename in files_to_read:
        file = fits.open(filename)
        
        header = file[0].header
        
        ''' Calibration frame correction '''
        data = (file[0].data) - bias #- dark)/flat 
        headerTime = header['JD']
            
        file.close()

        imagesData.append(data)
        imagesTimes.append(headerTime)
         
    '''make into array'''
    imagesData = np.array(imagesData, dtype='float64')
    
    '''reshape, make data type into floats'''
    if imagesData.shape[0] == 1:
        imagesData = imagesData[0]
        imagesData = imagesData.astype('float64')
        
    return imagesData, imagesTimes

######################################
# RCD reading section - MJM 20210827 #
######################################

# Function for reading specified number of bytes
def readxbytes(fid, numbytes):
    for i in range(1):
        data = fid.read(numbytes)
        if not data:
            break
    return data

# Function to read 12-bit data with Numba to speed things up
@nb.njit(nb.uint16[::1](nb.uint8[::1]),fastmath=True,parallel=True)
def nb_read_data(data_chunk):
    """data_chunk is a contigous 1D array of uint8 data)
    eg.data_chunk = np.frombuffer(data_chunk, dtype=np.uint8)"""
    #ensure that the data_chunk has the right length

    assert np.mod(data_chunk.shape[0],3)==0

    out=np.empty(data_chunk.shape[0]//3*2,dtype=np.uint16)
    image1 = np.empty((2048,2048),dtype=np.uint16)
    image2 = np.empty((2048,2048),dtype=np.uint16)

    for i in nb.prange(data_chunk.shape[0]//3):
        fst_uint8=np.uint16(data_chunk[i*3])
        mid_uint8=np.uint16(data_chunk[i*3+1])
        lst_uint8=np.uint16(data_chunk[i*3+2])

        out[i*2] =   (fst_uint8 << 4) + (mid_uint8 >> 4)
        out[i*2+1] = ((mid_uint8 % 16) << 8) + lst_uint8

    return out

def getSizeRCD(filenames):
    """ MJM - Get the size of the images and number of frames """
    filename_first = filenames[0]
    frames = len(filenames)

    width = 2048
    height = 2048

    return width, height, frames

# Function to split high and low gain images
def split_images(data,pix_h,pix_v,gain):
    interimg = np.reshape(data, [2*pix_v,pix_h])

    if gain == 'low':
        image = interimg[::2]
    else:
        image = interimg[1::2]

    return image
 
# Function to read RCD file data
def readRCD(filename):

    hdict = {}

    with open(filename, 'rb') as fid:

        # Go to start of file
        fid.seek(0,0)

        # Serial number of camera
        fid.seek(63,0)
        hdict['serialnum'] = readxbytes(fid, 9)

        # Timestamp
        fid.seek(152,0)
        hdict['timestamp'] = readxbytes(fid, 29).decode('utf-8')

        # Load data portion of file
        fid.seek(384,0)

        table = np.fromfile(fid, dtype=np.uint8, count=12582912)

    return table, hdict

def importFramesRCD(parentdir, filenames, start_frame, num_frames, bias, gain = 'high'):
    """ reads in frames from .rcd files starting at frame_num
    input: parent directory (minute), list of filenames to read in, starting frame number, how many frames to read in, 
    bias image (2D array of fluxes)
    returns: array of image data arrays, array of header times of these images"""
    
    imagesData = []    #array to hold image data
    imagesTimes = []   #array to hold image times
    
    hnumpix = 2048
    vnumpix = 2048
    
    imgain = gain
    
    '''list of filenames to read between starting and ending points'''
    files_to_read = [filename for i, filename in enumerate(filenames) if i >= start_frame and i < start_frame + num_frames]
    
    for filename in files_to_read:


        data, header = readRCD(filename)
        headerTime = header['timestamp']

        images = nb_read_data(data)
        image = split_images(images, hnumpix, vnumpix, imgain)
        image = np.subtract(image,bias)

        #change time if time is wrong (29 hours)
        hour = str(headerTime).split('T')[1].split(':')[0]
        fileMinute = str(headerTime).split(':')[1]
        dirMinute = str(parentdir).split('_')[1].split('.')[1]
      #  dirMinute = '30'
        
        #check if hour is bad, if so take hour from directory name and change header
        if int(hour) > 23:
            
            #directory name has local hour, header has UTC hour, need to convert (+4)
            #for red: local time is UTC time (don't need +4)
            newLocalHour = int(parentdir.name.split('_')[1].split('.')[0])
        
            if int(fileMinute) < int(dirMinute):
                newUTCHour = newLocalHour + 4 + 1     #add 1 if hour changed over during minute
               # newUTCHour = newLocalHour + 1         #FOR RED
            else:
                newUTCHour = newLocalHour + 4
               # newUTCHour = newLocalHour              #FOR RED
        
            #replace bad hour in timestamp string with correct hour
            newUTCHour = str(newUTCHour)
            newUTCHour = newUTCHour.zfill(2)
        
            replaced = str(headerTime).replace('T' + hour, 'T' + newUTCHour).strip('b').strip(' \' ')
        
            #encode into bytes
            headerTime = replaced


        imagesData.append(image)
        imagesTimes.append(headerTime)

    '''make into array'''
    imagesData = np.array(imagesData, dtype='float64')
    
    '''reshape, make data type into floats'''
    if imagesData.shape[0] == 1:
        imagesData = imagesData[0]
        imagesData = imagesData.astype('float64')
        
    return imagesData, imagesTimes

###################
# End RCD section #
###################

def getBias(filepath, numOfBiases, gain, isRCD):
    print('Bias')
    filepath = filepath.joinpath('Bias')
    """ get median bias image from a set of biases (length =  numOfBiases) from filepath
    input: bias image directory (path object), number of bias images to take median from (int), gain level ('low' or 'high')
    return: median bias image """

    if not isRCD:     #for .fits files
        '''get list of bias images to combine'''
        biasFileList = sorted(filepath.glob('*.fits'))
        biases = []   #list to hold bias data
        
        '''append data from each bias image to list of biases'''
        for i in biasFileList:
             biases.append(fits.getdata(i))
     
        
    else: #for rcd files:
        '''get list of images to combine'''
        biasFileList = sorted(filepath.glob('*.rcd'))
        
        #import images, using array of zeroes as bias
        biases = importFramesRCD(filepath, biasFileList, 0, numOfBiases, np.zeros((2048,2048)), gain)[0]
    
    '''take median of bias images'''
    biasMed = np.median(biases, axis=0)
    
    return biasMed

def getDark(filepath, numOfDarks, isRCD):
    print('Dark')
    filepath = filepath.joinpath('Dark')
    if not isRCD:
        darkFilelist = sorted(filepath.glob('*.fits'))
        darks = []
        for i in range(numOfDarks):
            darks.append(fits.getdata(darkFilelist[i]))

    else:
        darkFilelist = sorted(filepath.glob('*.rcd'))

        #import images, using array of zeroes as bias
        darks = importFramesRCD(filepath, darkFilelist, 0, numOfDarks, np.zeros((2048,2048)))[0]

    '''take median of dark images'''
    darkMed = np.median(darks, axis=0)
    
    return darkMed

def getFlat(filepath, numOfFlats, masterDark, isRCD):
    print('Flat')
    filepath = filepath.joinpath('Flat')
    if not isRCD:
        flatFilelist = sorted(filepath.glob('*.fits'))
        flats = []
        for i in range(numOfFlats):
            flats.append(fits.getdata(flatFilelist[i]))

    else:
        flatFilelist = sorted(filepath.glob('*.rcd'))

        #import images, using array of zeroes as flat
        flats = importFramesRCD(filepath, flatFilelist, 0, numOfFlats, np.zeros((2048,2048)))[0]

    '''create master flat'''
    flats = np.array([flat-masterDark for flat in flats])
    flats = np.array([np.array(i)/np.median(i) for i in flats])
    masterFlat = np.median(flats, axis=0)

    return masterFlat



def getDateTime(folder):
    """function to get date and time of folder, then make into python datetime object
    input: filepath 
    returns: datetime object"""
    
    #time is in format ['hour', 'minute', 'second', 'msec']
    #folderDate = str(folder.name).split('_')[0]                 #get date folder was created from its name
    folderDate = '20220607'
    #folderTime = str(folder.name).split('_')[1].split('.')
    folderTime = '19.30.00.000'
    folderTime = folderTime.split('.')
    folderDate = datetime.date(int(folderDate[:4]), int(folderDate[4:6]), int(folderDate[-2:]))  #convert to date object
    folderTime = datetime.time(int(folderTime[0]), int(folderTime[1]), int(folderTime[2]))       #convert to time object
    folderDatetime = datetime.datetime.combine(folderDate, folderTime)                     #combine into datetime object
    
    return folderDatetime

def makeBiasSet(filepath, numOfBiases, savefolder, gain):
    """ get set of median-combined biases for entire night that are sorted and indexed by time,
    these are saved to disk and loaded in when needed
    input: filepath (string) to bias image directories, number of biases images to combine for master
    return: array with bias image times and filepaths to saved biases on disk"""
    
    biasFolderList = [f for f in filepath.iterdir() if f.is_dir()]   #list of bias folders
    ''' create folder for results, save bias images '''
    bias_savepath = savefolder.joinpath(gain + '_masterBiases')

    if not bias_savepath.exists():
        bias_savepath.mkdir()      #make folder to hold master bias images in
        
    #make list of times and corresponding master bias images
    biasList = []
    
    #loop through each folder of biases
    for folder in biasFolderList:
        masterBiasImage = getBias(folder, numOfBiases, gain)      #get median combined image from this folder
        
        #save as .fits file if doesn't already exist
        hdu = fits.PrimaryHDU(masterBiasImage)
        biasFilepath = bias_savepath.joinpath(folder.name + '_' + gain + '_medbias.fits')

        
        if not os.path.exists(biasFilepath):
            hdu.writeto(biasFilepath)
        
        folderDatetime = getDateTime(folder)
        
        biasList.append((folderDatetime, biasFilepath))
    
    #package times and filepaths into array, sort by time
    biasList = np.array(biasList)
    ind = np.argsort(biasList, axis=0)
    #biasList = biasList[ind[:,0]]
    
    return biasList

def chooseBias(obs_folder, MasterBiasList):
    """ choose correct master bias by comparing time to the observation time
    input: filepath to current minute directory, 2D numpy array of [bias datetimes, bias filepaths]
    returns: bias image that is closest in time to observation"""
    
    #current hour of observations
    current_dt = getDateTime(obs_folder)
    
    '''make array of time differences between current and biases'''
    bias_diffs = np.array(abs(MasterBiasList[:,0] - current_dt))
    bias_i = np.argmin(bias_diffs)    #index of best match
    
    '''select best master bias using above index'''
    bias_image = MasterBiasList[bias_i][1]
    #load in new master bias image
    bias = fits.getdata(bias_image)
        
    return bias

 
def getLightcurves(folder, savefolder, ap_r, gain, telescope, detect_thresh, RCDfiles = False):
    """ formerly 'main'
    Detect possible occultation events in selected file and archive results 
    
    input: name of current folder (path object), folder to save results in (path object),
    aperture nadius [px], gain (low or high), telescope name (string), 
    star detection threshold (float)
    
    output: printout of processing tasks, .npy file with star positions (if doesn't exist), 
    .txt file for each occultation event with names of images to be saved, the time 
    of that image, flux of occulted star in image
    """
    
    minutefolder = folder.name
    
    print (datetime.datetime.now(), "Opening:", folder)
        
        
    '''make bias set and load in appropriate master bias image'''
    NumBiasImages = 10              #number of bias images to combine in median bias image

    #get 2d np array with bias datetimes and master bias filepaths
    #MasterBiasList = makeBiasSet(folder.joinpath('Bias'), NumBiasImages, savefolder, gain)
    #bias = chooseBias(folder, MasterBiasList)
    bias = getBias(folder, 10, 'high', RCDfiles)
    dark = getDark(folder, 10, RCDfiles)
    #flat = getFlat(folder, 60, dark, RCDfiles)
    flat = fits.getdata(pathlib.Path('/Volumes/1TB HD/Colibri_Obs/Red, 12 Jun 2022/flat_1.5sigmarejection.fits')) #use if masterflat has already been made
    #bias = np.zeros((2048,2048)) #use if no bias frames have been made

    ''' get list of image names to process'''
    if RCDfiles == True: # Option for RCD or fits import - MJM 20210901
        filenames = sorted(folder.glob('*.rcd'))
    else:
        filenames = sorted(folder.glob('*.fits'))


    for filename in filenames:
        if filename.name.split('_')[0] == '.':
            os.remove(filename)

    if RCDfiles == True: # Option for RCD or fits import - MJM 20210901
        filenames = sorted(folder.glob('*.rcd'))
    else:
        filenames = sorted(folder.glob('*.fits'))
    
    field_name = str(filenames[0].name).split('_')[0]               #which of 11 fields are observed
#    pier_side = str(filenames[0].name).split('-')[1].split('_')[0]  #which side of the pier was scope on
    
    ''' get 2d shape of images, number of image in directory'''
    if RCDfiles == True:
        x_length, y_length, num_images = getSizeRCD(filenames) 
    else:
        x_length, y_length, num_images = getSizeFITS(filenames)

    print (datetime.datetime.now(), "Imported", num_images, "frames")
    print(len(filenames), 'len filenames')

    ''' load/create star positional data'''
    if RCDfiles == True: # Choose to open rcd or fits - MJM
        first_frame = importFramesRCD(folder, filenames, 0, 1, bias, gain)
        headerTimes = [first_frame[1]] #list of image header times
        last_frame = importFramesRCD(folder, filenames, len(filenames)-1, 1, bias, gain)
    else:
        first_frame = importFramesFITS(folder, filenames, 0, 1, bias, dark, flat)      #data and time from 1st image
        headerTimes = [first_frame[1]] #list of image header times
        last_frame = importFramesFITS(folder, filenames, len(filenames)-1, 1, bias, dark, flat) #data and time from last image

    headerTimes = [first_frame[1]]                             #list of image header times
        
    #stack first few images to do star finding
    numtoStack = 10
    startIndex = 2          #don't include 1st image (vignetting)
    print('stacking images %i to %i\n' %(startIndex, numtoStack))
    stacked = stackImages(folder, savefolder, startIndex, startIndex + numtoStack, bias, dark, flat, gain, RCDfiles)
    
    
    #find stars in first image
 #   star_find_results = tuple(initialFindFITS(first_frame[0]))
    star_find_results = tuple(initialFindFITS(stacked, detect_thresh))

        
    #remove stars where centre is too close to edge of frame
    before_trim = len(star_find_results)
    star_find_results = tuple(x for x in star_find_results if x[0] + ap_r + 3 < x_length and x[0] - ap_r - 3 > 0)
    star_find_results = tuple(y for y in star_find_results if y[1] + ap_r + 3 < x_length and y[1] - ap_r - 3 > 0)
    after_trim = len(star_find_results)
    
    print('Number of stars cut because too close to edge: ', before_trim - after_trim)
        
      
    #check number of stars for bad first image
    i = 0  #counter used if several images are poor
    min_stars = 10  #minimum stars in an image
    while len(star_find_results) < min_stars:
        print('too few stars, moving to next image ', len(star_find_results))

        if RCDfiles == True:
            first_frame = importFramesRCD(folder, filenames, 1+i, 1, bias, gain)
            headerTimes = [first_frame[1]]
            star_find_results = tuple(initialFindFITS(first_frame[0], detect_thresh))
        else:
            first_frame = importFramesFITS(folder, filenames, 1+i, 1, bias, dark, flat)
            headerTimes = [first_frame[1]]
            star_find_results = tuple(initialFindFITS(first_frame[0], detect_thresh))

            # star_find_results = tuple(x for x in star_find_results if x[0] > 250)
        
             #remove stars where centre is too close to edge of frame
        star_find_results = tuple(x for x in star_find_results if x[0] + ap_r + 3 < x_length and x[0] - ap_r - 3 > 0)
        star_find_results = tuple(y for y in star_find_results if y[1] + ap_r + 3 < x_length and y[1] - ap_r - 3 > 0)
        i += 1
             #check if out of bounds
        if (1+i) >= num_images:
            print('no good images in minute: ', folder)
            print (datetime.datetime.now(), "Closing:", folder)
            print ("\n")
            return -1

    #print('star finding file index: ', i)
        
        #save radii and positions as global variables
    star_find_results = np.array(star_find_results)
    radii = star_find_results[:,-1]
    prev_star_pos = star_find_results[:,:-1]

   #load in initial star positions from last image of previous minute
    initial_positions = prev_star_pos   
	
    #remove stars that have drifted out of frame
    initial_positions = initial_positions[(x_length >= initial_positions[:, 0])]
    initial_positions = initial_positions[(y_length >= initial_positions[:, 1])]

    #save file with updated positions each minute
    posfile = savefolder.joinpath(field_name + '_' + minutefolder + '_' + gain + '_' + str(detect_thresh) + 'sig_pos.npy')
    np.save(posfile, initial_positions)
    
    num_stars = len(initial_positions)      #number of stars in image
    print(datetime.datetime.now(), 'number of stars found: ', num_stars) 
    
    

    ''' flux and time calculations with optional time evolution '''
      
    #image data (2d array with dimensions: # of images x # of stars)
    data = np.empty([num_images, num_stars], dtype=(np.float64, 5))
    
    #get first image data from initial star positions
    data[0] = tuple(zip(initial_positions[:,0], 
                        initial_positions[:,1], 
                        #sum_flux(first_frame[0], initial_positions[:,0], initial_positions[:,1], ap_r),
                        (sep.sum_circle(first_frame[0], initial_positions[:,0], initial_positions[:,1], ap_r, bkgann = (ap_r+2., ap_r + 4.))[0]).tolist(), 
                        (sep.sum_circle(first_frame[0], initial_positions[:,0], initial_positions[:,1], ap_r, bkgann = (ap_r+2., ap_r + 4.))[1]).tolist(),
                        np.ones(np.shape(np.array(initial_positions))[0]) * (Time(first_frame[1], precision=9, format = 'jd').unix)))

    GaussSigma = np.mean(radii * 2. / 2.35)
    first_drift = refineCentroid(*first_frame, initial_positions, GaussSigma)
    drift_pos = []  #array to hold frame positions
    drift_pos.append(first_drift[0])

    x_drifts, y_drifts = [],[]
    for t in trange(1, num_images):
        #import image
        if RCDfiles == True:
            imageFile = importFramesRCD(folder, filenames, t, 1, bias, gain)
            headerTimes.append(imageFile[1])  #add header time to list
        else:
            imageFile = importFramesFITS(folder, filenames, t, 1, bias, dark, flat)
            headerTimes.append(imageFile[1])  #add header time to list
        
        '''drift computation, changed to calculate drift in each frame'''
        current_drift = refineCentroid(*imageFile, drift_pos[-1], GaussSigma)
        
        x_drift , y_drift = np.mean([current_drift[0][i][0]-drift_pos[-1][i][0] for i in range(len(drift_pos))]),np.mean([current_drift[0][i][1]-drift_pos[-1][i][1] for i in range(len(drift_pos))])
        x_drifts.append(x_drift)
        y_drifts.append(y_drift)

        drift_pos.append([(i[0]+x_drift, i[1]+y_drift) for i in current_drift[0]])
        
        """end drift computation"""

        data[t] = timeEvolveFITS(*imageFile, drift_pos[-1], ap_r, num_stars, x_length, y_length)
                    

    # data is an array of shape: [frames, star_num, {0:star x, 1:star y, 2:star flux, 3: unix_time}]  

    
    results = []
    for star in range(0, num_stars):
        results.append(fluxCheck(data[:, star, 2], star)+(x_drifts, y_drifts))
        
    results = np.array(results, dtype = object)


    ''' data archival '''
    
    #make directory to save lightcurves in
    lightcurve_savepath = savefolder.joinpath(gain + '_' + str(detect_thresh) +  'sig_lightcurves')
    if not lightcurve_savepath.exists():
        lightcurve_savepath.mkdir()      #make folder to hold master bias images in
    
    for row in results:  # loop through each detected event
        star_coords = initial_positions[row[2]]     #coords of occulted star
        #if error, what message to save to file
        error_code = row[0]
        
        if error_code == -1:
            error = 'empty profile'
        elif error_code == -2:
            error = 'light curve too short'
        elif error_code == -3:
            error = 'tracking failure'
        elif error_code == -4:
            error = 'SNR too low'
        else:
            error = 'None'
            
        #text file to save results in
        #saved file format: 'star-#_date_time_telescope_xpos-ypos.txt'
        #columns: fits filename and path | header time (seconds) |  star flux
        savefile = lightcurve_savepath.joinpath('star' + str(row[2]) + "_" + str(minutefolder) + '_' + telescope + "_"+  str(int(star_coords[0])) + "-" + str(int(star_coords[1])) + ".txt")
        
        #open file to save results
        with open(savefile, 'w') as filehandle:
            
            #file header
            filehandle.write('#\n#\n#\n#\n')
            filehandle.write('#    First Image File: %s\n' %(filenames[0]))
            filehandle.write('#    Star Coords: %f %f\n' %(star_coords[0], star_coords[1]))
            filehandle.write('#    DATE-OBS (JD): %s\n' %(headerTimes[0]))
            filehandle.write('#    Telescope: %s\n' %(telescope))
            filehandle.write('#    Field: %s\n' %(field_name))
            filehandle.write('#    Error: %s\n' %(error))
            filehandle.write('#\n#\n#\n')
            filehandle.write('#filename     time      flux      x_drift     y_drift\n')
          
            data = row[1]
            x_drifts = [0]+row[3]
            y_drifts = [0]+row[4]
        
            files_to_save = filenames
            star_save_flux = data             #part of light curve to save
      
            #loop through each frame to be saved
            for i in range(0, len(files_to_save)):  
                filehandle.write('%s %f  %f  %f  %f\n' % (files_to_save[i], float(headerTimes[i][0]-headerTimes[0][0])*24*60, float(star_save_flux[i]), float(x_drifts[i]), float(y_drifts[i])))


    print ("\n")




if __name__ == '__main__':
    import pathlib
    path = pathlib.Path('/Volumes/1TB HD/Colibri_Obs/LSR J1835+3259/2022-06-18')
    getLightcurves(path, path.parents[0], 4, 'high', 'Red', 3, False)