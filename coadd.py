import numpy as np
from astropy.io import fits
from lightcurve_makerBen import importFramesFITS, importFramesRCD, getBias
import os
import pathlib
from tqdm import trange

def stackImages(folder, save_path, startIndex, numImages, bias, gain, RCDfiles = True):
    """make mean combined image of numImages in a directory
    input: current directory of images (path object), directory to save stacked image in (path object), starting index (int), 
    number of images to combine (int), bias image (2d numpy array), gain level ('low' or 'high')
    return: median combined bias subtracted image for star detection"""
    
    #for .fits files:
    if RCDfiles == False: 
        fitsimageFileList = sorted(folder.glob('*.fits'))
        fitsimageFileList.sort(key=lambda f: int(f.name.split('_')[2].split('.')[0]))
        fitsimages = []   #list to hold bias data
    
        '''append data from each image to list of images'''
    #    for i in range(startIndex, numImages):
    #        fitsimages.append(fits.getdata(fitsimageFileList[i]))
        
        fitsimages, time = importFramesFITS(fitsimageFileList, startIndex, numImages, bias)
        time = time[0]

    
        '''take median of images and subtract bias'''
        imageStacked = np.mean(fitsimages, axis=0)
        hdu = fits.PrimaryHDU(imageStacked)
    
    else:
        #for rcd files:
        '''get list of images to combine'''
        rcdimageFileList = sorted(folder.glob('*.rcd'))         #list of .rcd images
        imageStacked = np.zeros((2048,2048))
        if len(rcdimageFileList)<2400:
            return None
        for i in trange(numImages-startIndex):
                rcdimage, time = importFramesRCD(folder, rcdimageFileList, startIndex+i, 1, bias, gain)    #import images & subtract bias
                time = time[0].split('T')[1]
                time = time.split(':')
                time = time[0]+'_'+time[1]+'_'+time[2]
                
                imageStacked += rcdimage
                  #get average value
        imageStacked /= numImages-startIndex
        '''save median combined bias subtracted image as .fits'''
        hdu = fits.PrimaryHDU(imageStacked) 

    medFilepath = save_path.joinpath(gain +  str(time) +'.fits')     #save stacked image

    #if image doesn't already exist, save to path
    if not os.path.exists(medFilepath):
        hdu.writeto(medFilepath)
   
    return imageStacked


if __name__ == '__main__':
    
    folder = pathlib.Path('D:\\ColibriData\\20220731')

    save_path = pathlib.Path('D:\\BenOrRoman\\Test')
    bias_filepath = pathlib.Path('D:\\ColibriData\\20220731\\Bias\\20220731_04.44.41.609')

    gain = 'high'

    bias = getBias(filepath = bias_filepath, numOfBiases = 10, gain = gain, isRCD = True)
   

    for minute in os.listdir(folder):
        print(minute, '\n')
        stackImages(folder.joinpath(minute), save_path, 1, 1200, bias, gain)
        stackImages(folder.joinpath(minute), save_path, 1201, 2400, bias, gain)

