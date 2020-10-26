#!/usr/bin/env python3

##This is a flag.
Silent=True

##These are library imports. OS is for file manipulation
import os
##Numpy is for mathematical manipulation, especially with matrices
import numpy as np

##Matplotlib is for processing images into and out of python, and graphing
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.gridspec import GridSpec
##Imports latex labels
from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
from matplotlib.backends.backend_pdf import PdfPages


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Calibri'
##Scipy contains some elements used in image processing that are required here
#import scipy.spatial as scspat
import scipy.ndimage as ndim
##This imports the statistics module from scipy
#import scipy.stats as scstats

##Scikit image is a package which contains a large number of image processing functions
import skimage.io as skio
import skimage.morphology as skmorph
import skimage.filters as filt
import skimage.measure as skmeas
#import skimage.feature as skfea
import skimage.segmentation as skseg
import skimage.draw as skdr
#import skimage.color as skcol

#import csv
import math
#import random as rd

##Imports my Voronoi Split Algorithm into a Module
#import VorSplit
from math import log10, floor

#import seaborn as sns

#import itertools

import sys


Rat=25/(92*1000000)


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)


##Defines paths to my directories and save file locations
path1=os.getcwd()
path0=os.path.dirname(os.path.dirname(path1))
path2=path0+ '/Output_Barcode'
path3=path1 + '/Sample_Img'
ext='.eps'







res=100
SSpace=np.linspace(start=0,stop=20,num=res)
cval=10



##Set the scan Beta, Delta and Alpha cell count to zero, and make an image number string
SmI1=0
SmS1=0
SmG1=0
Dim=skio.imread(path3+'/DAPI.tif')
Gim=skio.imread(path3+'/GCG.tif')
Iim=skio.imread(path3+'/INS.tif')
Sim=skio.imread(path3+'/SST.tif')
Oim=skio.imread(path3+'/Olay.tif')
#print(Dim.shape)
##Get rid of the backround blood
Sim1=Sim[:,:,0]+Sim[:,:,1]+Sim[:,:,2]




##Take the Laplacian of the Stomatostatin
##Get rid of the scale bar
DiffSim1=filt.laplace(Sim1)
##Get rid of scale bar
DiffSim1[-50:,-200:]=0
B0Temp=np.zeros(res)
B1Temp=np.zeros(res)

##Smooth it
DiffSim2=filt.gaussian(DiffSim1,1)

#Filter it
DiffSim3=DiffSim2>filt.threshold_triangle(DiffSim2)
##Mask the original image with the blood removed
Stcked=np.stack([DiffSim3,DiffSim3,DiffSim3],axis=2)
Sim=np.multiply(Stcked,Sim)

##First thing is to find the islet shape. This can be done by adding the non Dapi, smoothing and thresholding
##Add together the three scans
IsltShp=Gim+Iim+Sim

##Add together the red blue and green to flatten the array
IsltShp0=IsltShp[:,:,0]+IsltShp[:,:,1]+IsltShp[:,:,2]
##Get rid of scale bar
IsltShp0[-50:,-200:]=0


for ssn,ss in enumerate(SSpace):


    ##Gaussian smooth followed by a triangle filter to determine the boundary of the islets
    IsltShp1=filt.gaussian(IsltShp0,10)
    
    ##Gaussian smooth followed by a triangle filter to determine the boundary of the islets
    IsltShp2=IsltShp1>filt.threshold_triangle(IsltShp1)
    ##Get rid of any small objects that may yet exist
    IsltShp3=skmorph.remove_small_holes(IsltShp2, connectivity=1, area_threshold=1000)
    IsltShp4=skmorph.remove_small_objects(IsltShp3, connectivity=1, min_size=10000)
    
    
    IsltLabs, B0 = skmeas.label(IsltShp4, return_num=1)
    IsltHoles, B1 = skmeas.label(skseg.clear_border(np.logical_not(IsltShp4)), return_num=1)
    B0Temp[ssn]=B0
    B1Temp[ssn]=B1
    print(ss)

f, ax = plt.subplots(figsize=(20,5))
plt.tick_params(labelsize=20)
ax.set_xlabel(r'Islet Gauss Blur Parameter Value', fontsize=25)
ax.set_ylabel(r'B0 value', rotation='vertical', fontsize=25)
ax.plot(SSpace, B0Temp, 'k-')
ax.axvline(cval,ls='--',c='r')
ax.set_xlim(SSpace[0],SSpace[-1])
ax.set_ylim(0,1.1*B0Temp.max())
plt.tight_layout()
plt.savefig(path2+'/Individ_B0_Gauss_Blur'+ext, interpolation='none')
plt.close()


f, ax = plt.subplots(figsize=(20,5))
plt.tick_params(labelsize=20)
ax.set_xlabel(r'Islet Gauss Blur Parameter Value', fontsize=25)
ax.set_ylabel(r'B1 value', rotation='vertical', fontsize=25)
ax.plot(SSpace, B1Temp, 'k-')
ax.axvline(cval,ls='--',c='r')
ax.set_xlim(SSpace[0],SSpace[-1])
ax.set_ylim(0,1.1*B1Temp.max())
plt.tight_layout()
plt.savefig(path2+'/Individ_B1_Gauss_Blur'+ext, interpolation='none')
plt.close()