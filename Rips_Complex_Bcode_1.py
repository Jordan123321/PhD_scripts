#!/usr/bin/env python3

##This is a flag.
Silent=True
import itertools
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
import skimage.color as skcol


#import csv
import math
#import random as rd

##Imports my Voronoi Split Algorithm into a Module
#import VorSplit
from math import log10, floor

#import seaborn as sns

#import itertools

import sys

import scipy.spatial.distance

Rat=25/(92*1000000)


def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)


##Defines paths to my directories and save file locations
path1=os.getcwd()
path0=os.path.dirname(path1)
path2=path0+ '/Output_Barcode_2'
path3=path1 + '/Barcodes/Sample_Img'
ext='.eps'







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





##Gaussian smooth followed by a triangle filter to determine the boundary of the islets
IsltShp1=filt.gaussian(IsltShp0,10)

##Gaussian smooth followed by a triangle filter to determine the boundary of the islets
IsltShp2=IsltShp1>filt.threshold_triangle(IsltShp1)
##Get rid of any small objects that may yet exist
IsltShp3=skmorph.remove_small_holes(IsltShp2, connectivity=1, area_threshold=1000)
IsltShp4=skmorph.remove_small_objects(IsltShp3, connectivity=1, min_size=10000)


Nuc=Dim[:,:,0]+Dim[:,:,1]+Dim[:,:,2]

##Get rid of scale bar
Nuc[-50:,-200:]=0

##Local threshold
NucPos=(Nuc>filt.threshold_local(Nuc, block_size=101))

##Mask the nuclei with the islet shape
IsltNuc=np.multiply(NucPos,IsltShp4)
##Fill the holes to make contingent nuclei
NucPos2=ndim.binary_fill_holes(IsltNuc)
##Remove small artefacts
NucPos3=skmorph.remove_small_holes(NucPos2, connectivity=1, area_threshold=1000)
NucPos4=skmorph.remove_small_objects(NucPos3, connectivity=1, min_size=100)

##Erode the image, to disentangle joint up DAPI stains
NucPos5=skmorph.erosion(NucPos4, selem=skmorph.disk(7))
##Label and number the nuclei
NucLabs,Cellnum = skmeas.label(NucPos5, return_num=1)

##Properties of regions assigned to an array
props=skmeas.regionprops(NucLabs)
##Make an array to store where the centers are
NPosArr=np.empty((Cellnum,2),dtype=int)
##Loop through and save this
for c, p in enumerate(props):
    NPosArr[c,0],NPosArr[c,1]=round(p.centroid[0]),round(p.centroid[1])


epsilon=70
##Create a distance matrix between the centroids
DistArr=scipy.spatial.distance.pdist(NPosArr)
DistArr1=scipy.spatial.distance.squareform(DistArr)
DistArr2=np.triu(DistArr1)

DistArr3= np.logical_and(DistArr2>0,DistArr2<epsilon)
DistArr4= np.logical_and(DistArr1>0,DistArr1<epsilon)

##Create blank image
canvas=np.full_like(Oim[:,:,0],255,dtype=np.uint8)

##Generate the twoplexes
twoplex=[]
for i in range(Cellnum):
    for j in range (i):
        if DistArr3[j,i] and i!=j:
            for k in range(j):
                if DistArr4[k,i] and DistArr4[k,j] and (k !=j and k!=i):
                    twoplex.append(sorted([i,j,k]))
twoplex.sort()
twoplex=list(twoplex for twoplex,_ in itertools.groupby(twoplex))


##Generate the threeplexes from the twoplexes
threeplex=[]
for i in twoplex:
    for j in range(Cellnum) :
        if (DistArr4[i[0],j] and DistArr4[i[1],j] and DistArr4[i[2],j])and (j not in i):
            threeplex.append(sorted([*i,j]))
    
threeplex.sort()
threeplex=list(threeplex for threeplex,_ in itertools.groupby(threeplex))

#sys.exit()

##Draw the twoplexes
for i in twoplex:
    tri_row=[NPosArr[i[0],0],NPosArr[i[1],0],NPosArr[i[2],0]]
    tri_col=[NPosArr[i[0],1],NPosArr[i[1],1],NPosArr[i[2],1]]
    rr,cc=skdr.polygon(r=tri_row, c=tri_col, shape=canvas.shape)
    canvas[rr,cc]=0



lst0,lst1=np.nonzero(DistArr3)
for i,j in zip(lst0,lst1):
    rr,cc=skdr.line(*NPosArr[i],*NPosArr[j])
    canvas[rr,cc]=0


##Track the B0 Components
label_B0 = skmeas.label(np.logical_not(canvas),connectivity=2)


colmap=plt.cm.hsv
colmap.set_under(color='k')

f, ax = plt.subplots()
ax.set_title(r'B0 Rips Components Coloured')
ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
ax.imshow(label_B0, interpolation='none',cmap=colmap,vmin=0.1)
plt.savefig(path2+'/Rips_Coloured_B0'+ext, interpolation='none')
plt.close()


f, ax = plt.subplots()
ax.set_title(r'Rips Complex')
ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
ax.imshow(canvas, interpolation='none',cmap='Greys')
plt.savefig(path2+'/Rips_Complex_Bcode_1'+ext, interpolation='none')
plt.close()


##Track the B1 components
label_B1 = skseg.clear_border(skmeas.label(canvas,connectivity=1))



f, ax = plt.subplots()
ax.set_title(r'B1 Rips Components Coloured')
ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
ax.imshow(label_B1, interpolation='none',cmap=colmap,vmin=0.1)
plt.savefig(path2+'/Rips_Coloured_B1'+ext, interpolation='none')
plt.close()


