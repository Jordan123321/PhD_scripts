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
#import matplotlib.mlab as mlab
#from matplotlib.gridspec import GridSpec
##Imports latex labels
#from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
#from matplotlib.backends.backend_pdf import PdfPages


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
#import math
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




##Create a distance matrix between the centroids
DistArr=scipy.spatial.distance.pdist(NPosArr)
DistArr1=scipy.spatial.distance.squareform(DistArr)
DistArr2=np.triu(DistArr1)
B0Barcode=[]
B1Barcode=[]

B0num=[]
B1num=[]

epsstp=1
epsmx=141
for epsilon in np.arange(0,epsmx,epsstp):
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
    
    
    
    ##Draw the twoplexes
    for i in twoplex:
        tri_row=[NPosArr[i[0],0],NPosArr[i[1],0],NPosArr[i[2],0]]
        tri_col=[NPosArr[i[0],1],NPosArr[i[1],1],NPosArr[i[2],1]]
        rr,cc=skdr.polygon(r=tri_row, c=tri_col, shape=canvas.shape)
        canvas[rr,cc]=0
    
    
    ##Draw the oneplexes
    lst0,lst1=np.nonzero(DistArr3)
    for i,j in zip(lst0,lst1):
        rr,cc=skdr.line(*NPosArr[i],*NPosArr[j])
        canvas[rr,cc]=0
    
    ##Draw the zeroplexes
    for i in NPosArr:
        canvas[i[0],i[1]]=0
    
    
    
    
    
    ##Track the B0 Components
    label_B0, B0_num = skmeas.label(np.logical_not(canvas),connectivity=2, return_num=True)
    B0props=skmeas.regionprops(label_B0)
    for reg in B0props:
        ##Return the coord list of the simplex
        clist=reg.coords.tolist()
        #returns the 0simplex points that constittue the simplex
        B0smpco=np.array([[*i,num] for num,i in enumerate(NPosArr) if list(i) in clist])
        ##Find the lowest value of zero smplex label
        B0label=np.amin(B0smpco,axis=0)[2]
        B0label2 = ''
        for i in B0smpco:
            B0label2 += str(i[2])+'_'
        B0label2=B0label2[:-1]
        ##append this and epsilon value to B0 barcode
        B0Barcode.append([epsilon,B0label,B0label2])
    B0num.append(B0_num)
    
    
    ##Track the B1 components
    label_B1_ , B1_num = skmeas.label(canvas,connectivity=1,return_num=True)
    label_B1=skseg.clear_border(label_B1_)
    B1props=skmeas.regionprops(label_B1)
    B1_num = len(B1props)
    for reg in B1props:
        ##Dilate region by a pixel to find bounding zero simplexes
        canvas1=np.zeros_like(Oim[:,:,0],dtype=bool)
        canvas1[reg.bbox[0]:reg.bbox[2],reg.bbox[1]:reg.bbox[3]]=reg.image
        canvas1=dilated=skmorph.dilation(canvas1, selem=skmorph.disk(3))
        ##Return the nonzero coordinates of the region
        l1,l2=canvas1.nonzero()
        clist=[[a,b] for a,b in zip(l1,l2)]
        #returns the 0simplex points that constittue the simplex
        try:
            B1smpco=np.array([[*i,num] for num,i in enumerate(NPosArr) if list(i) in clist])
            B1label=np.amin(B1smpco,axis=0)[2]
            B1avg=int(np.mean(B1smpco,axis=0)[2])
            ##append this and epsilon value to B1 barcode
            B1label2 = ''
            for i in B1smpco:
                B1label2 += str(i[2])+'_'
            B1label2=B1label2[:-1]
            B1Barcode.append([epsilon,B1label,B1avg,B1label2])
        except ValueError:
            print("Value Error at " + str(epsilon))
            
    B1num.append(B1_num)
            
    if False:
        print('Plotting here')
        ##plot out the canvases
        colmap=plt.cm.hsv
        colmap.set_under(color='k')
        epsistring='espsilon_'+str(epsilon)
        if not os.path.exists(path2+'/'+epsistring):
                os.makedirs(path2+'/'+epsistring)
        if np.amax(label_B0)>0:
            f, ax = plt.subplots()
            ax.set_title(r'B0 Rips Coloured. Epsilon is' + epsistring)
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(label_B0, interpolation='none',cmap=colmap,vmin=0.1)
            plt.savefig(path2+'/'+epsistring+'/Rips_Coloured_B0'+ext, interpolation='none')
            plt.close()
    
    
        f, ax = plt.subplots()
        ax.set_title(r'Rips Complex. Epsilon is' + epsistring)
        ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
        ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
        ax.imshow(canvas, interpolation='none',cmap='Greys')
        plt.savefig(path2+'/'+epsistring+'/BW_Rips_Complex'+ext, interpolation='none')
        plt.close()
    
        if np.amax(label_B1)>0:
            f, ax = plt.subplots()
            ax.set_title(r'B1 Rips Coloured. Epsilon is' + epsistring )
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(label_B1, interpolation='none',cmap=colmap,vmin=0.1)
            plt.savefig(path2+'/'+epsistring+'/Rips_Coloured_B1'+ext, interpolation='none')
            plt.close()


    print("done for " + str(epsilon))

f, ax = plt.subplots(figsize=(20,5))
plt.tick_params(labelsize=20)
ax.set_title('B0 Plots',fontsize=25)
ax.set_xlabel(r'Epsilon Value',fontsize=25)
ax.set_ylabel('B0 Number',fontsize=25)
ax.plot(B0num, 'k-')
ax.set_xlim([0,epsmx-1])
ax.set_ylim([0,1.1*max(B0num)])
f.tight_layout()
plt.savefig(path2+'/B0_lineplot'+ext, interpolation='none')
plt.close()


f, ax = plt.subplots(figsize=(20,5))
plt.tick_params(labelsize=20)
ax.set_title('B1 Plots',fontsize=25)
ax.set_xlabel(r'Epsilon Value',fontsize=25)
ax.set_ylabel('B1 Number',fontsize=25)
ax.plot(B1num, 'k-')
ax.set_xlim([0,epsmx-1])
ax.set_ylim([0,1.1*max(B1num)])
f.tight_layout()
plt.savefig(path2+'/B1_lineplot'+ext, interpolation='none')
plt.close()


B0Barcode.sort(key=lambda x: (x[1],x[0]))
B1Barcode.sort(key=lambda x: (x[1],x[2],x[0]))

B0Barcodesorted=[]
for n_i,i in enumerate(B0Barcode):
    if n_i==0:
        l1=[]
        l1.append(i[0])
    else:
        if i[1]==B0Barcode[n_i-1][1]:
            l1.append(i[0])
        else:
            B0Barcodesorted.append(l1)
            l1=[]
            l1.append(i[0])
B0Barcodesorted.sort(key=len, reverse=True)


##Find height of B1 barcode representation
B1Height=0
ind=0
B1Barcodesorted=[]
l1=[]
for n_i,i in enumerate(B1Barcode):
    if n_i==0:
        if i[1]!=B1Barcode[ind][1]:
            ##set start of bar to change of index
            ind=n_i
            ##iterate the height
            B1Height+=1
            l1=[]
            l1.append(i[0])
    else:
        ##if a change in barcode identifier
        if i[1]!=B1Barcode[ind][1] or i[0]<B1Barcode[n_i-1][0]:
            ##set start of bar to change of index
            ind=n_i
            ##iterate the height
            B1Height+=1
            B1Barcodesorted.append(l1)
            l1=[]
            l1.append(i[0])
        else:
            l1.append(i[0])
    if n_i==len(B1Barcode):
        B1Barcodesorted.append(l1)
B1Barcodesorted.sort(key=len, reverse=True)


B0Canvas=np.full(((Cellnum+2)*5,epsmx*5+epsstp), False, dtype=bool)
B0CanvasSorted=np.full(((Cellnum+2)*5,epsmx*5+epsstp), False, dtype=bool)
B1Canvas=np.full(((B1Height+2)*5,epsmx*5+epsstp), False, dtype=bool)
B1CanvasSorted=np.full(((B1Height+2)*5,epsmx*5+epsstp), False, dtype=bool)
GapCanvas=np.full((10,epsmx*5+epsstp), True, dtype=bool)


##Generate B0Barcode
ind=0
start=B0Barcode[ind][0]
indheight=1
for n,i in enumerate(B0Barcode):
    ##if a change in barcode identifier
    if i[1]!=B0Barcode[ind][1]:
        ##set bar limits
        xs,xe,ys,ye=max(B0Barcode[ind][0]*5-0.5*epsstp,0),B0Barcode[n-1][0]*5+0.5*epsstp,indheight*5+2,indheight*5+3
        ##draw the rectange
        rr,cc=skdr.rectangle(start=(ys,xs),end=(ye,xe),shape=B0Canvas.shape)
        ##pop this on the canvas
        B0Canvas[rr.astype(int),cc.astype(int)]=True
        ##set start of bar to change of index
        ind=n
        ##iterate the height
        indheight+=1
    else:
        pass
    
    
##Generate B0SortedBarcode
for n,i in enumerate(B0Barcodesorted):
    ##set bar limits
    xs,xe,ys,ye=max(i[0]*5-0.5*epsstp,0),i[-1]*5+0.5*epsstp,(n+1)*5+2,(n+1)*5+3
    ##draw the rectange
    rr,cc=skdr.rectangle(start=(ys,xs),end=(ye,xe),shape=B0Canvas.shape)
    ##pop this on the canvas
    B0CanvasSorted[rr.astype(int),cc.astype(int)]=True

##Generate B1SortedBarcode
for n,i in enumerate(B1Barcodesorted):
    ##set bar limits
    xs,xe,ys,ye=max(i[0]*5-0.5*epsstp,0),i[-1]*5+0.5*epsstp,(n+1)*5+2,(n+1)*5+3
    ##draw the rectange
    rr,cc=skdr.rectangle(start=(ys,xs),end=(ye,xe),shape=B1Canvas.shape)
    ##pop this on the canvas
    B1CanvasSorted[rr.astype(int),cc.astype(int)]=True

##Generate B1Barcode
ind=0
start=B1Barcode[ind][0]
indheight=1
for n,i in enumerate(B1Barcode):
    if n==0:
        if i[1]!=B1Barcode[ind][1]:
            ##set bar limits
            xs,xe,ys,ye=max(B1Barcode[ind][0]*5-0.5*epsstp,0),B1Barcode[n-1][0]*5+0.5*epsstp,indheight*5+2,indheight*5+3
            ##draw the rectange
            rr,cc=skdr.rectangle(start=(ys,xs),end=(ye,xe),shape=B1Canvas.shape)
            ##pop this on the canvas
            B1Canvas[rr.astype(int),cc.astype(int)]=True
            ##set start of bar to change of index
            ind=n
            ##iterate the height
            indheight+=1
    else:
        ##if a change in barcode identifier
        if i[1]!=B1Barcode[ind][1] or i[0]<B1Barcode[n-1][0]:
            ##set bar limits
            xs,xe,ys,ye=max(B1Barcode[ind][0]*5-0.5*epsstp,0),B1Barcode[n-1][0]*5+0.5*epsstp,indheight*5+2,indheight*5+3
            ##draw the rectange
            rr,cc=skdr.rectangle(start=(ys,xs),end=(ye,xe),shape=B1Canvas.shape)
            ##pop this on the canvas
            B1Canvas[rr.astype(int),cc.astype(int)]=True
            ##set start of bar to change of index
            ind=n
            ##iterate the height
            indheight+=1
        else:
            pass


joint=np.concatenate((B0Canvas,GapCanvas,B1Canvas), axis=0)

xtres=10

##plot the barcode
f, ax = plt.subplots()
ax.set_title(r'B0Barcode')
ax.set_xlabel(r'Epsilon Value ')
ax.set_ylabel(r'B0 Components', rotation='vertical')
ax.set_yticks([])
ax.set_xticks( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True)*5 )
ax.set_xticklabels( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True,dtype=int) )
ax.tick_params(axis='x', labelsize=6)
ax.imshow(B0Canvas, interpolation='none',cmap='Greys')
plt.savefig(path2+'/B0_Barcode'+ext, interpolation='none')
plt.close()

##plot the barcode
f, ax = plt.subplots()
ax.set_title(r'B0Barcode')
ax.set_xlabel(r'Epsilon Value ')
ax.set_ylabel(r'B0 Components', rotation='vertical')
ax.set_yticks([])
ax.set_xticks( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True)*5 )
ax.set_xticklabels( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True,dtype=int) )
ax.tick_params(axis='x', labelsize=6)
ax.imshow(B0CanvasSorted, interpolation='none',cmap='Greys')
plt.savefig(path2+'/B0_Barcode_Sorted'+ext, interpolation='none')
plt.close()

f, ax = plt.subplots()
ax.set_title(r'B1Barcode')
ax.set_xlabel(r'Epsilon Value ')
ax.set_ylabel(r'B1 Components', rotation='vertical')
ax.set_yticks([])
ax.set_xticks( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True)*5 )
ax.set_xticklabels( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True,dtype=int) )
ax.tick_params(axis='x', labelsize=6)
ax.imshow(B1Canvas, interpolation='none',cmap='Greys')
plt.savefig(path2+'/B1_Barcode'+ext, interpolation='none')
plt.close()

f, ax = plt.subplots()
ax.set_title(r'B1Barcode')
ax.set_xlabel(r'Epsilon Value ')
ax.set_ylabel(r'B1 Components', rotation='vertical')
ax.set_yticks([])
ax.set_xticks( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True)*5 )
ax.set_xticklabels( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True,dtype=int) )
ax.tick_params(axis='x', labelsize=6)
ax.imshow(B1CanvasSorted, interpolation='none',cmap='Greys')
plt.savefig(path2+'/B1_Barcode_Sorted'+ext, interpolation='none')
plt.close()



f, ax = plt.subplots()
ax.set_title(r'Overall Barcode')
ax.set_xlabel(r'Epsilon Value ')
ax.set_ylabel(r'', rotation='vertical')
ax.set_yticks([])
ax.set_xticks( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True)*5 )
ax.set_xticklabels( np.linspace(start=0,stop=epsmx, num=xtres,endpoint=True,dtype=int) )
ax.tick_params(axis='x', labelsize=6)
ax.imshow(joint, interpolation='none',cmap='Greys')
plt.savefig(path2+'/Overall_Barcode'+ext, interpolation='none')
plt.close()