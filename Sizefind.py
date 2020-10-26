#!/usr/bin/env python3

##This is a flag.
Silent=False

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
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

##Scipy contains some elements used in image processing that are required here
#import scipy.spatial as scspat
#import scipy.ndimage as ndim
##This imports the statistics module from scipy
#import scipy.stats as scstats

##Scikit image is a package which contains a large number of image processing functions
import skimage.io as skio
import skimage.morphology as skmorph
import skimage.filters as filt
import skimage.measure as skmeas
#import skimage.feature as skfea
#import skimage.segmentation as skseg
import skimage.draw as skdr
#import skimage.color as skcol

#import csv
import math
#import random as rd

##Imports my Voronoi Split Algorithm into a Module
#import VorSplit
from math import log10, floor

#import seaborn as sns
import basic_units

#import itertools

Rat=25/(92*1000000)


def normedmean(arr,ax):
    mean = np.mean(arr)
    variance = np.var(arr)
    sigma = np.sqrt(variance)
    x = np.linspace(min(arr), max(arr), 100)
    ax.plot(x, mlab.normpdf(x, mean, sigma))

def ExCirc(f_image,im_cent):
    #print(im_cent)
    #f, ax = plt.subplots()
    #ax.set_title(r'Test Islet Shape')
    #ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
    #ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
    #ax=plt.matshow(f_image, interpolation='none')
    #plt.plot(im_cent[0],im_cent[1], marker='o', markersize=3, color="red")
    #plt.savefig('Test'+ext, interpolation='none')
    #plt.close()


    maxind=math.ceil(math.sqrt(f_image.shape[0]*f_image.shape[0]+f_image.shape[1]*f_image.shape[1])/2)
    minind=math.floor(max(f_image.shape)/2)
    r_ex=minind
    for i in range(maxind,minind,-1):
        cicim=np.ones_like(f_image,dtype=bool)
        rr,cc=skdr.circle(im_cent[0],im_cent[1],i,f_image.shape)
        cicim[rr,cc]=0

        if (np.logical_and(cicim,f_image).any()):
            r_ex=i
            break
    #print(r_ex)
    return r_ex

def InCirc(f_image,im_cent):
    #print(im_cent)
    #f, ax = plt.subplots()
    #ax.set_title(r'Test Islet Shape')
    #ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
    #ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
    #ax=plt.matshow(f_image, interpolation='none')
    #plt.plot(im_cent[0],im_cent[1], marker='o', markersize=3, color="red")
    #plt.savefig('Test'+ext, interpolation='none')
    #plt.close()

    maxind=math.ceil(min(f_image.shape)/2)
    r_in=0
    for i in range(0,maxind):
        cicim=np.zeros_like(f_image,dtype=bool)
        rr,cc=skdr.circle(im_cent[0],im_cent[1],i,f_image.shape)
        cicim[rr,cc]=1

        if (np.logical_and(cicim,np.logical_not(f_image)).any()):
            r_in=i
            break
    return r_in

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)


##Defines paths to my directories and save file locations
path1=os.getcwd()
path0=os.path.dirname(path1)
path2=path0+ '/Output_Isize'
ext='.eps'



##Creates list of paths to diabetes patients
lst0= ["T1D/"+f for f in os.listdir(path0+"/T1D")if (f.startswith('T'))]
lst1= ["T2D/"+f for f in os.listdir(path0+"/T2D")if (f.startswith('T'))]
lst2= ["Young_Onset/" + f for f in os.listdir(path0+"/Young_Onset" ) if ('CONTROL' in f.upper()) or ('CASE') in f.upper()]






AvradArr=[]
ExRadArr=[]
InRadArr=[]
DiameterArr=[]
CircRatArr=[]
AreaArr=[]
PerimArr=[]
IsoparametricratioArr=[]
PerimRatArr=[]

AvradT1DArr=[]
ExRadT1DArr=[]
InRadT1DArr=[]
DiameterT1DArr=[]
CircRatT1DArr=[]
AreaT1DArr=[]
PerimT1DArr=[]
IsoparametricratioT1DArr=[]
PerimRatT1DArr=[]

AvradT1DCTLArr=[]
ExRadT1DCTLArr=[]
InRadT1DCTLArr=[]
DiameterT1DCTLArr=[]
CircRatT1DCTLArr=[]
AreaT1DCTLArr=[]
PerimT1DCTLArr=[]
IsoparametricratioT1DCTLArr=[]
PerimRatT1DCTLArr=[]


AvradT2DArr=[]
ExRadT2DArr=[]
InRadT2DArr=[]
DiameterT2DArr=[]
CircRatT2DArr=[]
AreaT2DArr=[]
PerimT2DArr=[]
IsoparametricratioT2DArr=[]
PerimRatT2DArr=[]

AvradT2DCTLArr=[]
ExRadT2DCTLArr=[]
InRadT2DCTLArr=[]
DiameterT2DCTLArr=[]
CircRatT2DCTLArr=[]
AreaT2DCTLArr=[]
PerimT2DCTLArr=[]
IsoparametricratioT2DCTLArr=[]
PerimRatT2DCTLArr=[]


AvradYOArr=[]
ExRadYOArr=[]
InRadYOArr=[]
DiameterYOArr=[]
CircRatYOArr=[]
AreaYOArr=[]
PerimYOArr=[]
IsoparametricratioYOArr=[]
PerimRatYOArr=[]

AvradYOCTLArr=[]
ExRadYOCTLArr=[]
InRadYOCTLArr=[]
DiameterYOCTLArr=[]
CircRatYOCTLArr=[]
AreaYOCTLArr=[]
PerimYOCTLArr=[]
IsoparametricratioYOCTLArr=[]
PerimRatYOCTLArr=[]

def Histplot(Histlist,namelist,name,units=None):
    plt.figure(figsize=(20 ,20))
    gs = GridSpec(4,2,width_ratios=[1,1],height_ratios=[1,1,1,1],wspace=0.3,hspace=0.3)
    axlist=[]
    axlist.append( plt.subplot(gs[0,:]))
    view_rat=0.1
    for i in [1,2,3]:
        axlist.append(plt.subplot(gs[i,0]))
        axlist.append(plt.subplot(gs[i,1]))
    for n, (i,j,k) in enumerate(zip(Histlist,axlist,namelist)):
        if n ==0:
            y,x,_=j.hist(i,bins=50)
            xmin=max(0,x.min()-view_rat*(x.max()-x.min()))
            xmax=x.max()+view_rat*(x.max()-x.min())
            ymax=(1+view_rat)*y.max()
        elif n % 2 == 1:
            y,x,_=j.hist(i,bins=10)
            y1,x1,_=axlist[n+1].hist(Histlist[n+1],bins=10)
            ymax=max((1+view_rat)*y.max(),(1+view_rat)*y1.max())
        j.set_title(k + r', n = ' +str(len(i)))
        j.set_ylabel(r'frequency count', rotation='vertical')
        if units:
            j.set_xlabel(name+r' of the islet (' +units+')')
        else:
            j.set_xlabel(name+r' of the islet')
        j.set_xlim(xmin,xmax)
        j.set_ylim(0,ymax)
        j.set_xticklabels(["{:.2e}".format(t) for t in j.get_xticks()])
    #plt.tight_layout()
    plt.savefig(path2+'/'+name+ext, interpolation='none')
    plt.close()


##This Loops Through the Patients
for PNum, f0 in enumerate(lst0+lst1+lst2):
    ##Sets the input and output paths for the count
    ppaths=path0+'/'+f0
    opaths=path2+'/'+f0
    #print(ppaths)



    if 'YOUNG' in ppaths.upper() and 'CONTROL' in ppaths.upper():
        DAPIlst=[]
        GCGlst=[]
        INSlst=[]
        SSTlst=[]
        OALLlst=[]
        for i in [f for f in os.listdir(ppaths)if f.startswith('Image')]:
            DAPIlst+=[i+'/'+f for f in os.listdir(ppaths+'/'+i) if (f.endswith('.tif') and 'DAPI' in f.upper())]
            GCGlst+=[i+'/'+f for f in os.listdir(ppaths+'/'+i) if (f.endswith('.tif') and 'GCG' in f.upper())]
            INSlst+=[i+'/'+f for f in os.listdir(ppaths+'/'+i) if (f.endswith('.tif') and 'INS' in f.upper())]
            SSTlst+=[i+'/'+f for f in os.listdir(ppaths+'/'+i) if (f.endswith('.tif') and 'SST' in f.upper())]
            OALLlst+=[i+'/'+f for f in os.listdir(ppaths+'/'+i) if (f.endswith('.tif') and 'OVERLAY' in f.upper())]
    elif 'YOUNG' in ppaths.upper() and 'CASE' in ppaths.upper():
        DAPIlst=[]
        GCGlst=[]
        INSlst=[]
        SSTlst=[]
        OALLlst=[]
        for i in [f for f in os.listdir(ppaths+'/ICI images/')if f.startswith('Image')]:
            DAPIlst+=['/ICI images/'+i+'/'+f for f in os.listdir(ppaths+'/ICI images/'+i) if (f.endswith('.tif') and 'DAPI' in f.upper())]
            GCGlst+=['/ICI images/'+i+'/'+f for f in os.listdir(ppaths+'/ICI images/'+i) if (f.endswith('.tif') and 'GCG' in f.upper())]
            INSlst+=['/ICI images/'+i+'/'+f for f in os.listdir(ppaths+'/ICI images/'+i) if (f.endswith('.tif') and 'INS' in f.upper())]
            SSTlst+=['/ICI images/'+i+'/'+f for f in os.listdir(ppaths+'/ICI images/'+i) if (f.endswith('.tif') and 'SST' in f.upper())]
            OALLlst+=['/ICI images/'+i+'/'+f for f in os.listdir(ppaths+'/ICI images/'+i) if (f.endswith('.tif') and 'OVERLAY' in f.upper())]

#        for i in [f for f in os.listdir(ppaths+'/IDI images/')if f.startswith('Image')]:
#            DAPIlst+=['/IDI images/'+i+'/'+f for f in os.listdir(ppaths+'/IDI images/'+i) if (f.endswith('.tif') and 'DAPI' in f.upper())]
#            GCGlst+=['/IDI images/'+i+'/'+f for f in os.listdir(ppaths+'/IDI images/'+i) if (f.endswith('.tif') and 'GCG' in f.upper())]
#            INSlst+=['/IDI images/'+i+'/'+f for f in os.listdir(ppaths+'/IDI images/'+i) if (f.endswith('.tif') and 'INS' in f.upper())]
#            SSTlst+=['/IDI images/'+i+'/'+f for f in os.listdir(ppaths+'/IDI images/'+i) if (f.endswith('.tif') and 'SST' in f.upper())]
#            OALLlst+=['/IDI images/'+i+'/'+f for f in os.listdir(ppaths+'/IDI images/'+i) if (f.endswith('.tif') and 'OVERLAY' in f.upper())]

    else:
        ##These return lists of the scans that have the different hormone stains
        DAPIlst=[f for f in os.listdir(ppaths)if (f.endswith('.tif') and 'DAPI' in f.upper())]
        GCGlst=[f for f in os.listdir(ppaths)if (f.endswith('.tif') and 'GCG' in f.upper())]
        INSlst=[f for f in os.listdir(ppaths)if (f.endswith('.tif') and 'INS' in f.upper())]
        SSTlst=[f for f in os.listdir(ppaths)if (f.endswith('.tif') and 'SST' in f.upper())]
        OALLlst=[f for f in os.listdir(ppaths)if (f.endswith('.tif') and 'OVERLAY' in f.upper())]
        ##Creates string detailing the patients
    pnt=f0.split("/")[-1]
    #print(OALLlst)
    #print('\n \n')
    #continue
    DAPIlst.sort(key=lambda x: x.split()[-1])
    GCGlst.sort(key=lambda x: x.split()[-1])
    INSlst.sort(key=lambda x: x.split()[-1])
    SSTlst.sort(key=lambda x: x.split()[-1])
    OALLlst.sort(key=lambda x: x.split()[-1])
    ##Makes sure all the lists are consistent
    if (any(len(lst) != len(DAPIlst) for lst in [GCGlst, INSlst, SSTlst])) or (DAPIlst == []):
        print(f0)
        print([[lst, len(lst)] for lst in [DAPIlst, GCGlst, INSlst, SSTlst, OALLlst]])
        continue
    if not os.path.exists(opaths):
        os.makedirs(opaths)




    ##This loops through the Scans for the various patients
    for num,(D,G,I,S,O) in enumerate(zip(DAPIlst,GCGlst, INSlst, SSTlst,OALLlst)):
        tst=D.split()[-1]
        if ((D.split()[-1] != G.split()[-1]) or (D.split()[-1] != I.split()[-1]) or (D.split()[-1] != I.split()[-1]) or (D.split()[-1] != O.split()[-1]) ):
            print(D.split()[-1] ,G.split()[-1] ,I.split()[-1] ,S.split()[-1] ,O.split()[-1] )
            continue
        numbr=tst[:-4]
        #print(tst, numbr)

        ##Set the scan Beta, Delta and Alpha cell count to zero, and make an image number string
        SmI1=0
        SmS1=0
        SmG1=0
        Dim=skio.imread(ppaths+'/'+D)
        Gim=skio.imread(ppaths+'/'+G)
        Iim=skio.imread(ppaths+'/'+I)
        Sim=skio.imread(ppaths+'/'+S)
        Oim=skio.imread(ppaths+'/'+O)
        #print(Dim.shape)
        ##If silet is off, then plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Original overlay. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(Oim, interpolation='none')
            plt.savefig(opaths+'/'+numbr+'_Orig_Img'+ext, interpolation='none')
            plt.close()
        opathsSST=opaths+'/SST'
        if not os.path.exists(opathsSST):
            os.makedirs(opathsSST)
        ##Get rid of the backround blood
        Sim1=Sim[:,:,0]+Sim[:,:,1]+Sim[:,:,2]

        ##Take the Laplacian of the Stomatostatin
        ##Get rid of the scale bar
        DiffSim1=filt.laplace(Sim1)
        ##If silet is off, then plot this out
        ##Get rid of scale bar
        DiffSim1[-50:,-200:]=0
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Laplacian map. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(DiffSim1, interpolation='none')
            plt.savefig(opathsSST+'/'+numbr+'_1_Lplce'+ext, interpolation='none')
            plt.close()
        ##Smooth it
        DiffSim2=filt.gaussian(DiffSim1,1)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Laplacian map+Gauss blur. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(DiffSim2, interpolation='none')
            plt.savefig(opathsSST+'/'+numbr+'_2_LPCE_GSS'+ext, interpolation='none')
            plt.close()
        #Filter it
        DiffSim3=DiffSim2>filt.threshold_triangle(DiffSim2)



        ##If silet is off, then plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'SST Mask. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(DiffSim3, interpolation='none')
            plt.savefig(opathsSST+'/'+numbr+'_3_BloodFilt'+ext, interpolation='none')
            plt.close()


        ##Mask the original image with the blood removed
        Stcked=np.stack([DiffSim3,DiffSim3,DiffSim3],axis=2)
        Sim=np.multiply(Stcked,Sim)
        ##If silet is off, then plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Filt_SST. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(Sim, interpolation='none')
            plt.savefig(opathsSST+'/'+numbr+'_4_SomatoStatin'+ext, interpolation='none')
            plt.close()
        opathsIslt=opaths+'/Islt'
        if not os.path.exists(opathsIslt):
            os.makedirs(opathsIslt)
        ##First thing is to find the islet shape. This can be done by adding the non Dapi, smoothing and thresholding
        ##Add together the three scans
        IsltShp=Gim+Iim+Sim
        ##Add together the red blue and green to flatten the array
        IsltShp0=IsltShp[:,:,0]+IsltShp[:,:,1]+IsltShp[:,:,2]
        ##Get rid of scale bar
        IsltShp0[-50:,-200:]=0
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Hormone_Grayscale. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ irection (pixels)', rotation='vertical')
            ax=plt.imshow(IsltShp0, interpolation='none')
            plt.savefig(opathsIslt+'/'+numbr+'_1_GScale'+ext, interpolation='none')
            plt.close()
        ##Gaussian smooth followed by a triangle filter to determine the boundary of the islets
        IsltShp1=filt.gaussian(IsltShp0,10)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Gaussian Islet shape. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(IsltShp1, interpolation='none')
            plt.savefig(opathsIslt+'/'+numbr+'_2_G_smooth'+ext, interpolation='none')
            plt.close()
        IsltShp2=IsltShp1>filt.threshold_triangle(IsltShp1)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Triangle Filt shape. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(IsltShp2, interpolation='none')
            plt.savefig(opathsIslt+'/'+numbr+'_3_Tfilt'+ext, interpolation='none')
            plt.close()
        ##Get rid of any small objects that may yet exist
        IsltShp3=skmorph.remove_small_holes(IsltShp2, connectivity=1, area_threshold=1000)
        IsltShp4=skmorph.remove_small_objects(IsltShp3, connectivity=1, min_size=10000)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Objects_removed. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(IsltShp4, interpolation='none')
            plt.savefig(opathsIslt+'/'+numbr+'_4_obj_rem'+ext, interpolation='none')
            plt.close()

        IsltLabs = skmeas.label(IsltShp4, return_num=0)
        props=skmeas.regionprops(IsltLabs)
        Oim1=Oim.copy()
        for p in props:
            if(
                p.bbox[0]==0 or
                p.bbox[1]==0 or
                p.bbox[2]==IsltShp4.shape[0] or
                 p.bbox[3]==IsltShp4.shape[1]
            ):
                continue
            else:
                Avrad=(0.5*p.minor_axis_length+0.5*p.major_axis_length)*Rat
                ExRad=ExCirc(p.image,(p.centroid[0]-p.bbox[0],p.centroid[1]-p.bbox[1]))
                InRad=InCirc(p.image,(p.centroid[0]-p.bbox[0],p.centroid[1]-p.bbox[1]))
                Diameter=2*ExRad*Rat
                CircRat=InRad/ExRad
                Area=p.area*Rat*Rat
                Perim=p.perimeter*Rat
                Isoparametricratio=4*math.pi*Area/(Perim*Perim)
                PerimRat=Perim/(2*math.pi*Avrad)
                AvradArr.append(Avrad)
                ExRadArr.append(ExRad)
                InRadArr.append(InRad)
                DiameterArr.append(Diameter)
                CircRatArr.append(CircRat)
                AreaArr.append(Area)
                PerimArr.append(Perim)
                IsoparametricratioArr.append(Isoparametricratio)
                PerimRatArr.append(PerimRat)

                for rad_rng in range(-1,2):
                    if (rad_rng+InRad)>=0:
                        ri, ci = skdr.circle_perimeter(int(p.centroid[0]),int(p.centroid[1]), radius=int(InRad+rad_rng), shape=Oim1.shape)
                    re, ce= skdr.circle_perimeter(int(p.centroid[0]),int(p.centroid[1]),  radius=int(ExRad+rad_rng), shape=Oim1.shape)
                    Oim1[ri,ci,:]=(255,255,0)
                    Oim1[re,ce,:]=(0,255,255)


                if 'YOUNG' in f0.upper():
                    if 'CONTROL' in f0.upper():
                        AvradYOCTLArr.append(Avrad)
                        ExRadYOCTLArr.append(ExRad)
                        InRadYOCTLArr.append(InRad)
                        DiameterYOCTLArr.append(Diameter)
                        CircRatYOCTLArr.append(CircRat)
                        AreaYOCTLArr.append(Area)
                        PerimYOCTLArr.append(Perim)
                        IsoparametricratioYOCTLArr.append(Isoparametricratio)
                        PerimRatYOCTLArr.append(PerimRat)
                    else:
                        AvradYOArr.append(Avrad)
                        ExRadYOArr.append(ExRad)
                        InRadYOArr.append(InRad)
                        DiameterYOArr.append(Diameter)
                        CircRatYOArr.append(CircRat)
                        AreaYOArr.append(Area)
                        PerimYOArr.append(Perim)
                        IsoparametricratioYOArr.append(Isoparametricratio)
                        PerimRatYOArr.append(PerimRat)
                elif 'T2D' in f0.upper():
                    if 'CONTROL' in f0.upper():
                        AvradT2DCTLArr.append(Avrad)
                        ExRadT2DCTLArr.append(ExRad)
                        InRadT2DCTLArr.append(InRad)
                        DiameterT2DCTLArr.append(Diameter)
                        CircRatT2DCTLArr.append(CircRat)
                        AreaT2DCTLArr.append(Area)
                        PerimT2DCTLArr.append(Perim)
                        IsoparametricratioT2DCTLArr.append(Isoparametricratio)
                        PerimRatT2DCTLArr.append(PerimRat)
                    else:
                        AvradT2DArr.append(Avrad)
                        ExRadT2DArr.append(ExRad)
                        InRadT2DArr.append(InRad)
                        DiameterT2DArr.append(Diameter)
                        CircRatT2DArr.append(CircRat)
                        AreaT2DArr.append(Area)
                        PerimT2DArr.append(Perim)
                        IsoparametricratioT2DArr.append(Isoparametricratio)
                        PerimRatT2DArr.append(PerimRat)
                else:
                    if 'CONTROL' in f0.upper():
                        AvradT1DCTLArr.append(Avrad)
                        ExRadT1DCTLArr.append(ExRad)
                        InRadT1DCTLArr.append(InRad)
                        DiameterT1DCTLArr.append(Diameter)
                        CircRatT1DCTLArr.append(CircRat)
                        AreaT1DCTLArr.append(Area)
                        PerimT1DCTLArr.append(Perim)
                        IsoparametricratioT1DCTLArr.append(Isoparametricratio)
                        PerimRatT1DCTLArr.append(PerimRat)
                    else:
                        AvradT1DArr.append(Avrad)
                        ExRadT1DArr.append(ExRad)
                        InRadT1DArr.append(InRad)
                        DiameterT1DArr.append(Diameter)
                        CircRatT1DArr.append(CircRat)
                        AreaT1DArr.append(Area)
                        PerimT1DArr.append(Perim)
                        IsoparametricratioT1DArr.append(Isoparametricratio)
                        PerimRatT1DArr.append(PerimRat)
            opathsOlay=opaths+'/Olay'
            if not os.path.exists(opathsOlay):
                os.makedirs(opathsOlay)
            if Silent == False:
                f, ax = plt.subplots()
                ax.set_title(r'ExCirc/InCirc Overlay. Islet is: %s and patient is %s' %(numbr,pnt))
                ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
                ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
                ax=plt.imshow(Oim1, interpolation='none')
                plt.savefig(opathsOlay+'/'+numbr+'In_Ex_Ovlay'+ext, interpolation='none')
                plt.close()

    print('Done for patient %s, %s of %s' %(f0, str(PNum),str(len(lst0+lst2+lst1)-1) ) )

stringlist= ['Overall','YO','YO_Control','T1D','T1D_Control','T2D','T2D_Control']
Histplot([AvradArr,
                  AvradYOArr, AvradYOCTLArr,
                  AvradT1DArr, AvradT1DCTLArr,
                  AvradT2DArr, AvradT2DCTLArr],
                  stringlist,'Average_Radius', 'Metres')

Histplot([DiameterArr,
                  DiameterYOArr, DiameterYOCTLArr,
                  DiameterT1DArr, DiameterT1DCTLArr,
                  DiameterT2DArr, DiameterT2DCTLArr],
                  stringlist,'Diameter', 'Metres')

Histplot([CircRatArr,
                  CircRatYOArr, CircRatYOCTLArr,
                  CircRatT1DArr, CircRatT1DCTLArr,
                  CircRatT2DArr, CircRatT2DCTLArr],
                  stringlist,'Ratio_Between_In_and_Excirlces')

Histplot([AreaArr,
                  AreaYOArr, AreaYOCTLArr,
                  AreaT1DArr, AreaT1DCTLArr,
                  AreaT2DArr, AreaT2DCTLArr],
                  stringlist,'Area','Metres^2')

Histplot([PerimArr,
                  PerimYOArr, PerimYOCTLArr,
                  PerimT1DArr, PerimT1DCTLArr,
                  PerimT2DArr, PerimT2DCTLArr],
                  stringlist,'Perimeter','Metres')

Histplot([IsoparametricratioArr,
                  IsoparametricratioYOArr, IsoparametricratioYOCTLArr,
                  IsoparametricratioT1DArr, IsoparametricratioT1DCTLArr,
                  IsoparametricratioT2DArr, IsoparametricratioT2DCTLArr],
                  stringlist,'Isoparametric_ratio')


Histplot([PerimRatArr,
                  PerimRatYOArr, PerimRatYOCTLArr,
                  PerimRatT1DArr, PerimRatT1DCTLArr,
                  PerimRatT2DArr, PerimRatT2DCTLArr],
                  stringlist,'Perimeter_to_radius_ratio')
