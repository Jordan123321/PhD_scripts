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
path2=path0+ '/Output_Barcode/Chap2'
ext='.eps'



##Creates list of paths to diabetes patients
lst0= ["T1D/"+f for f in os.listdir(path0+"/T1D")if (f.startswith('T'))]
lst1= ["T2D/"+f for f in os.listdir(path0+"/T2D")if (f.startswith('T'))]
lst2= ["Young_Onset/" + f for f in os.listdir(path0+"/Young_Onset" ) if ('CONTROL' in f.upper()) or ('CASE') in f.upper()]




res=21
SSpace=np.linspace(start=0,stop=20,num=res)
cval=7


CellTot=np.zeros(res)

CellTotYO=np.zeros(res)
CellTotYOCTL=np.zeros(res)


CellTotT1D=np.zeros(res)
CellTotT1DCTL=np.zeros(res)

CellTotT2D=np.zeros(res)
CellTotT2DCTL=np.zeros(res)

Totnum=0
YOnum=0
YOCTLnum=0
T1Dnum=0
T1DCTLnum=0
T2Dnum=0
T2DCTLnum=0

def RngPlot(Lnlist,namelist,graphname,varname):
    plt.figure(figsize=(20 ,20))
    gs = GridSpec(4,2,width_ratios=[1,1],height_ratios=[1,1,1,1],wspace=0.3,hspace=0.3)
    axlist=[]
    axlist.append( plt.subplot(gs[0,:]))
    ymax=0
    for i in Lnlist:
        ymax=max(ymax,np.amax(i))
    for i in [1,2,3]:
        axlist.append(plt.subplot(gs[i,0]))
        axlist.append(plt.subplot(gs[i,1]))
    for n, (i,j,k) in enumerate(zip(Lnlist,axlist,namelist)):
        j.plot(SSpace,i, 'k-')
        j.set_title(k)
        j.set_ylabel(r'Mean ' + graphname, rotation='vertical',fontsize=25)
        j.set_xlabel(varname + r' Parameter Value',fontsize=25)
        j.set_xlim(SSpace[0],SSpace[-1])
        j.set_ylim(0,ymax)
        j.axvline(cval,ls='--',c='r')
        j.tick_params(labelsize=20)
        #j.set_xticklabels(["{:.2e}".format(t) for t in j.get_xticks()])
    plt.tight_layout()
    plt.savefig(path2+'/'+varname+'_'+graphname+ext, interpolation='none')
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
        opathsSST=opaths+'/SST'
        if not os.path.exists(opathsSST):
            os.makedirs(opathsSST)
        ##Get rid of the backround blood
        Sim1=Sim[:,:,0]+Sim[:,:,1]+Sim[:,:,2]

        ##Take the Laplacian of the Stomatostatin
        ##Get rid of the scale bar
        DiffSim1=filt.laplace(Sim1)
        ##Get rid of scale bar
        DiffSim1[-50:,-200:]=0
        CellTemp=np.zeros(res)

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
        for ssn,ss in enumerate(SSpace):
            ##Erode the image, to disentangle joint up DAPI stains
            NucPos5=skmorph.erosion(NucPos4, selem=skmorph.disk(ss))
        
            ##Label and number the nuclei
            NucLabs,Cellnum = skmeas.label(NucPos5, return_num=1)
            CellTemp[ssn]=Cellnum
        print("done for scan")
        CellTot+=CellTemp
        Totnum+=1
        if 'YOUNG' in f0.upper():
            if 'CONTROL' in f0.upper():
                CellTotYOCTL+=CellTemp
                YOCTLnum+=1
            else:
                CellTotYO+=CellTemp
                YOnum+=1
        elif 'T2D' in f0:
            if 'CONTROL' in f0.upper():
                CellTotT2DCTL+=CellTemp
                T2DCTLnum+=1
            else:
                CellTotT2D+=CellTemp
                T2Dnum+=1
        else:
            if 'CONTROL' in f0.upper():
                CellTotT1DCTL+=CellTemp
                T1DCTLnum+=1
            else:
                CellTotT1D+=CellTemp
                T1Dnum+=1
    print("done for patient")

CellTotnorm=CellTot/Totnum


CellTotYOnorm=CellTotYO/YOnum
CellTotYOCTLnorm=CellTotYOCTL/YOCTLnum

CellTotT1Dnorm=CellTotT1D/T1Dnum
CellTotT1DCTLnorm=CellTotT1DCTL/T1DCTLnum


CellTotT2Dnorm=CellTotT2D/T2Dnum
CellTotT2DCTLnorm=CellTotT2DCTL/T2DCTLnum

stringlist= ['Overall','YO','YO_Control','T1D','T1D_Control','T2D','T2D_Control']

RngPlot([CellTotnorm,
         CellTotYOnorm, CellTotYOCTLnorm,
         CellTotT1Dnorm, CellTotT1DCTLnorm,
         CellTotT2Dnorm, CellTotT2DCTLnorm],
         stringlist,'CellCount', 'Erosion Disc Size Block Size')
