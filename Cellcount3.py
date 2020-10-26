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
from matplotlib import rc
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)
from matplotlib.backends.backend_pdf import PdfPages


plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['font.family'] = 'Calibri'

##Scipy contains some elements used in image processing that are required here
import scipy.spatial as scspat
import scipy.ndimage as ndim
##This imports the statistics module from scipy
import scipy.stats as scstats

##Scikit image is a package which contains a large number of image processing functions
import skimage.io as skio
import skimage.morphology as skmorph
import skimage.filters as filt
import skimage.measure as skmeas
import skimage.feature as skfea
import skimage.segmentation as skseg
import skimage.draw as skdr
import skimage.color as skcol

import csv
import math
import random as rd

##Imports my Voronoi Split Algorithm into a Module
import VorSplit
from math import log10, floor

import seaborn as sns

import itertools


def normedmean(arr,ax):
    mean = np.mean(arr)
    variance = np.var(arr)
    sigma = np.sqrt(variance)
    x = np.linspace(min(arr), max(arr), 100)
    ax.plot(x, mlab.normpdf(x, mean, sigma))

def round_sig(x, sig=3):
    return round(x, sig-int(floor(log10(abs(x))))-1)


##Defines paths to my directories and save file locations
path1=os.getcwd()
path0=os.path.dirname(path1)
path2=path0+ '/Output'
ext='.svg'



##Creates list of paths to diabetes patients
lst0= ["T1D/"+f for f in os.listdir(path0+"/T1D")if (f.startswith('T'))]
lst1= ["T2D/"+f for f in os.listdir(path0+"/T2D")if (f.startswith('T'))]
lst2= ["Young_Onset/" + f for f in os.listdir(path0+"/Young_Onset" ) if ('CONTROL' in f.upper()) or ('CASE') in f.upper()]




##Zeroes all the sums of the types of patients (Sm)=sum     (I)=Beta    (T1D)=Type 1 Diabetes       (CTL)=Control
##                                                          (S)=Delta   (T2D)=Type 2 Diabetes
##                                                          (G)=Alpha   (YO)=Young Onset

SmIYO=0
SmIT1=0
SmIT2=0
SmSYO=0
SmST1=0
SmST2=0
SmGYO=0
SmGT1=0
SmGT2=0

SmIYOCTL=0
SmIT1CTL=0
SmIT2CTL=0
SmSYOCTL=0
SmST1CTL=0
SmST2CTL=0
SmGYOCTL=0
SmGT1CTL=0
SmGT2CTL=0

##Create empty arrays of the percentage of alpha and delta cells that are delta cells
ArrYO=[]
ArrT1D=[]
ArrT2D=[]

ArrIYO=[]
ArrSYO=[]
ArrGYO=[]

ArrIT1D=[]
ArrST1D=[]
ArrGT1D=[]

ArrIT2D=[]
ArrST2D=[]
ArrGT2D=[]

ArrYO0=[]
ArrT1D0=[]
ArrT2D0=[]

ArrIYO0=[]
ArrSYO0=[]
ArrGYO0=[]

ArrIT1D0=[]
ArrST1D0=[]
ArrGT1D0=[]

ArrIT2D0=[]
ArrST2D0=[]
ArrGT2D0=[]

YO_sm=[]
T1D_sm=[]
T2D_sm=[]



ArrYOCTL=[]
ArrT1DCTL=[]
ArrT2DCTL=[]

ArrIYOCTL=[]
ArrSYOCTL=[]
ArrGYOCTL=[]

ArrIT1DCTL=[]
ArrST1DCTL=[]
ArrGT1DCTL=[]

ArrIT2DCTL=[]
ArrST2DCTL=[]
ArrGT2DCTL=[]

ArrYOCTL0=[]
ArrT1DCTL0=[]
ArrT2DCTL0=[]

ArrIYOCTL0=[]
ArrSYOCTL0=[]
ArrGYOCTL0=[]

ArrIT1DCTL0=[]
ArrST1DCTL0=[]
ArrGT1DCTL0=[]

ArrIT2DCTL0=[]
ArrST2DCTL0=[]
ArrGT2DCTL0=[]

YOCTL_sm=[]
T1DCTL_sm=[]
T2DCTL_sm=[]

CCount_T1D0=[]
CCount_T1DCTL0=[]
CCount_T2D0=[]
CCount_T2DCTL0=[]
CCount_YO0=[]
CCount_YOCTL0=[]

PathCW=path2+'/Cell_Counts.csv'
if os.path.exists(PathCW):
    os.remove(PathCW)
with open(PathCW,mode='a+',newline='') as csvfile:
    Wrtr=csv.writer(csvfile, delimiter=',')
    Wrtr.writerow(['Patient Type','Patient Identifier','Islet number','Total Cell Count','Alpha Cell Count', 'Beta Cell Count', 'Delta Cell Count'])


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

    ##Set the patient Beta, Delta and Alpha cell count to zero
    SmI0=0
    SmS0=0
    SmG0=0

    ArrG1=[]
    ArrI1=[]
    ArrS1=[]

    P_sm=[]



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
            ax=plt.imshow(DiffSim1, interpolation='none')
            plt.savefig(opathsSST+'/'+numbr+'_1_Lplce'+ext, interpolation='none')
            plt.close()
        ##Smooth it
        DiffSim2=filt.gaussian(DiffSim1,1)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Laplacian map+Gauss blur. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(DiffSim2, interpolation='none')
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
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
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

        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Histogram of original overlay. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'Pixel Intensity')
            ax.set_ylabel(r'number iof pixels', rotation='vertical')
            ax.hist(IsltShp1.ravel(), bins=IsltShp0.ravel().max())
            plt.savefig(opaths+'/'+numbr+'_1_Hist_Hormone'+ext, interpolation='none')
            plt.close()



        opathsNuc=opaths+'/Nuc'
        if not os.path.exists(opathsNuc):
            os.makedirs(opathsNuc)
        ##Next is to find the Nuclei positions
        ##Flatten the arrays
        Nuc=Dim[:,:,0]+Dim[:,:,1]+Dim[:,:,2]

        ##Get rid of scale bar
        Nuc[-50:,-200:]=0

        ##Local threshold
        NucPos=(Nuc>filt.threshold_local(Nuc, block_size=101))
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Local filter. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(NucPos, interpolation='none')
            plt.savefig(opathsNuc+'/'+numbr+'_1_Local_Filt'+ext, interpolation='none')
            plt.close()

        ##Mask the nuclei with the islet shape
        IsltNuc=np.multiply(NucPos,IsltShp4)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Isltshape mask. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(IsltNuc, interpolation='none')
            plt.savefig(opathsNuc+'/'+numbr+'_2_Masked'+ext, interpolation='none')
            plt.close()
        ##Fill the holes to make contingent nuclei
        NucPos2=ndim.binary_fill_holes(IsltNuc)
        ##Remove small artefacts
        NucPos3=skmorph.remove_small_holes(NucPos2, connectivity=1, area_threshold=1000)
        NucPos4=skmorph.remove_small_objects(NucPos3, connectivity=1, min_size=100)


        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Nuc Small obj rem. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(NucPos4, interpolation='none')
            plt.savefig(opathsNuc+'/'+numbr+'_3_Rem'+ext, interpolation='none')
            plt.close()

        ##Erode the image, to disentangle joint up DAPI stains
        NucPos5=skmorph.erosion(NucPos4, selem=skmorph.disk(9))
        ##If silent is off, plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Nucleus position. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(NucPos5, interpolation='none')
            plt.savefig(opathsNuc+'/'+numbr+'_Nuclei'+ext, interpolation='none')
            plt.close()


        ##Label and number the nuclei
        NucLabs,NucNum = skmeas.label(NucPos5, return_num=1)
        P_sm.append(NucNum)
        ##Properties of regions assigned to an array
        props=skmeas.regionprops(NucLabs)
        ##Make an array to store where the centers are
        NPosArr=np.empty((NucNum,2))
        ##Loop through and save this
        for c, p in enumerate(props):
            NPosArr[c,0],NPosArr[c,1]=p.centroid
        ##Skip the rest of this loop in the unlikely event that there are no nuclei
        if NucNum < 4:
            continue
        ##Do a voronoi split on the Islet and the points of the nuclei centres
        canvas= VorSplit.VorSplt(NPosArr, IsltShp4.copy())

        ##Plot this if silent is off
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Cell Shape. Islet is: %s and patient is %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(canvas, interpolation='none')
            plt.savefig(opaths+'/'+numbr+'_Cells'+ext, interpolation='none',cmap='gray')
            ax.scatter(x=NPosArr[:,1],y=NPosArr[:,0], s=20, c='k')
            plt.savefig(opaths+'/'+numbr+'_Cells_2'+ext, interpolation='none',cmap='gray')
            plt.close()

        ##Labels the parts of the voronoi algorithm for the proporties array. 1 connectivity denotes the single line separating the cells
        CellLab,CellNum = skmeas.label(canvas, return_num=1,connectivity=1,background=0)
        ##Get properties of the regions
        props=skmeas.regionprops(CellLab)
        ##Flatten the scan arrays
        IimF=Iim[:,:,0]+Iim[:,:,1]+Iim[:,:,2]
        GimF=Gim[:,:,0]+Gim[:,:,1]+Gim[:,:,2]
        SimF=Sim[:,:,0]+Sim[:,:,1]+Sim[:,:,2]
        ##Create an 8bit rbg canvas
        Canvas1=np.zeros((Oim.shape[0],Oim.shape[1],3),dtype=np.uint8)
        for n,p in enumerate(props):
            ##Mask the insulin with the cells as defined by the Voronoi algorithm
            CellI=np.multiply(p.image,IimF[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3]])
            CellG=np.multiply(p.image,GimF[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3]])
            CellS=np.multiply(p.image,SimF[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3]])

            ##Calculate the total intensity of the scans within the defined cells
            totI=np.sum(CellI)
            totG=np.sum(CellG)
            totS=3*np.sum(CellS)
            ##Compare these cells and evaluate accordingly
            ##blue is insulin
            if totI>totG and totI>totS:
                Canvas1[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3],2]+=p.image
                SmI0+=1
                SmI1+=1
                if 'YOUNG' in f0.upper():
                    if 'CONTROL' in f0.upper():
                        SmIYOCTL+=1
                    else:
                        SmIYO+=1
                elif 'T2D' in f0:
                    if 'CONTROL' in f0.upper():
                        SmIT2CTL+=1
                    else:
                        SmIT2+=1
                else:
                    if 'CONTROL' in f0.upper():
                        SmIT1CTL+=1
                    else:
                        SmIT1+=1
            ##red is glucagon
            elif totG>totS:
                Canvas1[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3],0]+=p.image
                SmG0+=1
                SmG1+=1
                if 'YOUNG' in f0.upper():
                    if 'CONTROL' in f0.upper():
                        SmGYOCTL+=1
                    else:
                        SmGYO+=1
                elif 'T2D' in f0:
                    if 'CONTROL' in f0.upper():
                        SmGT2CTL+=1
                    else:
                        SmGT2+=1
                else:
                    if 'CONTROL' in f0.upper():
                        SmGT1CTL+=1
                    else:
                        SmGT1+=1
            ##Green is delta cells
            else:
                Canvas1[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3],1]+=p.image
                SmS0+=1
                SmS1+=1
                if 'YOUNG' in f0.upper():
                    if 'CONTROL' in f0.upper():
                        SmSYOCTL+=1
                    else:
                        SmSYO+=1
                elif 'T2D' in f0:
                    if 'CONTROL' in f0.upper():
                        SmST2CTL+=1
                    else:
                        SmST2+=1
                else:
                    if 'CONTROL' in f0.upper():
                        SmST1CTL+=1
                    else:
                        SmST1+=1
        ArrG1.append(SmG1)
        ArrI1.append(SmI1)
        ArrS1.append(SmS1)
        Canvas1=Canvas1*255
        with open(PathCW,mode='a+',newline='') as csvfile:
            Wrtr=csv.writer(csvfile, delimiter=',')
            Wrtr.writerow([f0.split("/")[0],pnt,str(numbr),str(SmS1+SmG1+SmI1),str(SmG1), str(SmI1), str(SmS1)])

        if (SmG1+SmS1)!=0:
            if 'YOUNG' in f0.upper():
                if 'CONTROL' in f0.upper():
                    ArrYOCTL.append(SmS1/(SmS1+SmG1))
                    ArrIYOCTL.append(SmI1)
                    ArrSYOCTL.append(SmS1)
                    ArrGYOCTL.append(SmG1)

                else:
                    ArrYO.append(SmS1/(SmS1+SmG1))
                    ArrIYO.append(SmI1)
                    ArrSYO.append(SmS1)
                    ArrGYO.append(SmG1)
            elif 'T2D' in f0:
                if 'CONTROL' in f0.upper():
                    ArrT2DCTL.append(SmS1/(SmS1+SmG1))
                    ArrIT2DCTL.append(SmI1)
                    ArrST2DCTL.append(SmS1)
                    ArrGT2DCTL.append(SmG1)
                else:
                    ArrT2D.append(SmS1/(SmS1+SmG1))
                    ArrIT2D.append(SmI1)
                    ArrST2D.append(SmS1)
                    ArrGT2D.append(SmG1)
            else:
                if 'CONTROL' in f0.upper():
                    ArrT1DCTL.append(SmS1/(SmS1+SmG1))
                    ArrIT1DCTL.append(SmI1)
                    ArrST1DCTL.append(SmS1)
                    ArrGT1DCTL.append(SmG1)
                else:
                    ArrT1D.append(SmS1/(SmS0+SmG1))
                    ArrIT1D.append(SmI1)
                    ArrST1D.append(SmS1)
                    ArrGT1D.append(SmG1)



        ##If Silent is off, then plot these cells
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Assigned cells. Alpha cells are red. Beta cells are blue.' + '\n' + r'Delta cells are green. Islet is: %s and patient is: %s' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(Canvas1, interpolation='none')
            plt.savefig(opaths+'/'+numbr+'_Assgned_Cells'+ext, interpolation='none',cmap='gray')
            plt.close()

        fracs=[SmI1,SmS1,SmG1]
        str1, str2, str3=str(SmI1)+' Beta Cells',str(SmS1)+' Delta Cells',str(SmG1)+' Alpha Cells'
        labels=[str1, str2, str3]
        colors=['red','blue','green']
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Pie Chart Displaying the amount of cells.' +'\n' + r'Islet is: %s and patient is: %s' %(numbr,pnt))
            ax=plt.pie(fracs, labels=labels, colors=colors)
            plt.savefig(opaths+'/'+numbr+'IsltPChart'+ext, interpolation='none')
            plt.close()




    if 'T1D' in f0:
        ArrT1DCTL0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0 ])
        ArrIT1DCTL0.append(ArrI1)
        ArrST1DCTL0.append(ArrS1)
        ArrGT1DCTL0.append(ArrG1)
        T1DCTL_sm.append(P_sm)
    elif 'T2D' in f0:
        ArrT2DCTL0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1)if (x+y) !=0 ])
        ArrIT2DCTL0.append(ArrI1)
        ArrST2DCTL0.append(ArrS1)
        ArrGT2DCTL0.append(ArrG1)
        T2DCTL_sm.append(P_sm)
    else:
        ArrYOCTL0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1)if (x+y) !=0 ])
        ArrIYOCTL0.append(ArrI1)
        ArrSYOCTL0.append(ArrS1)
        ArrGYOCTL0.append(ArrG1)
        YOCTL_sm.append(P_sm)


    if 'YOUNG' in f0.upper():
        if 'CONTROL' in f0.upper():
            ArrYOCTL0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0])
            ArrIYOCTL0.append(ArrI1)
            ArrSYOCTL0.append(ArrS1)
            ArrGYOCTL0.append(ArrG1)
            CCount_T1D0.append([x+y+z for x,y,z in zip(ArrS1,ArrG1,ArrI1)])
        else:
            ArrYO0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0])
            ArrIYO0.append(ArrI1)
            ArrSYO0.append(ArrS1)
            ArrGYO0.append(ArrG1)
            CCount_T1DCTL0.append([x+y+z for x,y,z in zip(ArrS1,ArrG1,ArrI1)])
    elif 'T2D' in f0:
        if 'CONTROL' in f0.upper():
            ArrT2DCTL0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0])
            ArrIT2DCTL0.append(ArrI1)
            ArrST2DCTL0.append(ArrS1)
            ArrGT2DCTL0.append(ArrG1)
            CCount_T2D0.append([x+y+z for x,y,z in zip(ArrS1,ArrG1,ArrI1)])
        else:
            ArrT2D0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0])
            ArrIT2D0.append(ArrI1)
            ArrST2D0.append(ArrS1)
            ArrGT2D0.append(ArrG1)
            CCount_T2DCTL0.append([x+y+z for x,y,z in zip(ArrS1,ArrG1,ArrI1)])
    else:
        if 'CONTROL' in f0.upper():
            ArrT1DCTL0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0])
            ArrIT1DCTL0.append(ArrI1)
            ArrST1DCTL0.append(ArrS1)
            ArrGT1DCTL0.append(ArrG1)
            CCount_YO0.append([x+y+z for x,y,z in zip(ArrS1,ArrG1,ArrI1)])
        else:
            ArrT1D0.append([x/(x+y) for x,y in zip(ArrS1,ArrG1) if (x+y) !=0])
            ArrIT1D0.append(ArrI1)
            ArrST1D0.append(ArrS1)
            ArrGT1D0.append(ArrG1)
            CCount_YOCTL0.append([x+y+z for x,y,z in zip(ArrS1,ArrG1,ArrI1)])


    fracs=[SmI0,SmS0,SmG0]
    str1, str2, str3=str(SmI0)+' Beta Cells',str(SmS0)+' Delta Cells',str(SmG0)+' Alpha Cells'
    labels=[str1, str2, str3]
    colors=['red','blue','green']
    if Silent == False:
        f, ax = plt.subplots()
        ax.set_title(r'Pie Chart Displaying the amount of cells in patient %s' %(pnt))
        ax=plt.pie(fracs, labels=labels, colors=colors)
        plt.savefig(opaths+'/CellCountPChart'+ext, interpolation='none')
        plt.close()
    print(len(ArrYOCTL),len(ArrYO),len(ArrT1DCTL),len(ArrT1D),len(ArrT2DCTL),len(ArrT2D))
    print('Done for patient %s, %s of %s' %(f0, str(PNum),str(len(lst0+lst2+lst1)-1) ) )



fracs=[SmIT1CTL,SmST1CTL,SmGT1CTL]
str1, str2, str3=str(SmIT1CTL)+' Beta Cells',str(SmST1CTL)+' Delta Cells',str(SmGT1CTL)+' Alpha Cells'
labels=[str1, str2, str3]
colors=['red','blue','green']
if Silent == False:
    f, ax = plt.subplots()
    ax.set_title(r'Pie Chart Displaying the Amount of Cells in Type 1 Diabetes Controls')
    ax=plt.pie(fracs, labels=labels, colors=colors)
    plt.savefig(path2+'/CellCountPChart_T1D'+ext, interpolation='none')
    plt.close()

fracs=[SmIT2CTL,SmST2CTL,SmGT2CTL]
str1, str2, str3=str(SmIT2CTL)+' Beta Cells',str(SmST2CTL)+' Delta Cells',str(SmGT2CTL)+' Alpha Cells'
labels=[str1, str2, str3]
colors=['red','blue','green']
if Silent == False:
    f, ax = plt.subplots()
    ax.set_title(r'Pie Chart Displaying the Amount of Cells in Type 2 Diabetes Controls')
    ax=plt.pie(fracs, labels=labels, colors=colors)
    plt.savefig(path2+'/CellCountPChart_T2D'+ext, interpolation='none')
    plt.close()

fracs=[SmIYOCTL,SmSYOCTL,SmGYOCTL]
str1, str2, str3=str(SmIYOCTL)+' Beta Cells',str(SmSYOCTL)+' Delta Cells',str(SmGYOCTL)+' Alpha Cells'
labels=[str1, str2, str3]
colors=['red','blue','green']
if Silent == False:
    f, ax = plt.subplots()
    ax.set_title(r'Pie Chart Displaying the Amount of Cells in the Young Onset Control Sample')
    ax=plt.pie(fracs, labels=labels, colors=colors)
    plt.savefig(path2+'/CellCountPChart_T1Control'+ext, interpolation='none')
    plt.close()

K2_T1DCTL, P_T1DCTL=scstats.normaltest(ArrT1DCTL)
K2_T2DCTL, P_T2DCTL=scstats.normaltest(ArrT2DCTL)
K2_YOCTL, P_YOCTL=scstats.normaltest(ArrYOCTL)
K2_T1D, P_T1D=scstats.normaltest(ArrT1D)
K2_T2D, P_T2D=scstats.normaltest(ArrT2D)
K2_YO, P_YO=scstats.normaltest(ArrYO)


K2_GT1DCTL, P_GT1DCTL=scstats.normaltest(ArrGT1DCTL)
K2_GT2DCTL, P_GT2DCTL=scstats.normaltest(ArrGT2DCTL)
K2_GYOCTL, P_GYOCTL=scstats.normaltest(ArrGYOCTL)
K2_GT1D, P_GT1D=scstats.normaltest(ArrGT1D)
K2_GT2D, P_GT2D=scstats.normaltest(ArrGT2D)
K2_GYO, P_GYO=scstats.normaltest(ArrGYO)


K2_IT1D, P_IT1D=scstats.normaltest(ArrIT1D)
K2_IT2D, P_IT2D=scstats.normaltest(ArrIT2D)
K2_IYO, P_IYO=scstats.normaltest(ArrIYO)
K2_IT1DCTL, P_IT1DCTL=scstats.normaltest(ArrIT1DCTL)
K2_IT2DCTL, P_IT2DCTL=scstats.normaltest(ArrIT2DCTL)
K2_IYOCTL, P_IYOCTL=scstats.normaltest(ArrIYOCTL)


K2_ST1D, P_ST1D=scstats.normaltest(ArrST1D)
K2_ST2D, P_ST2D=scstats.normaltest(ArrST2D)
K2_SYO, P_SYO=scstats.normaltest(ArrSYO)
K2_ST1DCTL, P_ST1DCTL=scstats.normaltest(ArrST1DCTL)
K2_ST2DCTL, P_ST2DCTL=scstats.normaltest(ArrST2DCTL)
K2_SYOCTL, P_SYOCTL=scstats.normaltest(ArrSYOCTL)

pathSt=path2+'/StatText.txt'
if os.path.exists(pathSt):
    os.remove(pathSt)
StT=open(pathSt,'w+')
str0='Type 1 Control [#delta]/[#delta+#alpha] cell percentage normalness has a P val of: ' + str(P_T1DCTL)
str2='Young Onset Control [#delta]/[#delta+#alpha] cell percentage normalness has a P val of: ' + str(P_YOCTL)
str4='Type 2 Control [#delta]/[#delta+#alpha] cell percentage normalness has a P val of: ' + str(P_T1DCTL)
str1='Type 1 [#delta]/[#delta+#alpha] cell percentage normalness has a P val of: ' + str(P_T1D)
str3='Young Onset [#delta]/[#delta+#alpha] cell percentage normalness has a P val of: ' + str(P_YO)
str5='Type 2 [#delta]/[#delta+#alpha] cell percentage normalness has a P val of: ' + str(P_T1D)



strG0='Type 1 Control alpha cell count normalness has a P val of: ' + str(P_GT1DCTL)
strG2='Young Onset Control alpha cell count normalness has a P val of: ' + str(P_GYOCTL)
strG4='Type 2 Control alpha cell count normalness has a P val of: ' + str(P_GT2DCTL)
strG1='Type 1 alpha cell count normalness has a P val of: ' + str(P_GT1D)
strG3='Young Onset alpha cell count normalness has a P val of: ' + str(P_GYO)
strG5='Type 2 alpha cell count normalness has a P val of: ' + str(P_GT2D)

strI0='Type 1 Control beta cell count normalness has a P val of: ' + str(P_IT1DCTL)
strI2='Young Onset Control beta cell count normalness has a P val of: ' + str(P_IYOCTL)
strI4='Type 2 Control beta cell count normalness has a P val of: ' + str(P_IT2DCTL)
strI1='Type 1 beta cell count normalness has a P val of: ' + str(P_IT1D)
strI3='Young Onset beta cell count normalness has a P val of: ' + str(P_IYO)
strI5='Type 2 beta cell count normalness has a P val of: ' + str(P_IT2D)


strS0='Type 1 Control delta cell count normalness has a P val of: ' + str(P_ST1DCTL)
strS2='Young Onset Control delta cell count normalness has a P val of: ' + str(P_SYOCTL)
strS4='Type 2 Control delta cell count normalness has a P val of: ' + str(P_ST2DCTL)
strS1='Type 1 Control delta cell count normalness has a P val of: ' + str(P_ST1D)
strS3='Young Onset Control delta cell count normalness has a P val of: ' + str(P_SYO)
strS5='Type 2 Control delta cell count normalness has a P val of: ' + str(P_ST2D)


ST1D, PT1D=scstats.ttest_ind(ArrT1DCTL,ArrT1D, equal_var=False)
ST2D, PT2D=scstats.ttest_ind(ArrT2DCTL,ArrT2D, equal_var=False)
SYO, PYO=scstats.ttest_ind(ArrYOCTL,ArrYO, equal_var=False)

str10='Type 1 [#delta]/[#delta+#alpha] values chance of being different just by chance has a P val of: ' + str(PT1D)
str12='Type 2 [#delta]/[#delta+#alpha] values chance of being different just by chance has a P val of: ' + str(PT2D)
str11='Young onset [#delta]/[#delta+#alpha] values chance of being different just by chance has a P val of: ' + str(PYO)

DT1D, PvalT1D=scstats.ks_2samp(ArrT1DCTL,ArrT1D)
DT2D, PvalT2D=scstats.ks_2samp(ArrT2DCTL,ArrT2D)
DYO, PvalYO=scstats.ks_2samp(ArrYOCTL,ArrYO)

DGT1D, PGvalT1D=scstats.ks_2samp(ArrGT1DCTL,ArrGT1D)
DGT2D, PGvalT2D=scstats.ks_2samp(ArrGT2DCTL,ArrGT2D)
DGYO, PGvalYO=scstats.ks_2samp(ArrGYOCTL,ArrGYO)

DIT1D, PIvalT1D=scstats.ks_2samp(ArrIT1DCTL,ArrIT1D)
DIT2D, PIvalT2D=scstats.ks_2samp(ArrIT2DCTL,ArrIT2D)
DIYO, PIvalYO=scstats.ks_2samp(ArrIYOCTL,ArrIYO)

DST1D, PSvalT1D=scstats.ks_2samp(ArrST1DCTL,ArrST1D)
DST2D, PSvalT2D=scstats.ks_2samp(ArrST2DCTL,ArrST2D)
DSYO, PSvalYO=scstats.ks_2samp(ArrSYOCTL,ArrSYO)

str20='Chance of type 1 [#delta]/[#delta+#alpha] values and type 1 controls being drawn from the same distribution is: ' + str(PvalT1D)
str21='Chance of type 2 [#delta]/[#delta+#alpha] values and type 2 controls being drawn from the same distribution is: ' + str(PvalT2D)
str22='Chance of young onset [#delta]/[#delta+#alpha] values and young onset controls being drawn from the same distribution is: ' + str(PvalYO)


strG10='Chance of type 1 alpha cell count values and type 1 controls being drawn from the same distribution is: ' + str(PGvalT1D)
strG11='Chance of type 2 alpha cell count values and type 2 controls being drawn from the same distribution is: ' + str(PGvalT2D)
strG12='Chance of young onset alpha cell count values and young onset controls being drawn from the same distribution is: ' + str(PGvalYO)

strI10='Chance of type 1 beta cell count values and type 1 controls being drawn from the same distribution is: ' + str(PIvalT1D)
strI11='Chance of type 2 beta cell count values and type 2 controls being drawn from the same distribution is: ' + str(PIvalT2D)
strI12='Chance of young onset beta cell count values and young onset controls being drawn from the same distribution is: ' + str(PIvalYO)

strS10='Chance of type 1 delta cell count values and type 1 controls being drawn from the same distribution is: ' + str(PSvalT1D)
strS11='Chance of type 2 delta cell count values and type 2 controls being drawn from the same distribution is: ' + str(PSvalT2D)
strS12='Chance of young onset delta cell count values and young onset controls being drawn from the same distribution is: ' + str(PSvalYO)

DKT1D, PvalKT1D=scstats.mannwhitneyu(ArrT1DCTL,ArrT1D)
DKT2D, PvalKT2D=scstats.mannwhitneyu(ArrT2DCTL,ArrT2D)
DKYO, PvalKYO=scstats.mannwhitneyu(ArrYOCTL,ArrYO)


DKGT1D, PGvalKT1D=scstats.mannwhitneyu(ArrGT1DCTL,ArrGT1D)
DKGT2D, PGvalKT2D=scstats.mannwhitneyu(ArrGT2DCTL,ArrGT2D)
DKGYO, PGvalKYO=scstats.mannwhitneyu(ArrGYOCTL,ArrGYO)
DKIT1D, PIvalKT1D=scstats.mannwhitneyu(ArrIT1DCTL,ArrIT1D)
DKIT2D, PIvalKT2D=scstats.mannwhitneyu(ArrIT2DCTL,ArrIT2D)
DKIYO, PIvalKYO=scstats.mannwhitneyu(ArrIYOCTL,ArrIYO)
DKST1D, PSvalKT1D=scstats.mannwhitneyu(ArrST1DCTL,ArrST1D)
DKST2D, PSvalKT2D=scstats.mannwhitneyu(ArrST2DCTL,ArrST2D)
DKSYO, PSvalKYO=scstats.mannwhitneyu(ArrSYOCTL,ArrSYO)

str23='Chance of type 1 [#delta]/[#delta+#alpha] values and type 1 controls having the same median value is: ' + str(PvalKT1D)
str24='Chance of type 2 [#delta]/[#delta+#alpha] values and type 2 controls having the same median value is: ' + str(PvalKT2D)
str25='Chance of young onset [#delta]/[#delta+#alpha] values and young onset controls having the same median value is: ' + str(PvalKYO)


strG13='Chance of type 1 alpha cell count values and type 1 controls having the same median value is: ' + str(PGvalKT1D)
strG14='Chance of type 2 alpha cell count values and type 2 controls having the same median value is: ' + str(PGvalKT2D)
strG15='Chance of young onset alpha cell count values and young onset controls having the same median value is: ' + str(PGvalKYO)


strI13='Chance of type 1 beta cell count values and type 1 controls having the same median value is: ' + str(PIvalKT1D)
strI14='Chance of type 2 beta cell count values and type 2 controls having the same median value is: ' + str(PIvalKT2D)
strI15='Chance of young oonset beta cell count values and young onset controls having the same median value is: ' + str(PIvalKYO)



strS13='Chance of type 1 delta cell count values and type 1 controls having the same median value is: ' + str(PSvalKT1D)
strS14='Chance of type 2 delta cell count values and type 2 controls having the same median value is: ' + str(PSvalKT2D)
strS15='Chance of young onset delta cell count values and young onset controls having the same median value is: ' + str(PSvalKYO)

StT.write(str20+'\n'+str21+'\n'+str22+'\n'+strG10+'\n'+strG11+'\n'+strG12+'\n'+strI10+'\n'+strI11+'\n'+strI12+'\n'+strS10+'\n'+strS11+'\n'+strS12+'\n'+'\n' +"Kolomogorov-Smirnov statistic for testing sameness of distribution for the two samples was done using https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ks_2samp.html \n\n\n")
StT.write(str23+'\n'+str24+'\n'+str25+'\n'+strG13+'\n'+strG14+'\n'+strG15+'\n'+ strI13+'\n'+strI14+'\n'+strI15+'\n'+strS13+'\n'+strS14+'\n'+strS15+'\n'+'\n' +"Mann Whitney test for testing sameness of median for the two samples was done using https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.mannwhitneyu.html \n\n\n")

StT.write(str0+'\n' + '     Description of this distribution is: ')
StT.write(str(scstats.describe(ArrT1D)))
StT.write('\n\n' + str1+'\n'+ '     Description of this distribution is: ')
StT.write(str(scstats.describe(ArrT2D)))
StT.write('\n\n' + str2+'\n'+ '     Description of this distribution is: ')
StT.write(str(scstats.describe(ArrT1DCTL)))
StT.write('\n\n' + str3+'\n'+ '     Description of this distribution is: ')
StT.write(str(scstats.describe(ArrT2DCTL)))
StT.write('\n\n' + str4+'\n'+ '     Description of this distribution is: ')
StT.write(str(scstats.describe(ArrYOCTL)))
StT.write('\n\n' + str5+'\n'+ '     Description of this distribution is: ')
StT.write(str(scstats.describe(ArrYO)))


#StT.write('\n\n\n\n'+strG0+'\n' + '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrGT1D)))
#StT.write('\n\n' + strG1+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrGT2D)))
#StT.write('\n\n' + strG2+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrGT1DCTL)))
#StT.write('\n\n' + strG3+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrGT2DCTL)))

#StT.write('\n\n\n\n'+strI0+'\n' + '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrIT1D)))
#StT.write('\n\n' + strI1+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrIT2D)))
#StT.write('\n\n' + strI2+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrIT1DCTL)))
#StT.write('\n\n' + strI3+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrIT2DCTL)))

#StT.write('\n\n\n\n'+strS0+'\n' + '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrST1D)))
#StT.write('\n\n' + strS1+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrST2D))
#StT.write('\n\n' + strS2+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrST1DCTL)))
#StT.write('\n\n' + strS3+'\n'+ '     Description of this distribution is: ')
#StT.write(str(scstats.describe(ArrST2DCTL)))


#StT.write("\n \n \n"+"Normalness test was conducted using https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.normaltest.html#scipy.stats.normaltest \n"+"Data description was done using https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.describe.html#scipy.stats.describe \n \n \n")
#StT.write(str4+'\n'+str5+'\n'+'\n' +"Welch's T test conducted to determine whether results were difference was done using https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.ttest_ind.html \n\n\n")
#StT.write(str6+'\n'+str7+'\n'+strG6+'\n'+strG7+'\n'+strI6+'\n'+strI7+'\n'+strS6+'\n'+strS7+'\n'+'\n' +"Kolomogorov-Smirnov statistic for testing sameness of distribution for the two samples was done using https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.stats.ks_2samp.html \n\n\n")
#StT.write(str8+'\n'+str9+'\n'+strG8+'\n'+strG9+'\n'+ strI8+'\n'+strI9+'\n'+strS8+'\n'+strS9+'\n'+'\n' +"Kruskal-Wallis test for testing sameness of median for the two samples was done using https://docs.scipy.org/doc/scipy-0.19.1/reference/generated/scipy.stats.kruskal.html")

if Silent == True:
    f, (ax1,ax2,ax3) = plt.subplots(3)
    ax1.set_title(r'Histogram of the [#delta]/[#delta+#alpha] cells for Young onset control')
    ax2.set_xlabel(r'Percent of delta cells')
    ax3.set_xlabel(r'Percent of delta cells')
    ax1.set_xlabel(r'Percent of delta cells')
    ax1.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax1.hist(ArrYOCTL, normed=True, bins=20)
    normedmean(ArrYOCTL,ax1)
    ax2.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax2.set_title(r'Histogram of the [#delta]/[#delta+#alpha] cells for T1D Controls')
    ax2.hist(ArrT1DCTL, normed=True, bins=20)
    normedmean(ArrT1DCTL,ax2)
    ax3.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax2.set_title(r'Histogram of the [#delta]/[#delta+#alpha] cells for T2D Controls')
    ax2.hist(ArrT2DCTL, normed=True, bins=20)
    normedmean(ArrT2DCTL,ax3)

    plt.tight_layout()
    plt.savefig(path2+'/'+'Cell_ratio_histogram'+ext, interpolation='none')
    plt.close()


if Silent == True:
    f, (ax1,ax2,ax3)=plt.subplots(3)

    ax1.set_title(r'Histogram of the [#alpha] cells for young onset controls')
    ax1.set_xlabel(r'numbers of alpha cells per islet')
    ax1.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax1.hist(ArrGYOCTL, normed=True, bins=20)
    normedmean(ArrGYOCTL,ax1)

    ax2.set_title(r'Histogram of the [#alpha] cells for T1D controls')
    ax2.set_xlabel(r'numbers of alpha cells per islet')
    ax2.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax2.hist(ArrGT1DCTL, normed=True, bins=20)
    normedmean(ArrGT1DCTL,ax2)

    ax3.set_title(r'Histogram of the [#alpha] cells for T2D controls')
    ax3.set_xlabel(r'numbers of alpha cells per islet')
    ax3.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax3.hist(ArrGT2D, normed=True, bins=20)
    normedmean(ArrGT2D,ax3)
    plt.tight_layout()
    plt.savefig(path2+'/'+'Control_alpha_count_histogram'+ext, interpolation='none')
    plt.close()

if Silent == True:
    f, (ax1,ax2,ax3)=plt.subplots(3)

    ax1.set_title(r'Histogram of the [#beta] cells for young onset controls')
    ax1.set_xlabel(r'numbers of beta cells per islet')
    ax1.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax1.hist(ArrIYOCTL, normed=True, bins=20)
    normedmean(ArrIYOCTL,ax1)

    ax2.set_title(r'Histogram of the [#beta] cells for T1D controls')
    ax2.set_xlabel(r'numbers of beta cells per islet')
    ax2.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax2.hist(ArrIT1DCTL, normed=True, bins=20)
    normedmean(ArrIT1DCTL,ax2)

    ax3.set_title(r'Histogram of the [#beta] cells for T2D controls')
    ax3.set_xlabel(r'numbers of beta cells per islet')
    ax3.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax3.hist(ArrIT2D, normed=True, bins=20)
    normedmean(ArrIT2D,ax3)
    plt.tight_layout()
    plt.savefig(path2+'/'+'Control_beta_count_histogram'+ext, interpolation='none')
    plt.close()


if Silent == True:
    f, (ax1,ax2,ax3)=plt.subplots(3)

    ax1.set_title(r'Histogram of the [#delta] cells for young onset controls')
    ax1.set_xlabel(r'numbers of delta cells per islet')
    ax1.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax1.hist(ArrSYOCTL, normed=True, bins=20)
    normedmean(ArrSYOCTL,ax1)

    ax2.set_title(r'Histogram of the [#delta] cells for T1D controls')
    ax2.set_xlabel(r'numbers of delta cells per islet')
    ax2.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax2.hist(ArrST1DCTL, normed=True, bins=20)
    normedmean(ArrST1DCTL,ax2)

    ax3.set_title(r'Histogram of the [#delta] cells for T2D controls')
    ax3.set_xlabel(r'numbers of delta cells per islet')
    ax3.set_ylabel(r'Normed frequency count', rotation='vertical')
    ax3.hist(ArrST2D, normed=True, bins=20)
    normedmean(ArrST2D,ax3)
    plt.tight_layout()
    plt.savefig(path2+'/'+'Control_delta_count_histogram'+ext, interpolation='none')
    plt.close()






fracs=[SmIT2CTL,SmST2CTL,SmGT2CTL]
str1, str2, str3=str(SmIT2CTL)+' Beta Cells',str(SmST2CTL)+' Delta Cells',str(SmGT2CTL)+' Alpha Cells'
labels=[str1, str2, str3]
colors=['red','blue','green']
if Silent == False:
    f, ax = plt.subplots()
    ax.set_title(r'Pie Chart Displaying the Amount of Cells in the Type 2 Control Sample')
    ax=plt.pie(fracs, labels=labels, colors=colors)
    plt.savefig(path2+'/CellCountPChart_T2Control'+ext, interpolation='none')
    plt.close()

fracs=[SmIT1CTL,SmST1CTL,SmGT1CTL]
str1, str2, str3=str(SmIT1CTL)+' Beta Cells',str(SmST1CTL)+' Delta Cells',str(SmGT1CTL)+' Alpha Cells'
labels=[str1, str2, str3]
colors=['red','blue','green']
if Silent == False:
    f, ax = plt.subplots()
    ax.set_title(r'Pie Chart Displaying the Amount of Cells in the Type 1 Control Sample')
    ax=plt.pie(fracs, labels=labels, colors=colors)
    plt.savefig(path2+'/CellCountPChart_T1Control'+ext, interpolation='none')
    plt.close()

fracs=[SmIYOCTL,SmSYOCTL,SmGYOCTL]
str1, str2, str3=str(SmIYOCTL)+' Beta Cells',str(SmSYOCTL)+' Delta Cells',str(SmGYOCTL)+' Alpha Cells'
labels=[str1, str2, str3]
colors=['red','blue','green']
if Silent == False:
    f, ax = plt.subplots()
    ax.set_title(r'Pie Chart Displaying the Amount of Cells in the Young onset Control Sample')
    ax=plt.pie(fracs, labels=labels, colors=colors)
    plt.savefig(path2+'/CellCountPChart_T1Control'+ext, interpolation='none')
    plt.close()

if Silent== False:
    f, ax = plt.subplots()
    ax.set_title(r'Bar chart showing the cell count between the 3 cell types')
    N=3
    ind = np.arange(N)  # the x locations for the groups
    width = 0.15       # the width of the bars
    T1DCTL_means=[np.median(ArrGT1DCTL),np.median(ArrIT1DCTL),np.median(ArrST1DCTL)]
    T1DCTL_err=[[np.percentile(ArrGT1DCTL,25),np.percentile(ArrIT1DCTL,25),np.percentile(ArrST1DCTL,25)],[np.percentile(ArrGT1DCTL,75),np.percentile(ArrIT1DCTL,75),np.percentile(ArrST1DCTL,75)]]
    T2DCTL_means=[np.median(ArrGT2DCTL),np.median(ArrIT2DCTL),np.median(ArrST2DCTL)]
    T2DCTL_err=[[np.percentile(ArrGT2DCTL,25),np.percentile(ArrIT2DCTL,25),np.percentile(ArrST2DCTL,25)],[np.percentile(ArrGT2DCTL,75),np.percentile(ArrIT2DCTL,75),np.percentile(ArrST2DCTL,75)]]
    YOCTL_means=[np.median(ArrGYOCTL),np.median(ArrIYOCTL),np.median(ArrSYOCTL)]
    YOCTL_err=[[np.percentile(ArrGYOCTL,25),np.percentile(ArrIYOCTL,25),np.percentile(ArrSYOCTL,25)],[np.percentile(ArrGYOCTL,75),np.percentile(ArrIYOCTL,75),np.percentile(ArrSYOCTL,75)]]

    T1D_means2=[np.mean(ArrGT1D),np.mean(ArrIT1D),np.mean(ArrST1D)]
    T1D_err2=[np.std(ArrGT1D),np.std(ArrIT1D),np.std(ArrST1D)]
    T1DCTL_means2=[np.mean(ArrGT1DCTL),np.mean(ArrIT1DCTL),np.mean(ArrST1DCTL)]
    T1DCTL_err2=[np.std(ArrGT1DCTL),np.std(ArrIT1DCTL),np.std(ArrST1DCTL)]
    T2D_means2=[np.mean(ArrGT2D),np.mean(ArrIT2D),np.mean(ArrST2D)]
    T2D_err2=[np.std(ArrGT2D),np.std(ArrIT2D),np.std(ArrST2D)]
    T2DCTL_means2=[np.mean(ArrGT2DCTL),np.mean(ArrIT2DCTL),np.mean(ArrST2DCTL)]
    T2DCTL_err2=[np.std(ArrGT2DCTL),np.std(ArrIT2DCTL),np.std(ArrST2DCTL)]
    YO_means2=[np.mean(ArrGYO),np.mean(ArrIYO),np.mean(ArrSYO)]
    YO_err2=[np.std(ArrGYO),np.std(ArrIYO),np.std(ArrSYO)]
    YOCTL_means2=[np.mean(ArrGYOCTL),np.mean(ArrIYOCTL),np.mean(ArrSYOCTL)]
    YOCTL_err2=[np.std(ArrGYOCTL),np.std(ArrIYOCTL),np.std(ArrSYOCTL)]



    ax.set_xlabel(r'Cell type')
    ax.set_ylabel(r'Cell count', rotation='vertical')
    ax.set_xticks(ind + 2.5*width)
    ax.set_xticklabels(('Alpha Cells', 'Beta Cells', 'Delta Cells'))
    rects1=ax.bar(ind, T1D_means2, width, color='r', yerr=T1D_err2, ecolor='0.7', capsize=5)
    rects2=ax.bar(ind+width, T1DCTL_means2, width, color='c', yerr=T1DCTL_err2, ecolor='0.7', capsize=5)
    rects3=ax.bar(ind+2*width, T2D_means2, width, color='g', yerr=T2D_err2, ecolor='0.7', capsize=5)
    rects4=ax.bar(ind+3*width, T2DCTL_means2, width, color='m',  yerr=T2DCTL_err2, ecolor='0.7', capsize=5)
    rects5=ax.bar(ind+4*width, YO_means2, width, color='b', yerr=YO_err2, ecolor='0.7', capsize=5)
    rects6=ax.bar(ind+5*width, YOCTL_means2, width, color='y',  yerr=YOCTL_err2, ecolor='0.7', capsize=5)
    ax.legend((rects1[0], rects2[0], rects3[0], rects4[0], rects5[0], rects6[0]), ('Old Onset', 'Old Control','T2D', 'T2D Control','Young Onset', 'Young Control'))
    def round_sig(x, sig=3):
        return round(x, sig-int(floor(log10(abs(x))))-1)
    def errbardraw(arr1,arr2,xpos):
        for i in range(N):
            x1=xpos+i
            ax.plot([x1, x1], [arr1[i],arr2[i]], lw=1.5, c='0.7')

    #pertl=np.percentile(T1,75)
    #pertlCTL=np.percentile(T1CTL,75)
    #pertl=np.mean(T1)+np.std(T1)
    #pertlT1CTL=np.mean(T1CTL)+np.std(T1CTL)
    #pertlT2CTL=np.mean(T2CTL)+np.std(T1CTL)
    #y, h, col = max(pertl,pertlCTL), 4, 'k'
    #x1,x2=nn+0.0,nn+0.15
    #ax.plot([x1, x1, x2, x2], [pertl+2, y+h, y+h, pertlCTL+2], lw=1.5, c=col)
    #stat=scstats.mannwhitneyu(T1,T1CTL)[1]
    #if stat <0.001:
    #    statstr='***'
    #elif stat <0.01:
    #    statstr='**'
    #elif stat <0.05:
    #    statstr='*'
    #else:
    #    statstr=''
    #ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' +str(round_sig(stat)), ha='center', va='bottom', color=col)

    for nn, (i, j,T1,T1CTL) in enumerate(zip (rects1,rects2,[ArrGT1D,ArrIT1D,ArrST1D],[ArrGT1DCTL,ArrIT1DCTL,ArrST1DCTL])):
        #pertl=np.percentile(T1,75)
        #pertlCTL=np.percentile(T1CTL,75)
        pertl=np.mean(T1)+np.std(T1)
        pertlCTL=np.mean(T1CTL)+np.std(T1CTL)
        y, h, col = max(pertl,pertlCTL), 4, 'k'
        x1,x2=nn+0.0,nn+0.15
        ax.plot([x1, x1, x2, x2], [pertl+2, y+h, y+h, pertlCTL+2], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(T1,T1CTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' +str(round_sig(stat)), ha='center', va='bottom', color=col, fontsize=6)

    for nn, (i, j,T2,T2CTL) in enumerate(zip (rects3,rects4,[ArrGT2D,ArrIT2D,ArrST2D],[ArrGT2DCTL,ArrIT2DCTL,ArrST2DCTL])):
        #pertl=np.percentile(T2,75)
        #pertlCTL=np.percentile(T2CTL,75)
        pertl=np.mean(T2)+np.std(T2)
        pertlCTL=np.mean(T2CTL)+np.std(T2CTL)
        y, h, col = max(pertl,pertlCTL), 4, 'k'
        x1,x2=nn+0.3,nn+0.45
        ax.plot([x1, x1, x2, x2], [pertl+2, y+h, y+h, pertlCTL+2], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(T2,T2CTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' + str(round_sig(stat)), ha='center', va='bottom', color=col, fontsize=6)

    for nn, (i, j,YO,YOCTL) in enumerate(zip (rects5,rects6,[ArrGYO,ArrIYO,ArrSYO],[ArrGYOCTL,ArrIYOCTL,ArrSYOCTL])):
        #pertl=np.percentile(T2,75)
        #pertlCTL=np.percentile(T2CTL,75)
        pertl=np.mean(YO)+np.std(YO)
        pertlCTL=np.mean(YOCTL)+np.std(YOCTL)
        y, h, col = max(pertl,pertlCTL), 4, 'k'
        x1,x2=nn+0.6,nn+0.75
        ax.plot([x1, x1, x2, x2], [pertl+2, y+h, y+h, pertlCTL+2], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(YO,YOCTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' + str(round_sig(stat)), ha='center', va='bottom', color=col, fontsize=6)




    ax.set_ylim([0,180])
    plt.savefig(path2+'/Count_barchart_Celltype'+ext, interpolation='none')
    plt.close()


if Silent == True:
    f, ax = plt.subplots()
    N=3
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars




    flierprops = dict(marker='.', markerfacecolor='k', markersize=0)
    rects1=ax.boxplot(positions=ind, x=[ArrGT1D,ArrIT1D, ArrST1D], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects2=ax.boxplot(positions=ind+width, x=[ArrGT1DCTL,ArrIT1DCTL, ArrST1DCTL], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects3=ax.boxplot(positions=ind+2*width, x=[ArrGT2D,ArrIT2D, ArrST2D], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects4=ax.boxplot(positions=ind+3*width, x=[ArrGT2DCTL,ArrIT2DCTL, ArrST2DCTL], widths=width, vert=True, patch_artist=True, flierprops=flierprops)


    #sns.swarmplot(data=[np.zeros_like(ArrGT1D),ArrGT1D], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrIT1D), ArrIT1D], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrST1D)*2, ArrST1D], size=2, ax=ax)

    #sns.swarmplot(data=[np.zeros_like(ArrGT1DCTL)+width, ArrGT1DCTL], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrIT1DCTL)+width, ArrIT1DCTL], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrST1DCTL)*2+width, ArrST1DCTL], size=2, ax=ax)

    #sns.swarmplot(data=[np.zeros_like(ArrGT2D)+width*2, ArrGT2D], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrIT2D)+width*2, ArrIT2D], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrST2D)*2+width*2, ArrST2D], size=2, ax=ax)

    #sns.swarmplot(data=[np.zeros_like(ArrGT2DCTL)+width*3, ArrGT2DCTL], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrIT2DCTL)+width*3, ArrIT2DCTL], size=2, ax=ax)
    #sns.swarmplot(data=[np.ones_like(ArrST2DCTL)*2+width*3, ArrST2DCTL], size=2, ax=ax)
    #ax.legend((rects1['boxes'][0], rects2['boxes'][0], rects3['boxes'][0], rects4['boxes'][0]), ('T1D', 'T1Dcontrol','T2D', 'T2Dcontrol'))
    #ax.set_ylim([0,155])
    ax.set_title(r'Bar chart showing the cell count between the 3 cell types')
    ax.set_xlabel(r'Cell type')
    ax.set_ylabel(r'Cell count', rotation='vertical')
    ax.set_xticks(ind + 1.5*width)
    ax.set_xticklabels(('Alpha Cells', 'Beta Cells', 'Delta Cells'))


    xc=-0.2
    for num, i in enumerate([ArrGT1D,ArrGT1DCTL,ArrGT2D,ArrGT2DCTL,ArrIT1D,ArrIT1DCTL,ArrIT2D,ArrIT2DCTL,ArrST1D,ArrST1DCTL,ArrST1D,ArrST1DCTL]):
        #print(i)
        #print('\n')
        xc+=0.2
        if num % 4 == 0 and num != 0:
            xc+=0.2
        xx=[]
        for j in i:
            rand=rd.uniform(-0.03,0.03)
            xx.append(xc+rand)
        xx=np.array((xx),dtype=float)
        ii=np.array((i),dtype=float)
        #print(xx)
        #print(ii)
        #print('\n')
        #sns.swarmplot(data=[xx,ii], size=2, ax=ax)

        #print('bar number is'+str(num) + ' and xdist is ' + str(xc) + '\n')
        ax.scatter(x=xx,y=ii, s=2.5,c='k',linewidths=0, zorder=10)

        for i in rects1['boxes']:
            i.set_facecolor('c')
        for i in rects2['boxes']:
            i.set_facecolor('m')
        for i in rects3['boxes']:
            i.set_facecolor('y')
        for i in rects4['boxes']:
            i.set_facecolor('g')

    ax.legend((rects1['boxes'][0], rects2['boxes'][0], rects3['boxes'][0], rects4['boxes'][0]), ('T1D', 'T1Dcontrol','T2D', 'T2Dcontrol'))

    for nn, (T1,T1CTL) in enumerate(zip ([ArrGT1D,ArrIT1D,ArrST1D],[ArrGT1DCTL,ArrIT1DCTL,ArrST1DCTL])):
        pertl=np.max(T1)
        pertlCTL=np.max(T1CTL)
        y, h, col = max(pertl,pertlCTL), 10, 'k'
        x1,x2=nn+0.0,nn+0.2
        ax.plot([x1, x1, x2, x2], [pertl+5, y+h, y+h, pertlCTL+5], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(T1,T1CTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' +str(round_sig(stat)), ha='center', va='bottom', color=col)

    for nn, (T2,T2CTL) in enumerate(zip ([ArrGT2D,ArrIT2D,ArrST2D],[ArrGT2DCTL,ArrIT2DCTL,ArrST2DCTL])):
        pertl=np.max(T2)
        pertlCTL=np.max(T2CTL)
        y, h, col = max(pertl,pertlCTL), 10, 'k'
        x1,x2=nn+0.4,nn+0.6
        ax.plot([x1, x1, x2, x2], [pertl+5, y+h, y+h, pertlCTL+5], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(T2,T2CTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' + str(round_sig(stat)), ha='center', va='bottom', color=col)
    ax.set_xlim([-0.2,2.8])
    ax.set_ylim([0,275])
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='on') # labels along the bottom edge are on
    plt.savefig(path2+'/Count_BoxWhisker_Celltype'+ext, interpolation='none')
    plt.close()




CCount_T1D=np.asarray(ArrGT1D)+np.asarray(ArrIT1D)+np.asarray(ArrST1D)
CCount_T1DCTL=np.asarray(ArrGT1DCTL)+np.asarray(ArrIT1DCTL)+np.asarray(ArrST1DCTL)
CCount_T2D=np.asarray(ArrGT2D)+np.asarray(ArrIT2D)+np.asarray(ArrST2D)
CCount_T2DCTL=np.asarray(ArrGT2DCTL)+np.asarray(ArrIT2DCTL)+np.asarray(ArrST2DCTL)

CCount_YO=np.asarray(ArrGYO)+np.asarray(ArrIYO)+np.asarray(ArrSYO)
CCount_YOCTL=np.asarray(ArrGYOCTL)+np.asarray(ArrIYOCTL)+np.asarray(ArrSYOCTL)








if Silent == False:
    f, ax = plt.subplots()
    N=1
    ind = np.arange(N)  # the x locations for the groups
    width = 0.2       # the width of the bars
    flierprops = dict(marker='.', markerfacecolor='k', markersize=0)
    rects1=ax.boxplot(positions=ind, x=[CCount_T1D], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects2=ax.boxplot(positions=ind+width, x=[CCount_T1DCTL], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects3=ax.boxplot(positions=ind+2*width, x=[CCount_T2D], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects4=ax.boxplot(positions=ind+3*width, x=[CCount_T2DCTL], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects5=ax.boxplot(positions=ind+4*width, x=[CCount_YO], widths=width, vert=True, patch_artist=True, flierprops=flierprops)
    rects6=ax.boxplot(positions=ind+5*width, x=[CCount_YOCTL], widths=width, vert=True, patch_artist=True, flierprops=flierprops)


    for i in rects1['boxes']:
        i.set_facecolor('r')
    for i in rects2['boxes']:
        i.set_facecolor('c')
    for i in rects3['boxes']:
        i.set_facecolor('g')
    for i in rects4['boxes']:
        i.set_facecolor('m')
    for i in rects5['boxes']:
        i.set_facecolor('b')
    for i in rects6['boxes']:
        i.set_facecolor('y')

    ax.legend((rects1['boxes'][0], rects2['boxes'][0], rects3['boxes'][0], rects4['boxes'][0], rects5['boxes'][0], rects6['boxes'][0]), ('Old Onset', 'Old Control','T2D', 'T2D Control','Young Onset', 'Young Control'),ncol=3)
    #ax.set_ylim([0,155])
    ax.set_title(r'Box whisker plot showing the total number of cells per islet')
    ax.set_xlabel(r'Patient type')
    ax.set_ylabel(r'Cell count', rotation='vertical')

    for nn, (T1,T1CTL) in enumerate(zip ([CCount_T1D],[CCount_T1DCTL])):
        pertl=np.max(T1)
        pertlCTL=np.max(T1CTL)
        y, h, col = max(pertl,pertlCTL), 10, 'k'
        x1,x2=nn+0.0,nn+0.2
        ax.plot([x1, x1, x2, x2], [pertl+5, y+h, y+h, pertlCTL+5], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(T1,T1CTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' +str(round_sig(stat)), ha='center', va='bottom', color=col)

    for nn, (T2,T2CTL) in enumerate(zip ([CCount_T2D],[CCount_T2DCTL])):
        pertl=np.max(T2)
        pertlCTL=np.max(T2CTL)
        y, h, col = max(pertl,pertlCTL), 10, 'k'
        x1,x2=nn+0.4,nn+0.6
        ax.plot([x1, x1, x2, x2], [pertl+5, y+h, y+h, pertlCTL+5], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(T2,T2CTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' + str(round_sig(stat)), ha='center', va='bottom', color=col)

    for nn, (YO,YOCTL) in enumerate(zip ([CCount_YO],[CCount_YOCTL])):
        pertl=np.max(YO)
        pertlCTL=np.max(YOCTL)
        y, h, col = max(pertl,pertlCTL), 10, 'k'
        x1,x2=nn+0.8,nn+1.0
        ax.plot([x1, x1, x2, x2], [pertl+5, y+h, y+h, pertlCTL+5], lw=1.5, c=col)
        stat=scstats.mannwhitneyu(YO,YOCTL, alternative='two-sided')[1]
        if stat <0.001:
            statstr='***'
        elif stat <0.01:
            statstr='**'
        elif stat <0.05:
            statstr='*'
        else:
            statstr=''
        ax.text((x1+x2)*.5, y+h, statstr+'P_val = \n' + str(round_sig(stat)), ha='center', va='bottom', color=col)



    xc=-0.2
    for num, i in enumerate([CCount_T1D,CCount_T1DCTL,CCount_T2D,CCount_T2DCTL,CCount_YO,CCount_YOCTL]):
        xc+=0.2
        xx=[]
        for j in i:
            rand=rd.uniform(-0.03,0.03)
            xx.append(xc+rand)
        xx=np.array((xx),dtype=float)
        ii=np.array((i),dtype=float)
        ax.scatter(x=xx,y=ii, s=2.5,c='k',linewidths=0, zorder=10)
    ax.set_xlim([-0.2,1.2])
    ax.set_ylim([0,400])
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom='off',      # ticks along the bottom edge are off
        top='off',         # ticks along the top edge are off
        labelbottom='off') # labels along the bottom edge are on
    ax.set_xticklabels("")
    plt.savefig(path2+'/Overall_Count_BoxWhisker'+ext, interpolation='none')
    plt.close()


pathSt=path2+'/Data_desc.csv'
if os.path.exists(pathSt):
    os.remove(pathSt)
with open(pathSt,mode='w+',newline='') as csvfile:
    Wrtr=csv.writer(csvfile, delimiter=',')


#    Wrtr.writerow(['T1D_subject_islet']+['CellCount']+['']+['']+['T1DCTL_subject_islet']+['CellCount']+['']+['']+['T2D_subject_islet']+['CellCount']+['']+['']+['T2DCTL_subject_islet']+['CellCount'])


#    for pnum, (T1n,T1CTLn,T2n,T2CTLn) in enumerate(itertools.zip_longest(T1D_sm,T1DCTL_sm,T2D_sm,T2DCTL_sm, fillvalue=[])):

#        for inum,(T1nn,T1CTLnn,T2nn,T2CTLnn) in enumerate(itertools.zip_longest(T1n,T1CTLn,T2n,T2CTLn)):
#            if T1nn is None:
#                str1=''
#                str10=''
#            else:
#                str1= 'T1D_S%d_I%d' %(pnum+1,inum+1)
#                str10= str(T1nn)
#            if T1nn is None:
#                str2=''
#                str20=''
#            else:
#                str2= 'T1DCTL_S%d_I%d' %(pnum+1,inum+1)
#                str20= str(T1CTLnn)
#            if T2nn is None:
#                str3=''
#                str30=''
#            else:
#                str3= 'T2D_S%d_I%d' %(pnum+1,inum+1)
#                str30= str(T2nn)
#            if T1nn is None:
#                str4=''
#                str40=''
#            else:
#                str4= 'T2DCTL_S%d_I%d' %(pnum+1,inum+1)
#                str40= str(T2CTLnn)
#            Wrtr.writerow([str1]+[str10]+['']+['']+[str2]+[str20]+['']+['']+[str3]+[str30]+['']+['']+[str4]+[str40])

    for i in range(5):
        Wrtr.writerow([''])
    Wrtr.writerow(['Dataset_name','Dataset_Size','Data_Mean', 'Data_Variance', 'Data_Std_Dev', 'Data_Skewness', 'Data_Kurtotis', 'Data_normalness_test_stat', 'Data_normalness_probability','K-W interpersonal test for subject species test stat','K-W interpersonal test for subject species P_Val'])
    for i,j,nam in zip([ArrIT1D,ArrST1D,ArrGT1D,CCount_T1D,ArrIT1DCTL,ArrST1DCTL,ArrGT1DCTL, CCount_T1DCTL, ArrIT2D,ArrST2D,ArrGT2D, CCount_T2D, ArrIT2DCTL,ArrST2DCTL,ArrGT2DCTL,CCount_T2DCTL,ArrIYO,ArrSYO,ArrGYO,CCount_YO,ArrIYOCTL,ArrSYOCTL,ArrGYOCTL,CCount_YOCTL],[ArrIT1D0,ArrST1D0,ArrGT1D0,CCount_T1D0,ArrIT1DCTL0,ArrST1DCTL0,ArrGT1DCTL0, CCount_T1DCTL0, ArrIT2D0,ArrST2D0,ArrGT2D0, CCount_T2D0, ArrIT2DCTL0, ArrST2DCTL0, ArrGT2DCTL0, CCount_T2DCTL0, ArrIYO0, ArrSYO0, ArrGYO0, CCount_YO0, ArrIYOCTL0, ArrSYOCTL0, ArrGYOCTL0, CCount_YO0],['T1D patient [#beta] cells','T1D patient [#delta] cells','T1D patient [#alpha] cells','T1D patient total cells','T1D contol [#beta] cells','T1D control [#delta] cells','T1D control [#alpha] cells','T1D control total cells','T2D patient [#beta] cells','T2D patient [#delta] cells','T2D patient [#alpha] cells','T2D patient total cells','T2D control [#beta] cells', 'T2D control [#delta] cells','T2D control [#alpha] cells','T2D control total cells', 'Young onset patient [#beta] cells','Young onset [#delta] cells','Young onset [#alpha] cells','Young onset total cells','Young onset control [#beta] cells', 'Young onset control [#delta] cells','Young control [#alpha] cells','Young control total cells']):
        #print([nam]+[scstats.describe(i)[0]]+[scstats.describe(i)[2]]+[scstats.describe(i)[3]]+[math.sqrt(scstats.describe(i)[3])]+[scstats.describe(i)[4]]+[scstats.describe(i)[5]]+[scstats.normaltest(i)[0]]+[scstats.normaltest(i)[1]])
        #print(j)
        #print([scstats.kruskal(*j)[0]]+[scstats.kruskal(*j)[1]])

        Wrtr.writerow([nam]+[scstats.describe(i)[0]]+[scstats.describe(i)[2]]+[scstats.describe(i)[3]]+[math.sqrt(scstats.describe(i)[3])]+[scstats.describe(i)[4]]+[scstats.describe(i)[5]]+[scstats.normaltest(i)[0]]+[scstats.normaltest(i)[1]]+[scstats.kruskal(*j)[0]]+[scstats.kruskal(*j)[1]])


    for i in range(5):
        Wrtr.writerow([''])
    Wrtr.writerow(['Intersubject comparison'])
    Wrtr.writerow(['Dataset_1','Dataset_2','KS test stat','KS test P_val', 'MW U test stat', 'MW U test P_val','MW Effect size', 'MW inverse effect size'])
    for i,j,nami,namj in zip([ArrIT1D,ArrST1D,ArrGT1D,CCount_T1D,ArrIT2D,ArrST2D,ArrGT2D,CCount_T2D,ArrIYO,ArrSYO,ArrGYO,CCount_YO],[ArrIT1DCTL,ArrST1DCTL,ArrGT1DCTL,CCount_T1DCTL,ArrIT2DCTL,ArrST2DCTL,ArrGT2DCTL,CCount_T2DCTL,ArrIYOCTL,ArrSYOCTL,ArrGYOCTL,CCount_YOCTL],['T1D patient [#beta] cells','T1D patient [#delta] cells','T1D patient [#alpha] cells','T1D patient total cells','T2D patient [#beta] cells','T2D patient [#delta] cells','T2D patient [#alpha] cells','T2D patient total cells','YO patient [#beta] cells','YO patient [#delta] cells','YO patient [#alpha] cells','YO patient total cells'],['T1D contol [#beta] cells','T1D control [#delta] cells','T1D control [#alpha] cells','T1D control total cells','T2D control [#beta] cells', 'T2D control [#delta] cells','T2D control [#alpha] cells','T2D control total cells','YO control [#beta] cells', 'YO control [#delta] cells','YO control [#alpha] cells','YO control total cells']):
        Wrtr.writerow([nami]+[namj]+[scstats.ks_2samp(i,j)[0]]+[scstats.ks_2samp(i,j)[1]]+[scstats.mannwhitneyu(i,j, alternative='two-sided')[0]]+[scstats.mannwhitneyu(i,j, alternative='two-sided')[1]]+[scstats.mannwhitneyu(i,j, alternative='two-sided')[0]/(len(i)*len(j))]+[scstats.mannwhitneyu(j,i, alternative='two-sided')[0]/(len(i)*len(j))])
print("done")
