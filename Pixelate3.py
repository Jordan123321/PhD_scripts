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
##Imports latex labels
from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

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
import skimage.transform as sktrans



import csv
import math
import random as rd

##Imports my Voronoi Split Algorithm into a Module
import VorSplit
from math import log10, floor

import seaborn as sns

import itertools

##Sets the x and y pixel numbers for the mesh
xRes=100
yRes=100


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
path2=path0+ '/Output_Pixel_3'
ext='.eps'



##Creates list of paths to diabetes patients
lst0= ["T1D/"+f for f in os.listdir(path0+"/T1D")if (f.startswith('T'))]
lst1= ["T2D/"+f for f in os.listdir(path0+"/T2D")if (f.startswith('T'))]
lst2= ["Young_Onset/" + f for f in os.listdir(path0+"/Young_Onset" ) if ('CONTROL' in f.upper()) or ('CASE') in f.upper()]




##Zeroes all the sums of the types of patients (Sm)=sum     (I)=Beta    (T1D)=Type 1 Diabetes       (CTL)=Control
##                                                          (S)=Delta   (T2D)=Type 2 Diabetes
##                                                          (G)=Alpha   (YO)=Young Onset





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

    DAPIlst.sort(key=lambda x: x.split()[-1])
    GCGlst.sort(key=lambda x: x.split()[-1])
    INSlst.sort(key=lambda x: x.split()[-1])
    SSTlst.sort(key=lambda x: x.split()[-1])
    OALLlst.sort(key=lambda x: x.split()[-1])

    ##This loops through the Scans for the various patients
    for num,(D,G,I,S,O) in enumerate(zip(DAPIlst,GCGlst, INSlst, SSTlst,OALLlst)):
        ##Set the scan Beta, Delta and Alpha cell count to zero, and make an image number string
        SmI1=0
        SmS1=0
        SmG1=0
        tst=D.split()[-1]
        if ((D.split()[-1] != G.split()[-1]) or (D.split()[-1] != I.split()[-1]) or (D.split()[-1] != I.split()[-1]) or (D.split()[-1] != O.split()[-1]) ):
            print(D.split()[-1] ,G.split()[-1] ,I.split()[-1] ,S.split()[-1] ,O.split()[-1] )
            continue
        nmbr=tst[:-4]
        Dim=skio.imread(ppaths+'/'+D)
        Gim=skio.imread(ppaths+'/'+G)
        Iim=skio.imread(ppaths+'/'+I)
        Sim=skio.imread(ppaths+'/'+S)
        Oim=skio.imread(ppaths+'/'+O)
        #print(Dim.shape)
        ##If silet is off, then plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Original overlay. Islet is: %s and patient is %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(Oim, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'_Orig_Img'+ext, interpolation='none')
            plt.close()



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

        ##If silet is off, then plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Islet Shape. Islet is: %s and patient is %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(DiffSim3, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'_BloodFilt'+ext, interpolation='none')
            plt.close()


        ##First thing is to find the islet shape. This can be done by adding the non Dapi, smoothing and thresholding
        ##Add together the three scans
        IsltShp=Gim+Iim+Sim
        ##Add together the red blue and green to flatten the array
        IsltShp0=IsltShp[:,:,0]+IsltShp[:,:,1]+IsltShp[:,:,2]
        ##Get rid of scale bar
        IsltShp0[-50:,-200:]=0


        ##Gaussian smooth followed by a triangle filter to determine the boundary of the islets
        IsltShp1=filt.gaussian(IsltShp0,10)
        IsltShp2=IsltShp1>filt.threshold_triangle(IsltShp1)
        ##Get rid of any small objects that may yet exist
        IsltShp3=skmorph.remove_small_holes(IsltShp2, connectivity=1, min_size=1000)
        IsltShp4=skmorph.remove_small_objects(IsltShp3, connectivity=1, min_size=10000)




        ##If silent is off, then plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Islet Shape. Islet is: %s and patient is %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(IsltShp4, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'_IsltShape'+ext, interpolation='none')
            plt.close()

        ##Next is to find the Nuclei positions
        ##Flatten the arrays
        Nuc=Dim[:,:,0]+Dim[:,:,1]+Dim[:,:,2]

        ##Get rid of scale bar
        Nuc[-50:,-200:]=0

        ##Local threshold
        NucPos=(Nuc>filt.threshold_local(Nuc, block_size=101))
        ##Fill the holes to make contingent nuclei
        IsltNuc=np.multiply(NucPos,IsltShp4)
        NucPos2=ndim.binary_fill_holes(IsltNuc)
        ##Mask the nuclei with the islet shape

        ##Remove small artefacts
        NucPos3=skmorph.remove_small_holes(NucPos2, connectivity=1, min_size=1000)
        NucPos4=skmorph.remove_small_objects(NucPos3, connectivity=1, min_size=100)
        ##Erode the image, to disentangle joint up DAPI stains
        NucPos5=skmorph.erosion(NucPos4, selem=skmorph.disk(7))
        ##If silent is off, plot this out
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Nucleus position. Islet is: %s and patient is %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(NucPos5, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'_Nuclei'+ext, interpolation='none')
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
        if NucNum==0:
            continue
        ##Do a voronoi split on the Islet and the points of the nuclei centres
        canvas= VorSplit.VorSplt(NPosArr, IsltShp4.copy())

        ##Plot this if silent is off
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Cell Shape. Islet is: %s and patient is %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(canvas, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'Cells'+ext, interpolation='none',cmap='gray')
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
        test_Canv=np.ones((Oim.shape[0],Oim.shape[1],3),dtype=np.uint8)
        Canvas2=np.zeros((Oim.shape[0],Oim.shape[1],3),dtype=np.uint8)

        IPosArr=[]
        GPosArr=[]
        SPosArr=[]
        for p in props:
            ##Mask the insulin with the cells as defined by the Voronoi algorithm
            CellI=np.multiply(p.image,IimF[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3]])
            CellG=np.multiply(p.image,GimF[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3]])
            CellS=np.multiply(p.image,SimF[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3]])
            tmpArr=np.round(np.array([p.centroid[0]*xRes/Oim.shape[0],p.centroid[1]*yRes/Oim.shape[1]])).astype(np.uint8)
            ##Calculate the total intensity of the scans within the defined cells
            totI=np.sum(CellI)
            totG=np.sum(CellG)
            totS=3*np.sum(CellS)
            ##Compare these cells and evaluate accordingly
            ##Red is insulin
            if totI>totG and totI>totS:
                Canvas1[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3],2]+=p.image
                SmI0+=1
                SmI1+=1
                IPosArr.append(tmpArr)
            ##Red is glucagon
            elif totG>totS:
                Canvas1[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3],0]+=p.image
                SmG0+=1
                SmG1+=1
                GPosArr.append(tmpArr)
            ##Green is delta cells
            else:
                Canvas1[p.bbox[0]:p.bbox[2],p.bbox[1]:p.bbox[3],1]+=p.image
                SmS0+=1
                SmS1+=1
                SPosArr.append(tmpArr)
        ArrG1.append(SmG1)
        ArrI1.append(SmI1)
        ArrS1.append(SmS1)

        ##If Silent is off, then plot these cells
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Assigned cells. Alpha cells are red. Beta cells are blue.' + '\n' + r'Delta cells are green. Islet is: %s and patient is: %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(Canvas1*255, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'_Assgned_Cells'+ext, interpolation='none',cmap='gray')
            plt.close()

        fracs=[SmI1,SmS1,SmG1]
        str1, str2, str3=str(SmI1)+' Beta Cells',str(SmS1)+' Delta Cells',str(SmG1)+' Alpha Cells'
        labels=[str1, str2, str3]
        colors=['red','blue','green']
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Pie Chart Displaying the amount of cells.' +'\n' + r'Islet is: %s and patient is: %s' %(nmbr,pnt))
            ax=plt.pie(fracs, labels=labels, colors=colors)
            plt.savefig(opaths+'/'+nmbr+'_IsltPChart'+ext, interpolation='none')
            plt.close()
        ##Reshape This Array to xres and yres
        Canvas1=sktrans.resize(Canvas1,(xRes,yRes),order=0,mode="edge").astype(np.bool)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title(r'Pixelated cells. Alpha cells are red. Beta cells are blue.' + '\n' + r'Delta cells are green. Islet is: %s and patient is: %s' %(nmbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$ ')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax=plt.imshow(Canvas1, interpolation='none')
            plt.savefig(opaths+'/'+nmbr+'_Pixelated_Cells'+ext, interpolation='none',cmap='gray')
            plt.close()
        ##Need to find nearest neigbor pixels
        GCList=[['g'] for x in range(len(GPosArr)) ]
        SCList=[['s'] for x in range(len(SPosArr)) ]
        ICList=[['i'] for x in range(len(IPosArr)) ]
        for xx in range(xRes):
            for yy in range(yRes):
                #print(Canvas1[xx,yy])
                if (Canvas1[xx,yy]==[0,0,0]).all():
                    continue
                elif (Canvas1[xx,yy]==[1,0,0]).all():
                    mindist= xRes*yRes
                    for n,G in enumerate(GPosArr):
                        tmpdist = math.sqrt((xx-G[0])*(xx-G[0])+(yy-G[1])*(yy-G[1]))
                        if tmpdist<mindist:
                            mindist=tmpdist
                            ind=n
                    GCList[ind].append(xRes*yy+xx)
                elif (Canvas1[xx,yy]==[0,1,0]).all():
                    mindist= xRes*yRes
                    for n,S in enumerate(SPosArr):
                        tmpdist = math.sqrt((xx-S[0])*(xx-S[0])+(yy-S[1])*(yy-S[1]))
                        if tmpdist<mindist:
                            mindist=tmpdist
                            ind=n
                    SCList[ind].append(xRes*yy+xx)
                elif (Canvas1[xx,yy]==[0,0,1]).all():
                    mindist= xRes*yRes
                    for n,I in enumerate(IPosArr):
                        tmpdist = math.sqrt((xx-I[0])*(xx-I[0])+(yy-I[1])*(yy-I[1]))
                        if tmpdist<mindist:
                            mindist=tmpdist
                            ind=n
                    ICList[ind].append(xRes*yy+xx)
                else:
                    print("Error. Value of array at point is " + str(Canvas1[xx,yy]))

        for cenum, ce in enumerate(GCList):
            GCList[cenum][1:]=sorted(ce[1:])
        for cenum, ce in enumerate(SCList):
            SCList[cenum][1:]=sorted(ce[1:])
        for cenum, ce in enumerate(ICList):
            ICList[cenum][1:]=sorted(ce[1:])

        thefile = open(opaths+'/'+nmbr+'_Cells.txt', 'w+')
        First=True
        for item in GCList:
            if len(item)<2:
                continue
            if First:
                First=False
            else:
                thefile.write("\n")
            for char in item:
                thefile.write("%s " % char)

        for item in SCList:
            if len(item)<2:
                continue
            if First:
                First=False
            else:
                thefile.write("\n")
            for char in item:
                thefile.write("%s " % char)

        for item in ICList:
            if len(item)<2:
                continue
            if First:
                First=False
            else:
                thefile.write("\n")
            for char in item:
                thefile.write("%s " % char)
        thefile.close()
    print('Done for patient %s, %s of %s' %(f0, str(PNum),str(len(lst0+lst2+lst1)-1) ) )

print("done")
