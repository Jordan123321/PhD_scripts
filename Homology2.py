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
#import skimage.color as skcol

#import csv
import math
#import random as rd

##Imports my Voronoi Split Algorithm into a Module
import VorSplit
from math import log10, floor

import seaborn as sns

import itertools

import sys

Rat=25/(92*1000000)


def normedmean(arr,ax):
    mean = np.mean(arr)
    variance = np.var(arr)
    sigma = np.sqrt(variance)
    x = np.linspace(min(arr), max(arr), 100)
    ax.plot(x, mlab.normpdf(x, mean, sigma))




##Defines paths to my directories and save file locations
path1=os.getcwd()
path0=os.path.dirname(path1)
path2=path0+ '/Output_Homology'
ext='.eps'



##Creates list of paths to diabetes patients
lst0= ["T1D/"+f for f in os.listdir(path0+"/T1D")if (f.startswith('T'))]
lst1= ["T2D/"+f for f in os.listdir(path0+"/T2D")if (f.startswith('T'))]
lst2= ["Young_Onset/" + f for f in os.listdir(path0+"/Young_Onset" ) if ('CONTROL' in f.upper()) or ('CASE') in f.upper()]


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
        ##Get rid of the backround blood
        Sim1=Sim[:,:,0]+Sim[:,:,1]+Sim[:,:,2]

        ##Take the Laplacian of the Stomatostatin
        ##Get rid of the scale bar
        DiffSim1=filt.laplace(Sim1)
        ##If silet is off, then plot this out
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

        IsltShp2=IsltShp1>filt.threshold_triangle(IsltShp1)

        ##Get rid of any small objects that may yet exist
        IsltShp3=skmorph.remove_small_holes(IsltShp2, connectivity=1, area_threshold=1000)
        IsltShp4=skmorph.remove_small_objects(IsltShp3, connectivity=1, min_size=10000)
        ##Find the first and second betti Numbers for this image
        IsltShpB1, B1 = skmorph.label(IsltShp4,return_num=True)
        B2Regs=np.logical_not(IsltShp4)
        B2Regs2=skseg.clear_border(B2Regs)
        IsltShpB2, B2 = skmorph.label(B2Regs2,return_num=True)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title('IsletShape. Islet is:\n %s and patient is %s \n' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$'+'\n\n\n')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(IsltShp4, interpolation='none')
            f.text(0.5,.05, 'First Betti Number is :%i \n Second Betti Number is: %i' %(B1,B2),ha='right')
            f.tight_layout()
            plt.savefig(opaths+'/'+numbr+'_Betti'+ext, interpolation='none')
            plt.close()

        ##Next is to find the Nuclei positions
        ##Flatten the arrays
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
        NucPos3=skmorph.remove_small_holes(NucPos2, connectivity=1, min_size=1000)
        NucPos4=skmorph.remove_small_objects(NucPos3, connectivity=1, min_size=100)


        ##Erode the image, to disentangle joint up DAPI stains
        NucPos5=skmorph.erosion(NucPos4, selem=skmorph.disk(7))



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
        Canvas= VorSplit.VorSplt(NPosArr, IsltShp4.copy())
        
        
        tri=scspat.Delaunay(NPosArr)
        lines=[]
        #smpfordel=[]
        for smpnum, smpi in enumerate(tri.simplices):
            p0, p1, p2 = tri.points[smpi].astype(int)
            triangleimg=np.zeros_like(B2Regs,dtype=bool)
            rr,cc=skdr.line(*p0,*p1)
            triangleimg[rr,cc]=1
            if ((not np.logical_and(triangleimg,B2Regs).any()) and (IsltShpB1[p0[0],p0[1]]-IsltShpB1[p1[0],p1[1]])==0):
                a=[smpi[0],smpi[1]]
                lines.append(sorted(a))
            triangleimg=np.zeros_like(B2Regs,dtype=bool)
            rr,cc=skdr.line(*p1,*p2)
            triangleimg[rr,cc]=1
            if ((not np.logical_and(triangleimg,B2Regs).any()) and (IsltShpB1[p1[0],p1[1]]-IsltShpB1[p2[0],p2[1]])==0):
                a=[smpi[1],smpi[2]]
                lines.append(sorted(a))
            triangleimg=np.zeros_like(B2Regs,dtype=bool)
            rr,cc=skdr.line(*p2,*p0)
            triangleimg[rr,cc]=1
            if ((not np.logical_and(triangleimg,B2Regs).any()) and (IsltShpB1[p2[0],p2[1]]-IsltShpB1[p0[0],p0[1]])==0):
                a=[smpi[2],smpi[0]]
                lines.append(sorted(a))
        lines2=[]
        for uniq in lines:
            if uniq not in lines2:
                lines2.append(uniq)
        opathsGph=opaths+'/Graph'
        if not os.path.exists(opathsGph):
            os.makedirs(opathsGph)
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title('Delanay of Nuclei. Islet is: \n %s and patient is %s \n' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$'+'\n\n\n')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(Canvas, interpolation='none')
            for line in lines2:
                ax.plot([NPosArr[line[0]][1],NPosArr[line[1]][1]],[NPosArr[line[0]][0],NPosArr[line[1]][0]],'-',linewidth=0.5)
            #ax.triplot(NPosArr[:,1], NPosArr[:,0], tri.simplices)
            ax.plot(NPosArr[:,1], NPosArr[:,0], 'o', markersize=1)
            ax.set_xlim([0,NucPos4.shape[1]])
            ax.set_ylim([0,NucPos4.shape[0]])
            #plt.gca().invert_yaxis()
            f.text(0.5,.05, 'First Betti Number is :%i \n Second Betti Number is: %i' %(B1,B2),ha='right')
            f.tight_layout()
            plt.savefig(opathsGph+'/'+numbr+'_Strongly_Connected'+ext, interpolation='none')
            plt.close()
        
        
        lines=[]
        for smpnum, smpi in enumerate(tri.simplices):
            p0, p1, p2 = tri.points[smpi].astype(int)
            triangleimg=np.zeros_like(B2Regs2,dtype=bool)
            rr,cc=skdr.line(*p0,*p1)
            triangleimg[rr,cc]=1
            if ((not np.logical_and(triangleimg,B2Regs2).any()) and (IsltShpB1[p0[0],p0[1]]-IsltShpB1[p1[0],p1[1]])==0):
                a=[smpi[0],smpi[1]]
                lines.append(sorted(a))
            triangleimg=np.zeros_like(B2Regs2,dtype=bool)
            rr,cc=skdr.line(*p1,*p2)
            triangleimg[rr,cc]=1
            if ((not np.logical_and(triangleimg,B2Regs2).any()) and (IsltShpB1[p1[0],p1[1]]-IsltShpB1[p2[0],p2[1]])==0):
                a=[smpi[1],smpi[2]]
                lines.append(sorted(a))
            triangleimg=np.zeros_like(B2Regs2,dtype=bool)
            rr,cc=skdr.line(*p2,*p0)
            triangleimg[rr,cc]=1
            if ((not np.logical_and(triangleimg,B2Regs2).any()) and (IsltShpB1[p2[0],p2[1]]-IsltShpB1[p0[0],p0[1]])==0):
                a=[smpi[2],smpi[0]]
                lines.append(sorted(a))
        lines2=[]
        for uniq in lines:
            if uniq not in lines2:
                lines2.append(uniq)
        
        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title('Delanay of Nuclei. Islet is: \n %s and patient is %s \n' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$'+'\n\n\n')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(Canvas, interpolation='none')
            for line in lines2:
                ax.plot([NPosArr[line[0]][1],NPosArr[line[1]][1]],[NPosArr[line[0]][0],NPosArr[line[1]][0]],'-',linewidth=0.5)
            #ax.triplot(NPosArr[:,1], NPosArr[:,0], tri.simplices)
            ax.plot(NPosArr[:,1], NPosArr[:,0], 'o', markersize=1)
            ax.set_xlim([0,NucPos4.shape[1]])
            ax.set_ylim([0,NucPos4.shape[0]])
            #plt.gca().invert_yaxis()
            f.text(0.5,.05, 'First Betti Number is :%i \n Second Betti Number is: %i' %(B1,B2),ha='right')
            f.tight_layout()
            plt.savefig(opathsGph+'/'+numbr+'_Weakly_Connected'+ext, interpolation='none')
            plt.close()
        
        
        lines=[]
        Thresh=200
        #print(tri.simplices)
        for smpnum, smpi in enumerate(tri.simplices):
            p0, p1, p2 = tri.points[smpi].astype(int)
            
            
            triangleimg=B2Regs.copy()
            rr,cc=skdr.line(*p0,*p1)
            triangleimg[rr,cc]=0
            labelimg=skmeas.label(triangleimg,connectivity=1)
            areas=[]
            for smlreg in skmeas.regionprops(labelimg):
                areas.append(smlreg.area)
            
            for _ in range(B2+1):
                areas.remove(max(areas))
            


            if (sum(areas)<=Thresh and (IsltShpB1[p0[0],p0[1]]-IsltShpB1[p1[0],p1[1]])==0):
                a=[smpi[0],smpi[1]]
                lines.append(sorted(a))


            triangleimg=B2Regs.copy()
            rr,cc=skdr.line(*p1,*p2)
            triangleimg[rr,cc]=0
            labelimg=skmeas.label(triangleimg,connectivity=1)
            areas=[]
            for smlreg in skmeas.regionprops(labelimg):
                areas.append(smlreg.area)
            for _ in range(B2+1):
                areas.remove(max(areas))
            if (sum(areas)<=Thresh and (IsltShpB1[p1[0],p1[1]]-IsltShpB1[p2[0],p2[1]])==0):
                a=[smpi[1],smpi[2]]
                lines.append(sorted(a))


            triangleimg=B2Regs.copy()
            rr,cc=skdr.line(*p2,*p0)
            triangleimg[rr,cc]=0
            labelimg=skmeas.label(triangleimg,connectivity=1)
            areas=[]
            for smlreg in skmeas.regionprops(labelimg):
                areas.append(smlreg.area)
            for _ in range(B2+1):
                areas.remove(max(areas))
            if (sum(areas)<=Thresh and (IsltShpB1[p2[0],p2[1]]-IsltShpB1[p0[0],p0[1]])==0):
                a=[smpi[2],smpi[0]]
                lines.append(sorted(a))


        lines2=[]
        for uniq in lines:
            if uniq not in lines2:
                lines2.append(uniq)

        if Silent == False:
            f, ax = plt.subplots()
            ax.set_title('Delanay of Nuclei. Islet is: \n %s and patient is %s \n' %(numbr,pnt))
            ax.set_xlabel(r'X direction (pixels) $\rightarrow$'+'\n\n\n')
            ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
            ax.imshow(Canvas, interpolation='none')
            for line in lines2:
                ax.plot([NPosArr[line[0]][1],NPosArr[line[1]][1]],[NPosArr[line[0]][0],NPosArr[line[1]][0]],'-',linewidth=0.5)
            #ax.triplot(NPosArr[:,1], NPosArr[:,0], tri.simplices)
            ax.plot(NPosArr[:,1], NPosArr[:,0], 'o', markersize=1)
            ax.set_xlim([0,NucPos4.shape[1]])
            ax.set_ylim([0,NucPos4.shape[0]])
            #plt.gca().invert_yaxis()
            f.text(0.5,.05, 'First Betti Number is :%i \n Second Betti Number is: %i' %(B1,B2),ha='right')
            f.tight_layout()
            plt.savefig(opathsGph+'/'+numbr+'_Estimated_Connected'+ext, interpolation='none')
            plt.close()
        #sys.exit()