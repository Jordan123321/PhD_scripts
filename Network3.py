#!/usr/bin/env python3

##This is a flag.
Silent=False

##These are library imports. OS is for file manipulation
import os
##Numpy is for mathematical manipulation, especially with matrices
import numpy as np

import pandas as pd


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
import scipy.spatial as scspat
import scipy.ndimage as ndim

##Scikit image is a package which contains a large number of image processing functions
import skimage.io as skio
import skimage.morphology as skmorph
import skimage.filters as filt
import skimage.measure as skmeas
import skimage.segmentation as skseg
import skimage.draw as skdr
#import skimage.color as skcol

#import csv
import math
#import random as rd

##Imports my Voronoi Split Algorithm into a Module
import VorSplit
import networkx as nx
import sys
Rat=25/(92*1000000)


def normedmean(arr,ax):
    mean = np.mean(arr)
    variance = np.var(arr)
    sigma = np.sqrt(variance)
    x = np.linspace(min(arr), max(arr), 100)
    ax.plot(x, mlab.normpdf(x, mean, sigma))


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
            j.set_xlabel(name+r' within the islet (' +units+')')
        else:
            j.set_xlabel(name+r' within the islet')
        j.set_xlim(xmin,xmax)
        j.set_ylim(0,ymax)
        #j.set_xticklabels(["{:.2e}".format(t) for t in j.get_xticks()])
    #plt.tight_layout()
    plt.savefig(path2+'/'+name+ext, interpolation='none')
    plt.close()



##Defines paths to my directories and save file locations
path1=os.getcwd()
path0=os.path.dirname(path1)
path2=path0+ '/Output_Network_3'
ext='.svg'



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
    opathsGph=opaths+'/Graph'
    if not os.path.exists(opathsGph):
        os.makedirs(opathsGph)

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

        for islt_n,islt in enumerate(skmeas.regionprops(IsltShpB1)):
            lcount=0
            islt_ns=str(islt_n)
            ##Mask the nuclei with the islet shape
            IsltNuc=np.multiply(NucPos[islt.bbox[0]:islt.bbox[2],islt.bbox[1]:islt.bbox[3]],islt.image)
            IsltShpB1_1, B1_1 = skmorph.label(islt.image,return_num=True)
            B2Regs_1=np.logical_not(islt.image)
            B2Regs2_1=skseg.clear_border(B2Regs_1)
            IsltShpB2_1, B2_1 = skmorph.label(B2Regs2_1,return_num=True)
            ##Fill the holes to make contingent nuclei
            NucPos2=ndim.binary_fill_holes(IsltNuc)
            ##Remove small artefacts
            NucPos3=skmorph.remove_small_holes(NucPos2, connectivity=1, area_threshold=1000)
            NucPos4=skmorph.remove_small_objects(NucPos3, connectivity=1, min_size=100)
    
    
            ##Erode the image, to disentangle joint up DAPI stains
            NucPos5=skmorph.erosion(NucPos4, selem=skmorph.disk(9))
    
    
    
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
            NPosdf=pd.DataFrame({'x':NPosArr[:,0],'y':NPosArr[:,1]})
            NPosdf['c']='k'
            NPosdf['Connected']=0
            NPosdf['To_Alpha']=0
            NPosdf['To_Beta']=0
            NPosdf['To_Delta']=0
            ##Skip the rest of this loop in the unlikely event that there are no nuclei
            if (NucNum <4):
                continue
            ##Do a voronoi split on the Islet and the points of the nuclei centres
            Canvas= VorSplit.VorSplt(NPosArr, islt.image.copy())
    
            tri=scspat.Delaunay(NPosArr)
    
            ##Labels the parts of the voronoi algorithm for the proporties array. 1 connectivity denotes the single line separating the cells
            CellLab,CellNum = skmeas.label(Canvas, return_num=1,connectivity=1,background=0)
            ##Get properties of the regions
            props=skmeas.regionprops(CellLab)
            ##Flatten the scan arrays
            IimF=Iim[:,:,0]+Iim[:,:,1]+Iim[:,:,2]
            GimF=Gim[:,:,0]+Gim[:,:,1]+Gim[:,:,2]
            SimF=Sim[:,:,0]+Sim[:,:,1]+Sim[:,:,2]
            for row in NPosdf.itertuples():
                NoNuc=True
                for n,p in enumerate(props):
    
                    if ( (math.floor(row[1])in p.coords[:,0] or math.ceil(row[1])in p.coords[:,0])
                        and (math.floor(row[2]) in p.coords[:,1] or math.ceil(row[2]) in p.coords[:,1])):
                        NoNuc=False
                        cell_x,cell_y=row[1],row[2]
                        tmp_ind=row[0]
                        break
                if NoNuc:
                    print('huh')
                    continue
    
                ##Mask the insulin with the cells as defined by the Voronoi algorithm
                xmin,xmax,ymin,ymax=islt.bbox[0]+p.bbox[0],islt.bbox[0]+p.bbox[2],islt.bbox[1]+p.bbox[1],islt.bbox[1]+p.bbox[3]
                CellI=np.sum(np.multiply(p.image,IimF[xmin:xmax,ymin:ymax]))
                CellG=np.sum(np.multiply(p.image,GimF[xmin:xmax,ymin:ymax]))
                CellS=np.sum(np.multiply(p.image,SimF[xmin:xmax,ymin:ymax]))
    
    
    
                if (CellI>CellG and CellI>CellS):
                    #Insulin (beta) blue
                    NPosdf.loc[tmp_ind,'c']='b'
                elif CellG>CellS:
                    #Glucagon (alpha) red
                    NPosdf.loc[tmp_ind,'c']='r'
                else:
                    #Somatostatin (delta) green
                    NPosdf.loc[tmp_ind,'c']='g'
    
    
            lines=[]
            ##Set area threshold for enclosed space to work out whether connected or not
            Thresh=200
            Lenthresh=150
            for smpnum, smpi in enumerate(tri.simplices):
                p0, p1, p2 = tri.points[smpi].astype(int)
    
                dist0=scspat.distance.euclidean(p0,p1)
                if dist0<=Lenthresh:
                    triangleimg=B2Regs_1.copy()
                    rr,cc=skdr.line(*p0,*p1)
                    triangleimg[rr,cc]=0
                    labelimg=skmeas.label(triangleimg,connectivity=1)
                    skseg.clear_border(labelimg, in_place=True)
                    
                    
                    areas=[]
                    for smlreg in skmeas.regionprops(labelimg):
                        areas.append(smlreg.area)
                    for _ in range(B2_1):
                        areas.remove(max(areas))
                    if (sum(areas)<=Thresh and (IsltShpB1_1[p0[0],p0[1]]-IsltShpB1_1[p1[0],p1[1]])==0):
                        a=[smpi[0],smpi[1]]
                        if sorted(a) not in lines:
                            lines.append(sorted(a))
    
                dist0=scspat.distance.euclidean(p1,p2)
                if dist0<=Lenthresh:
                    triangleimg=B2Regs_1.copy()
                    rr,cc=skdr.line(*p1,*p2)
                    triangleimg[rr,cc]=0
                    #f, ax = plt.subplots()
                    #ax.set_title('test1')
                    #ax.imshow(triangleimg, interpolation='none')
                    #f.tight_layout()
                    #lcounts=str(lcount)
                    #plt.savefig(opathsGph+'/'+numbr+'_'+islt_ns+'_'+lcounts+'_test1'+ext, interpolation='none')
                    #plt.close()
                    
                    labelimg=skmeas.label(triangleimg,connectivity=1)
                    skseg.clear_border(labelimg,in_place=True)
                    
                    
                    #f, ax = plt.subplots()
                    #ax.set_title('test2')
                    #ax.imshow(labelimg, interpolation='none')
                    #f.tight_layout()
                    #plt.savefig(opathsGph+'/'+numbr+'_'+islt_ns+'_'+lcounts+'_test2'+ext, interpolation='none')
                    #plt.close()
                    
                    #lcount+=1
                    
                    
                    areas=[]
                    for smlreg in skmeas.regionprops(labelimg):
                        areas.append(smlreg.area)
                    for _ in range(B2_1):
                        areas.remove(max(areas))
                    if (sum(areas)<=Thresh and (IsltShpB1_1[p1[0],p1[1]]-IsltShpB1_1[p2[0],p2[1]])==0):
                        a=[smpi[1],smpi[2]]
                        if sorted(a) not in lines:
                            lines.append(sorted(a))
    
    
    
                dist0=scspat.distance.euclidean(p2,p0)
                if dist0<=Lenthresh:
                    triangleimg=B2Regs_1.copy()
                    rr,cc=skdr.line(*p2,*p0)
                    triangleimg[rr,cc]=0
                    labelimg=skmeas.label(triangleimg,connectivity=1)
                    skseg.clear_border(labelimg,in_place=True)
                    
                    areas=[]
                    for smlreg in skmeas.regionprops(labelimg):
                        areas.append(smlreg.area)
                    for _ in range(B2_1):
                        areas.remove(max(areas))
                    if (sum(areas)<=Thresh and (IsltShpB1_1[p2[0],p2[1]]-IsltShpB1_1[p0[0],p0[1]])==0):
                        a=[smpi[2],smpi[0]]
                        if sorted(a) not in lines:
                            lines.append(sorted(a))
    
    
            lines2=pd.DataFrame(columns=['x0','y0','x1','y1','c','p0','p1'])
            SmRR=0
            SmRB=0
            SmRG=0
            SmGG=0
            SmGB=0
            SmBB=0
            if len(lines)==0:
                continue
    
            for ln in lines:
                NPosdf.iloc[ln[0],3]+=1
                NPosdf.iloc[ln[1],3]+=1
                P_a=NPosdf.loc[ln[0]]
                P_b=NPosdf.loc[ln[1]]
                bob=P_a[2]
                if (P_a[2] == 'r' and P_b[2] == 'r'):
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'r',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmRR+=1
                    NPosdf.iloc[ln[0],4]+=1
                    NPosdf.iloc[ln[1],4]+=1
    
                elif (P_a[2] == 'r' and P_b[2] == 'g'):
                    NPosdf.iloc[ln[0],6]+=1
                    NPosdf.iloc[ln[1],4]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'y',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmRG+=1
    
                elif (P_a[2] == 'g' and P_b[2] == 'r'):
                    NPosdf.iloc[ln[0],4]+=1
                    NPosdf.iloc[ln[1],6]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'y',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmRG+=1
    
                elif (P_a[2] == 'r' and P_b[2] == 'b'):
                    NPosdf.iloc[ln[0],5]+=1
                    NPosdf.iloc[ln[1],4]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'m',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmRB+=1
    
                elif (P_a[2] == 'b' and P_b[2] == 'r'):
                    NPosdf.iloc[ln[0],4]+=1
                    NPosdf.iloc[ln[1],5]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'m',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmRB+=1
    
                elif (P_a[2] == 'g' and P_b[2] == 'g'):
                    NPosdf.iloc[ln[0],6]+=1
                    NPosdf.iloc[ln[1],6]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'g',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmGG+=1
    
                elif (P_a[2] == 'b' and P_b[2] == 'g'):
                    NPosdf.iloc[ln[0],6]+=1
                    NPosdf.iloc[ln[1],5]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'c',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmGB+=1
    
                elif (P_a[2] == 'g' and P_b[2] == 'b'):
                    NPosdf.iloc[ln[0],5]+=1
                    NPosdf.iloc[ln[1],6]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'c',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmGB+=1
    
                elif (P_a[2] == 'b' and P_b[2] == 'b'):
                    NPosdf.iloc[ln[0],5]+=1
                    NPosdf.iloc[ln[1],5]+=1
                    lines2=lines2.append(pd.Series([P_a[0],P_a[1],P_b[0],P_b[1],'b',ln[0],ln[1]],index=['x0','y0','x1','y1','c','p0','p1']), ignore_index=True)
                    SmBB+=1
    
                else:
                    print("Huh")
                    
            if Silent == False:
                f, ax = plt.subplots()
                ax.set_title('Delanay of Nuclei. Islet is: \n %s and patient is %s \n' %(numbr,pnt))
                ax.set_xlabel(r'X direction (pixels) $\rightarrow$'+'\n\n\n')
                ax.set_ylabel(r'$\leftarrow$ Y direction (pixels)', rotation='vertical')
                ax.imshow(Canvas, interpolation='none',cmap='gray')
                for ll in lines2.itertuples():
                    ax.plot([ll[2],ll[4]],[ll[1],ll[3]],'-',color=ll[5],linewidth=0.5)
                #ax.triplot(NPosArr[:,1], NPosArr[:,0], tri.simplices)
                xcoord,ycoord=NPosdf['y'].tolist(),NPosdf['x'].tolist()
                ax.scatter(x=xcoord, y=ycoord, marker='o',color=NPosdf['c'].tolist(),s=2)
                #for i, txt in enumerate(NPosdf['Connected'].tolist()):
                    #ax.annotate(str(txt), (xcoord[i], ycoord[i]),fontsize=4, bbox=dict(facecolor='white',lw=0,boxstyle="round", alpha=0.5))
                ax.set_xlim([0,NucPos4.shape[1]])
                ax.set_ylim([0,NucPos4.shape[0]])
                #plt.gca().invert_yaxis()
                f.text(0.5,.05, 'First Betti Number is :%i \n Second Betti Number is: %i' %(B1,B2),ha='right')
                f.tight_layout()
                plt.savefig(opathsGph+'/'+numbr+'_'+islt_ns+'_Spat_Netw'+ext, interpolation='none')
                plt.close()
    
            G=nx.Graph()
            for row in NPosdf.iterrows():
                key=row[0]
                tmpdata=row[1]
                G.add_node(key,pos=(tmpdata[0],tmpdata[1]),c=tmpdata[2])
            for row in lines2.iterrows():
                key=row[0]
                tmpdata=row[1]
                G.add_edge(tmpdata[5],tmpdata[6], c=tmpdata[4],pos1=(tmpdata[0],tmpdata[1]),pos2=(tmpdata[2],tmpdata[3]))
            
            posdict=NPosdf[['y','x']].T.to_dict('list')
            sp_pos=nx.spring_layout(G,pos=posdict, k=100)
            ci_pos=nx.circular_layout(G)
            pl_pos=nx.planar_layout(G)
            kk_pos=nx.kamada_kawai_layout(G,pos=posdict)
            spec_pos=nx.spectral_layout(G)
            ecol=[G[f][t]['c'] for f,t in G.edges()]
            ncol=NPosdf['c'].tolist()
            
            
            f, ax = plt.subplots()
            ax.set_title('Network representation of islet. Layout is truepos \n Islet is: %s, Image is: %s and patient is %s' %(islt_ns,numbr,pnt))
            ax.axis('off')
            nx.draw_networkx(G,ax=ax,pos=posdict,with_labels=False,node_color=ncol,edge_color=ecol,node_size=50,width=1.5)
            plt.savefig(opathsGph+'/'+numbr+'_'+islt_ns+'_Original_Network'+ext, interpolation='none')
            plt.tight_layout()
            plt.close()
            
            f, ax = plt.subplots()
            ax.set_title('Network representation of islet. Layout is spring \n Islet is: %s, Image is: %s and patient is %s' %(islt_ns,numbr,pnt))
            ax.axis('off')
            nx.draw_networkx(G,ax=ax,pos=sp_pos,with_labels=False,node_color=ncol,edge_color=ecol,node_size=50,width=1.5)
            plt.savefig(opathsGph+'/'+numbr+'_'+islt_ns+'_Spring_Network'+ext, interpolation='none')
            plt.tight_layout()
            plt.close()
            
            f, ax = plt.subplots()
            ax.set_title('Network representation of islet. Layout is circular \n Islet is: %s, Image is: %s and patient is %s' %(islt_ns,numbr,pnt))
            ax.axis('off')
            nx.draw_networkx(G,ax=ax,pos=ci_pos,with_labels=False,node_color=ncol,edge_color=ecol,node_size=50,width=1.5)
            plt.savefig(opathsGph+'/'+numbr+'_'+islt_ns+'_Circular_Network'+ext, interpolation='none')
            plt.tight_layout()
            plt.close()
            


