#!/usr/bin/env python3

##Scipy contains some elements used in image processing that are required here
import scipy.spatial as scspat
import numpy as np
import skimage.draw as skdr

def VorSplt(points, canv):
    canvas=canv
    vor=scspat.Voronoi(points)
    ##The following segment of code arranges the points from the voronoi algorithm so they can easily plotted
    ##Adapted from https://docs.scipy.org/doc/scipy-0.18.1/reference/generated/scipy.spatial.voronoi_plot_2d.html
    center = vor.points.mean(axis=0)
    ptp_bound = vor.points.ptp(axis=0)
    finite_segments = []
    infinite_segments = []
    for pointidx, simplex in zip(vor.ridge_points, vor.ridge_vertices):
        simplex = np.asarray(simplex)
        if np.all(simplex >= 0):
            finite_segments.append(vor.vertices[simplex])
        else:
            ii = simplex[simplex >= 0][0]  # finite end Voronoi vertex
            
            t = vor.points[pointidx[1]] - vor.points[pointidx[0]]  # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal
            
            midpoint = vor.points[pointidx].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[ii] + direction * ptp_bound.max()
            
            infinite_segments.append([vor.vertices[ii], far_point])
    
    
    
    ##This plots out the points as indexed in the previous section of code
    ##Loops through the line points as defined above. It then uses skimage draw to put single pixel lines on the canvas.
    for j in finite_segments:
        x0, y0=j[1]
        x1, y1=j[0]
        rr,cc=skdr.line(int(x0),int(y0),int(x1),int(y1))
        for jj in range(len(rr)):
            if rr[jj]<0 or cc[jj] <0 or rr[jj]>=canvas.shape[0] or cc[jj]>=canvas.shape[1]:
                continue
            canvas[rr[jj],cc[jj]]=0
    for j in infinite_segments:
        x0,y0=j[1]
        x1,y1=j[0]
        rr,cc=skdr.line(int(x0),int(y0),int(x1),int(y1))
        for jj in range(len(rr)):
            if rr[jj]<0 or cc[jj] <0 or rr[jj]>=canvas.shape[0] or cc[jj]>=canvas.shape[1]:
                continue
            canvas[rr[jj],cc[jj]]=0
    return canvas
