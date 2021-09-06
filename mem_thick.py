#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  6 19:58:25 2021

@author: huangjianxiang
"""

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np

u=mda.Universe("nh3.pdb")

center=np.array([207.31954956054688, 208.61289978027344, 205.3704071044922])


pos=u.select_atoms("all").positions

extend=np.tile(center,7200)

center=extend.reshape(np.shape(pos))

xy=pos[:,0:2] 

cxy=center[:,0:2] 

cos_theta=(pos[:,0]-center[:,0])/np.linalg.norm(xy-cxy,axis=1)
cos_phi=np.linalg.norm(xy-cxy,axis=1)/np.linalg.norm(pos-center,axis=1)

theta_degree=np.zeros((7200,2))

y_vec=pos[:,1]-center[:,1]
z_vec=pos[:,2]-center[:,2]

for i in range(7200):
    if y_vec[i]>=0:
        theta=np.degrees(np.arccos(cos_theta[i]))
    else:
        theta=360-np.degrees(np.arccos(cos_theta[i]))
    if z_vec[i]>=0:
        phi=np.degrees(np.arccos(cos_phi[i]))
    else:
        phi=-1*np.degrees(np.arccos(cos_phi[i]))
    
    theta_degree[i:,]=[theta,phi]
    
xedges=np.arange(0,361,20)
#yedges=np.arange(-90,91,15)
yedges=np.array([-90, -60, -45, -35, -25, -10,   0,  10,  25,  35,  45,  60,  90])
    
bins=[xedges,yedges]

H,edges=np.histogramdd(theta_degree,bins)  

data=[]

for i in range(18):
#    record=np.arange(0,7200)
#    record1=record[theta_degree[:,0]>xedges[i]]
##    print()
##    record=theta_degree[:,0][record1]
#    record2=record[theta_degree[:,0]<=xedges[i+1]]
#    print(np.size(record2))
    for j in range(12):
#        pass
#        theta_degree[:,1]record2
        data1=[]
        for k in range(7200):
            if theta_degree[k,0]>xedges[i] and theta_degree[k,0]<=xedges[i+1] and theta_degree[k,1]>yedges[j] and theta_degree[k,1]<=yedges[j+1]:
                data1.append(k)
        data.append(data1)
            

###divide the beads into two parts

center=np.array([207.31954956054688, 208.61289978027344, 205.3704071044922])

head=u.select_atoms("all")

temp=np.zeros(7200)
for i in range(18):
#    record=np.arange(0,7200)
#    record1=record[theta_degree[:,0]>xedges[i]]
##    print()
##    record=theta_degree[:,0][record1]
#    record2=record[theta_degree[:,0]<=xedges[i+1]]
#    print(np.size(record2))

    for j in range(12):
        extend1=np.tile(center,len(data[12*i+j]))
        center1=extend1.reshape((len(data[12*i+j]),3))
#        print(pos[data[12*i+j]])
        radius=np.linalg.norm(pos[data[12*i+j]]-center1,axis=1)
        ave=np.mean(radius)
        
        half1=radius[radius>ave]
        half2=radius[radius<=ave]
#        print(half2,half1)
#        print(12*i+j,ave,np.mean(half1)-np.mean(half2))
        thickness=(np.mean(half1)-np.mean(half2))*0.1
        
        list1=data[12*i+j]
        for k in range(len(list1)):
            temp[list1[k]]=thickness
print(temp)

head.tempfactors=temp
print(head.tempfactors)
head.write("nh3_thickness.pdb")
        
        
        
        
    
    
    
    
    
    
    
    
    