# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 15:37:26 2018

@author: qwzhou

This is a heatmap visualization script for asm on repeat/histone peaks.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import cm 
from matplotlib import axes

y = []
z = []
k = []

def readfile(filename, y, z, k):
    nline=0
    with open(filename, 'r') as fig:
        for line in fig:
            data = line.split()
            if nline == 0:
                y.append(list(map(float,data[1:])))
            elif nline == 1:
                z.append(list(map(float,data[1:])))
            elif nline == 2:
                k.append(list(map(float,data[1:])))
            nline=nline+1
            
#readfile("./IMR90/hg38.fa.out.repeat.start.sites.Methy.1.txt.Aver", y, z, k)
readfile("./hg38.fa.out.repeat.LINE.sites.Methy.1.txt.Aver", y, z, k) 
readfile("./hg38.fa.out.repeat.SINE.sites.Methy.1.txt.Aver", y, z, k) 
readfile("./hg38.fa.out.repeat.LTR.sites.Methy.1.txt.Aver", y, z, k)
readfile("./hg38.fa.out.repeat.Others.sites.Methy.1.txt.Aver", y, z, k)



filename="ASMonTE.heatmap"
filename2=filename + ".pdf"
pdf = PdfPages(filename2)
xlabels=[]
#ylabels=[ 'H3k27ac', 'H3k4me3', 'H3k9ac' ,  'H3K4me1' , 'H3k27me3', 'H3k9me3', 'H3k36me3' ]
ylabels=[ 'LINE', 'SINE', 'LTR', 'Others']
ylabels=['repeat']
title="ASM distribution"

def draw_heatmap(data,title, xlabels,ylabels):
    #cmap=cm.Blues    
    cmap=cm.get_cmap('RdYlBu_r')
    figure=plt.figure(facecolor='w')
    ax=figure.add_subplot(1,1,1,position=[0.1,0.15,0.8,0.8])
    ax.set_yticks(range(len(ylabels)))
    ax.set_yticklabels(ylabels)
    ax.set_xticks(range(len(xlabels)))
    ax.set_xticklabels(xlabels)
    vmax=data[0][0]
    vmin=data[0][0]
    for i in data:
        for j in i:
            if j>vmax:
                vmax=j
            if j<vmin:
                vmin=j
    map=ax.imshow(data,interpolation='nearest',cmap=cmap,aspect='auto',vmin=vmin,vmax=vmax-0.32)
    cb=plt.colorbar(mappable=map,cax=None,ax=None,shrink=0.5)
    ax.set_title(title,size=15)
    #plt.show()
    #plt.savefig(filename3, bbox_inches='tight')
    pdf.savefig()
    
            
#a=np.random.rand(10,10)

draw_heatmap(y,title, xlabels,ylabels) 
pdf.close()

