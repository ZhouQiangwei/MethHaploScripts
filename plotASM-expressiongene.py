# -*- coding: utf-8 -*-
"""
Created on Fri Sep 21 10:29:50 2018

@author: qwzhou
"""

"""
=======================================
plot line and dash
=======================================
ASM across the site
"""
import numpy as np
import matplotlib.pyplot as plt
import sys
from matplotlib.backends.backend_pdf import PdfPages
#x = np.linspace(0, 10, 500)

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
 
#              
#         
#readfile("/Users/qiangweizhou/K562.asmonGeneend.Methy.1.txt", y,z, k)
#readfile("/Users/qiangweizhou/IMR90.asmonGeneend.Methy.1.txt", y,z, k)
#readfile("/Users/qiangweizhou/HepG2.asmonGeneend.Methy.1.txt", y,z, k)
#readfile("/Users/qiangweizhou/A549.asmonGeneend.Methy.1.txt", y,z, k)
#readfile("/Users/qiangweizhou/HUES64.asmonGeneend.Methy.1.txt", y,z, k)
#readfile("/Users/qiangweizhou/GM12878.asmonGeneend.Methy.1.txt", y,z, k)

readfile("/Users/qiangweizhou/A549.asmonExpressionGene.Methy.1.txt.Aver", y,z, k)
readfile("/Users/qiangweizhou/A549.asmonUnGene.Methy.1.txt.Aver", y,z, k)
x = np.linspace(1, len(y[0]), len(y[0]))

label=['Expression', 'Unexpreesion']
#label=['K562', 'IMR90', 'HepG2', 'A549', 'HUES64', 'GM12878']
filename="ASMonGeneExpression.all"
filename2=filename + ".pdf"
nsample=2
legend=1
percentage=1
cutoff=6
######################################################
def find_martrix_max_value(data_matrix):  
    new_data=[]  
    for i in range(len(data_matrix)):  
        new_data.append(max(data_matrix[i]))
    return max(new_data)

xlen=len(y[0])
print(xlen, xlen/2)
#######################################################
def plotline(x, y, title, label, nsample, legend, filename):
    prosamp = 0
    fig, ax = plt.subplots()
    while prosamp < nsample:
        y[prosamp] = [i*percentage for i in y[prosamp]]
        #for i,item in enumerate(y[prosamp]):
        #    if item >6:
        #        y[prosamp][i]=6
        ax.plot(x, y[prosamp], label=label[prosamp]) #,color="dodgerblue"
        prosamp = prosamp +1
    #dashes = [10, 5, 100, 5]
    #line1.set_dashes(dashes)   #   dash line
    # Remove the plot frame lines. They are unnecessary here.
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    ax.xaxis.set_major_formatter(plt.FuncFormatter('{:.0f}'.format))
    #ax.yaxis.set_major_formatter(plt.FuncFormatter('{:.1f}%'.format))
    #plt.grid(True, 'major', 'y', ls='--', lw=.5, c='k', alpha=.3)
    plt.tick_params(axis='both', which='both', bottom=False, top=False,
                labelbottom=True, left=False, right=False, labelleft=True)
    #ax.axes.get_xaxis().set_visible(False)
    if legend == 1:
        plt.legend(loc='best', prop={'size': 12})   #   legend , loc is the legend location
    #plt.axvline(x=xlen/2-1, ls="--", color='black')
    
    plt.axhline(y=0, xmin=0.05, xmax=0.5, linewidth=8, color='gray')
    plt.axhline(y=0, xmin=0.5, xmax=0.505, linewidth=8, color='k' )
    plt.axhline(y=0, xmin=0.505, xmax=0.95, linewidth=8, color='gray') 
    scale_ls = [1,len(x)/2,len(x)]
    index_ls = ['-200bp','Start', "+200bp"]
    plt.xticks(scale_ls,index_ls,color='k', size=15)
    ax.set_title(title,size=15)
    ax.set_ylabel('ASM distribution',size=15)
    #ax.set_ylabel('Methylation Level',size=15)
    maxy = 100
    maxy = find_martrix_max_value(y) * 1.1
    
    ax.set_ylim(0.0, maxy)
    #plt.show()
    
    #filename2=filename + ".png"
    #plt.savefig(filename2, bbox_inches='tight')

#label = ['IMR90', 'A549', 'H1', 'GM12878', 'encodeA549']


pdf = PdfPages(filename2)

plotline(x, y, "ASM distribution", label, nsample, legend, filename+".CG")

#plotline(x, y, "CG methylation distribution", label, nsample, legend, filename+".CG")
legend=0
#plotline(x, z, "CHG methylation distribution", label, nsample, legend, filename+".CHG")
#plotline(x, k, "CHH methylation distribution", label, nsample, legend, filename+".CHH")
 
pdf.savefig()
pdf.close()

'''
    fig, ax = plt.subplots()
    line1, = ax.plot(x, k1, label='IMR90')
    #dashes = [10, 5, 100, 5]
    #line1.set_dashes(dashes)   #   dash line
    line2, = ax.plot(x, k2, label='A549')
    
    #  several dash line example
    #line3, = ax.plot(x, y3, ':', label='..style')
    #line4, = ax.plot(x,-np.sin(x)/2, '-.', label='-.style')
    #line5, = ax.plot(x,np.sin(x)/4, '--', label='--style')
    #line6, = ax.plot(x,-np.sin(x)/4, '^', label='--style')
    
    #plt.axis('off')
    #plt.xticks([])
    #plt.yticks([])
    #ax.axes.get_yaxis().set_visible(False)
    # Remove the plot frame lines. They are unnecessary here.
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    #ax.spines['left'].set_visible(False)
    
    #ax.axes.get_xaxis().set_visible(False)
    plt.legend(loc='center right')   #   legend , loc is the legend location
    plt.axhline(y=0, xmin=0.05, xmax=0.35, linewidth=8, color='gray')
    plt.axhline(y=0, xmin=0.65, xmax=0.35, linewidth=8, color='k' )
    plt.axhline(y=0, xmin=0.65, xmax=0.95, linewidth=8, color='gray') 
    scale_ls = [1,39,76,117]
    index_ls = ['upstream','Start','End', "downstream"]
    plt.xticks(scale_ls,index_ls,color='k', size=15)
    #ax.set_title('Box plot')
    ax.set_ylabel('Methylation Level',size=15)
    maxy=max(k1)
    if max(k2) > maxy:
        maxy = max(k2)*1.1
    else:
        maxy = maxy*1.1
    
    ax.set_ylim(0.0, maxy)
    #plt.savefig("test.png")
    plt.show()
 '''