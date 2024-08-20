#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 11:57:36 2024

@author: Piero Wemyss
"""

from matplotlib.figure import Figure

def RCMplot(x,comps,opts):
    
    try:
        n_it = opts.n_it*1
    except TypeError:
        n_it = 100
    except AttributeError:
        n_it = 100
    try:
        linewidth= opts.linewidth*1
    except TypeError:
        linewidth = 1.2 
    except AttributeError:
        linewidth = 1.2 
        
    arr_pos = int(n_it-(0.1*n_it))
    
    fig = Figure(figsize=(5, 4), dpi=200)
    ax = fig.add_subplot(111)
    for i in range(0,len(x[0,0,:])):
        line, = ax.plot(x[:, 0, i], x[:, 1, i],linewidth=linewidth)
        color = line.get_color()
        if len(x[:,0,0]) > 2:
            ax.annotate('', xy=(x[arr_pos,0,i],x[arr_pos,1,i]),
                     xytext=(x[arr_pos+1,0,i],x[arr_pos+1,1,i]), 
                     arrowprops=dict(facecolor=color, edgecolor=color, 
                                     shrink=0.5, headwidth=3, headlength=3)) 
    ax.plot([0,1],[1,0],'k-',linewidth=1.5)
    ax.plot([0,0],[0,1],'k-',linewidth=1.5)
    ax.plot([0,1],[0,0],'k-',linewidth=1.5)
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    ax.axis('off')
    ax.annotate(comps[2],[0,0],[-0.025,-0.04],fontsize=7)
    ax.annotate(comps[1],[0,1],[-0.02,1.02],fontsize=7)
    ax.annotate(comps[0],[1,0],[1.0,-0.02],fontsize=7)
    ax.annotate('Generated using freeRCM by Piero Wemyss',[0.7,1.0],[0.7,1.0],fontsize=3,color='tab:gray')

    return fig, ax
