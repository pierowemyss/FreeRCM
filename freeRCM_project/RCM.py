#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 21:42:34 2024

@author: Piero Wemyss
"""

import numpy as np
from orgProps import orgProps
from antoineCalc import antoineCalc
from stripVLE import stripVLE
from scipy.optimize import root
from dict2struct import dict2struct

def RCM(comps,selected_comps,P,allProps,opts,x0n,genOpt):

    try:
        dxi = opts.dxi*1
    except TypeError:
        dxi = 0.05
    except AttributeError:
        dxi = 0.05
        
    try:
        n_it = opts.n_it*1
    except TypeError:
        n_it = 100
    except AttributeError:
        n_it = 100

    NComps = len(selected_comps)    
    if opts.antMethod == 1:
        props = orgProps(1,comps,selected_comps,allProps)
    else:
        props = orgProps(3,comps,selected_comps,allProps)

    if opts.activity == 2:
        props2 = orgProps(2,comps,selected_comps,allProps)
        for key in props2:
            props[key] = props2[key]

    if opts.activity == 3:
        props2 = orgProps(2,comps,selected_comps,allProps)
        props3 = orgProps(4,comps,selected_comps,allProps)
        for key in props2:
            props[key] = props2[key]
        for key in props3:
            props[key] = props3[key]

    if genOpt == 2:
        xPlot = np.zeros([2*n_it,NComps,1])
        yPlot = np.zeros([2*n_it,NComps,1])
        TPlot = np.zeros([2*n_it,1,1])
        x0 = x0n

        x = np.zeros([2*n_it,NComps])
        y = np.zeros([2*n_it,NComps])
        T = np.zeros([2*n_it,1])

        x[0,:] = x0
        x0Max = np.max(x0)
        majComp = int(np.where(x0 == x0Max)[0][0])
  
        if opts.antMethod == 1:
            condAntC = orgProps(1,comps,np.array([selected_comps[majComp]]),allProps)
            T0 = (condAntC.antoine[0][1]/(condAntC.antoine[0][0]-np.log10(P)) - condAntC.antoine[0][2])
            T0 = root(lambda T: antoineCalc(T,np.array([selected_comps[majComp]]),props,opts.antMethod)-P,T0,method='lm')
            Y0 = np.append(x0, T0[0])
        else:
            try:
                Tn = root(lambda T: antoineCalc(T,np.array([selected_comps[majComp]]),props,opts.antMethod)-P,50+(3*P),method='lm')
                T0 = Tn.x
                Y0 = np.append(x0, T0[0])
            except Tn.success == False:
                inc = 20
                while Tn.success == False:
                    Tn = root(lambda T: antoineCalc(T,np.array([selected_comps[majComp]]),props,opts.antMethod)-P,inc+(3*P),method='lm')
                    inc = inc + 0.5

                T0 = Tn.x
                Y0 = np.append(x0, T0[0])

        for i in range(0,n_it):
            
            Yn = root(lambda Y: stripVLE(Y,x[i,:],P,selected_comps,props,opts),Y0,method='lm',tol=1e-9,options=opts.lmopts)
            y[i,:] = Yn.x[0:NComps]
            T[i] = Yn.x[NComps]
            Y0 = np.append(y[i,:], T[i])
            # print('\nDEBUG: Yn (1st loop) = ', Y0,'\n')
            if Yn.success == False:
                print('\nWARNING: FALSE EXIT FLAG on iteration',i,'of line\n')
  
            x[i+1,:] = x[i,:] + dxi*(x[i,:] - y[i,:])
  
        x[0:n_it] = np.flip(x[0:n_it],0)
        y[0:n_it] = np.flip(y[0:n_it],0)
  
        Y0 = np.append(x0, T0[0])
  
        for i in range(n_it-1,2*n_it-1):
    
            Yn = root(lambda Y: stripVLE(Y,x[i,:],P,selected_comps,props,opts),Y0,method='lm',tol=1e-9,options=opts.lmopts)
            y[i,:] = Yn.x[0:NComps]
            T[i] = Yn.x[NComps]
            Y0 = np.append(y[i,:], T[i])
            if Yn.success == False:
                print('\nWARNING: FALSE EXIT FLAG on iteration',i,'\n')
    
            x[i+1,:] = x[i,:] - dxi*(x[i,:] - y[i,:])
    
        xPlot[:,:,0] = x
        yPlot[:,:,0] = y
        TPlot[:,:,0] = T

    else:
        xPlot = np.zeros([2*n_it,NComps,opts.lines])
        yPlot = np.zeros([2*n_it,NComps,opts.lines])
        TPlot = np.zeros([2*n_it,1,opts.lines])
    
        # x_bank = np.linspace(0.02, 0.64, opts.lines)
        x_bank = np.linspace(0.27, 0.49, opts.lines)
    
        for k in range(0, opts.lines):
            
            if NComps == 3:
                # x0 = np.array([x_bank[k],0.95-1.5*x_bank[k],0.05]+x_bank[k])
                x0 = np.array([x_bank[k],1.5-3*x_bank[k],-0.5+2*x_bank[k]])
            
            if NComps == 4:
                if opts.split == 1:
                    x0 = np.array([float(k)/float(opts.lines), 0.95-float(k)/float(opts.lines), 0.025, 0.025],dtype=float)
                
                if opts.split == 3:
                    x0 = np.array([0.025, float(k)/float(opts.lines), 0.95-float(k)/float(opts.lines), 0.025],dtype=float)
                
                if opts.split == 3:
                    x0 = np.array([0.025, 0.025, float(k)/float(opts.lines), 0.95-float(k)/float(opts.lines)],dtype=float)
                
            
            if np.sum(x0) >= 1.00001:
                raise Exception('x0 does not sum to 1 (sum is ',np.sum(abs(x0)),'), input valid value for opts.lines')
            
            x = np.zeros([2*n_it,NComps])
            y = np.zeros([2*n_it,NComps])
            T = np.zeros([2*n_it,1])
        
            x[0,:] = x0
            x0Max = np.max(x0)
            majComp = int(np.where(x0 == x0Max)[0][0])
      
            if opts.antMethod == 1:
                condAntC = orgProps(1,comps,np.array([selected_comps[majComp]]),allProps)
                T0 = (condAntC.antoine[0][1]/(condAntC.antoine[0][0]-np.log10(P)) - condAntC.antoine[0][2])
                T0 = root(lambda T: antoineCalc(T,np.array([selected_comps[majComp]]),props,opts.antMethod)-P,T0,method='lm')
                Y0 = np.append(x0, T0[0])
            else:
                try:
                    Tn = root(lambda T: antoineCalc(T,np.array([selected_comps[majComp]]),props,opts.antMethod)-P,50+(3*P),method='lm')
                    T0 = Tn.x
                    Y0 = np.append(x0, T0[0])
                except Tn.success == False:
                    inc = 20
                    while Tn.success == False:
                        Tn = root(lambda T: antoineCalc(T,np.array([selected_comps[majComp]]),props,opts.antMethod)-P,inc+(3*P),method='lm')
                        inc = inc + 0.5

                    T0 = Tn.x
                    Y0 = np.append(x0, T0[0])


            
            for i in range(0,n_it):
                
                Yn = root(lambda Y: stripVLE(Y,x[i,:],P,selected_comps,props,opts),Y0,method='lm',tol=1e-9,options=opts.lmopts)
                y[i,:] = Yn.x[0:NComps]
                T[i] = Yn.x[NComps]
                Y0 = np.append(y[i,:], T[i])
                if Yn.success == False:
                    print('\nWARNING: FALSE EXIT FLAG on iteration',i,'of line',k+1,'\n')
      
                x[i+1,:] = x[i,:] + dxi*(x[i,:] - y[i,:])
            
            x[0:n_it] = np.flip(x[0:n_it],0)
            y[0:n_it] = np.flip(y[0:n_it],0)
      
            Y0 = np.append(x0, T0[0])
      
            for i in range(n_it-1,2*n_it-1):
        
                Yn = root(lambda Y: stripVLE(Y,x[i,:],P,selected_comps,props,opts),Y0,method='lm',tol=1e-9,options=opts.lmopts)
                y[i,:] = Yn.x[0:NComps]
                T[i] = Yn.x[NComps]
                Y0 = np.append(y[i,:], T[i])
                if Yn.success == False:
                    print('\nWARNING: FALSE EXIT FLAG on iteration',i,'of line',k+1,'\n')
        
                x[i+1,:] = x[i,:] - dxi*(x[i,:] - y[i,:])
      
            xPlot[:,:,k] = x
            yPlot[:,:,k] = y
            TPlot[:,:,k] = T
      
    curves = dict2struct()
    curves.x = xPlot
    curves.y = yPlot
    curves.T = TPlot

    return curves
