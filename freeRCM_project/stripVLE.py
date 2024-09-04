#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 00:32:00 2024

@author: Piero Wemyss
"""

import numpy as np
import nifco as ni
from antoineCalc import antoineCalc

def stripVLE(Y, x, P, comps, props,  opts):
    NComps = len(comps)
    y = np.matrix(Y[0:NComps])
    T = Y[NComps]
    Psat = antoineCalc(T,comps,props,opts.antMethod)

    if opts.activity == 1:
        gamma = np.ones([NComps,1])
    if opts.activity >= 2:
        gamma = ni.nrtl(x,T,props.NRTL_aij,props.NRTL_bij,props.NRTL_cij,NComps)
        
    if opts.activity >= 3:
        phi = ni.srk(y,T,P,props.TcCel,props.Pc,props.omega,NComps)
    else:
        phi = np.ones([NComps,1])
    
    raoult = np.zeros([NComps,1])
    dalton = np.zeros([NComps,1])
    y = np.array(y)[0]

    for i in range(0,NComps):
        raoult[i,:] = x[i]*gamma[i]*Psat[i] - y[i]*phi[i]*P
        dalton[i] = y[i]*phi[i]/(gamma[i]*Psat[i])

    if NComps <= 3:
        tray = np.append(raoult,[[np.sum(np.abs(y)) - float(1.00), np.sum(y) - float(1.00)]])
        
    if NComps > 3:
        tray = [raoult,
                # 1/np.sum(dalton) - P,
                np.sum((y[1:NComps])) - 1]
    
    return tray
