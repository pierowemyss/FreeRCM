#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  7 23:55:06 2024

@author: Piero Wemyss

%% Instructions
% 
% method:       1 - Regular Antoine Equation 
%               2 - Extended Antoine Equation (PLXANT, Aspen Plus)
% 
%               HIGHLY RECOMMENDED TO USE PLXANT
% 
%% 

"""

import numpy as np

def antoineCalc(T,comps,props,method):

    if method == 1:
        # props = fetchProps(1,comps)
        antoineCoeffs =  props.antoine
        logP = np.empty(len(comps))
        for i in range(0,len(comps)):
            logP[i] = antoineCoeffs[i,0] - antoineCoeffs[i,1]/(antoineCoeffs[i,2]+T)
    
        Psat = 10**logP
    
    if method == 2:
        # props = fetchProps(3,comps)
        antXCoeffs =  props.PLXANT
        T = T + 273.5      # Both T and P are absolute for this method
        lnP = np.empty(len(comps))
        for i in range(0,len(comps)):
            lnP[i] = (antXCoeffs[i,0] + antXCoeffs[i,1]/(antXCoeffs[i,2]+T)
                    + antXCoeffs[i,3]*T + antXCoeffs[i,4]*np.log(T) +
                      antXCoeffs[i,5]*T**antXCoeffs[i,6])
        
        Psat = np.exp(lnP)
        
    return Psat
