#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul  6 22:34:08 2024

@author: Piero Wemyss


function props = fetchProps(params2fetch,comps)

"""

import numpy as np
from dict2struct import dict2struct

def orgProps(params2fetch,comps,selected_comps,allProps):

    compList = comps[:,np.newaxis]
    
    if params2fetch == 1:

        antoine = allProps.antoine

    if params2fetch == 2:
        
        NRTL_aij = allProps.NRTL_aij
        
        NRTL_bij = allProps.NRTL_bij
        
        NRTL_cij = allProps.NRTL_cij
    
    if params2fetch == 3:

        PLXANT = allProps.PLXANT
    
    if params2fetch == 4:
        crit = np.column_stack((allProps.TcCel, allProps.Pc, allProps.omega))
    
    compInd = np.zeros([len(selected_comps),1],dtype=int)
    for i in range(0,len(selected_comps)):
        compInd[i] = np.where(compList == selected_comps[i])[0][0]
    
    compInd = compInd[:,0]
    
    if params2fetch == 1:
        pantoine = np.zeros([len(selected_comps),3])
        for i in range(0,len(compInd)):
            pantoine[i,:] = antoine[compInd[i],:];
        props = dict2struct({'antoine':pantoine})
    
    if params2fetch == 2:
        pNRTL_aij = np.zeros([len(NRTL_aij),len(selected_comps)])
        pNRTL_bij = np.zeros([len(NRTL_bij),len(selected_comps)])
        pNRTL_cij = np.zeros([len(NRTL_cij),len(selected_comps)])
        for i in range(0,len(compInd)):
            pNRTL_aij[:,i] = NRTL_aij[:,compInd[i]]
            pNRTL_bij[:,i] = NRTL_bij[:,compInd[i]]
            pNRTL_cij[:,i] = NRTL_cij[:,compInd[i]]
                  
        for i in range(0,len(compInd)):
            pNRTL_aij[i,:] = pNRTL_aij[compInd[i],:]
            pNRTL_bij[i,:] = pNRTL_bij[compInd[i],:]
            pNRTL_cij[i,:] = pNRTL_cij[compInd[i],:]
            
        pNRTL_aij = pNRTL_aij[0:len(compInd),:]
        pNRTL_bij = pNRTL_bij[0:len(compInd),:]
        pNRTL_cij = pNRTL_cij[0:len(compInd),:]
        
        props = dict2struct({
            'NRTL_aij':pNRTL_aij,
            'NRTL_bij':pNRTL_bij,
            'NRTL_cij':pNRTL_cij,
            })

    if params2fetch == 3:
        pPLXANT = np.zeros([len(selected_comps),7])
        for i in range(0,len(compInd)):
            pPLXANT[i,:] = PLXANT[compInd[i],:];
        props = dict2struct({'PLXANT':pPLXANT})

    
    if params2fetch == 4:
        pcrit = np.zeros([len(selected_comps),3])
        for i in range(0,len(compInd)):
            pcrit[i,:] = crit[compInd[i],:];
        props = dict2struct({'TcCel':pcrit[:,0],
                             'Pc':pcrit[:,1],
                             'omega':pcrit[:,2]})

    return props
