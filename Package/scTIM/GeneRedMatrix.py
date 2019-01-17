# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:24:05 2018

@author: Zhanying Feng
"""
import numpy as np
eps = 0.00000000001

def GeneRedMatrix(x):#,y):
    xt = x.T + eps
    xt = xt/np.sum(xt,0)
    Redu = []
    for i in range(xt.shape[1]):
        xtc = np.array([list((xt.T)[i])]).T;
        c = xt + xtc
        JSD_i = 0.5 * np.sum(xtc * np.log(2 * xtc / c) + xt * np.log( 2 * xt / c),0)
        Redu.append(JSD_i)
    return np.array(Redu)
    