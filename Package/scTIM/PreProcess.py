# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:44:58 2018

@author: Zhanying Feng
"""
import numpy as np

def PreProcess(inputFile,logornot='y'):
    f = open(inputFile,'r');
    data = f.readlines(); f.close();
    del data[0]
    GeneSet = []
    for i in range(len(data)-1,-1,-1):
        data[i] = data[i].strip('\n');data[i] = data[i].split('\t');
        s = 0
        for j in range(len(data[i])-1):
            s += float(data[i][j+1])
        if s >= 2:
            GeneSet.insert(0,data[i][0])
            del data[i][0]
            for k in range(len(data[i])):
                data[i][j] = float(data[i][j])
        else:
            del data[i]
    data = np.array(data).astype('float')
    if logornot == 'y':
        data = np.log2(np.array(data) + 1)
    return data, GeneSet