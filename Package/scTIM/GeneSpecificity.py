# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:25:38 2018

@author: Zhanying Feng
"""
import numpy as np
import math
eps = 0.00000000001

def JSD(x,y):
	return 0.5*np.sum(x*np.log(2*x/(y+x))+y*np.log(2*y/(x+y)),0)
    
def get_gs(x):
    y = x.copy();z = x.copy()
    y = np.delete(y,0);z = np.delete(z,z.shape[0]-1)
    gap = y-z
    max_gap = np.argwhere(gap == max(gap))[0][0];
    vector = (np.array(list(range(x.shape[0]))) > max_gap).astype('int') + eps
    vector = (vector / np.sum(vector));          
    return -np.log2(JSD(x,vector))*((max_gap+1)/x.shape[0])
def EuclideanDistances(A, B):
    BT = B.transpose()
    vecProd = np.dot(A,BT)
    SqA =  A**2
    sumSqA = np.matrix(np.sum(SqA, axis=1))
    sumSqAEx = np.tile(sumSqA.transpose(), (1, vecProd.shape[1]))

    SqB = B**2
    sumSqB = np.sum(SqB, axis=1)
    sumSqBEx = np.tile(sumSqB, (vecProd.shape[0], 1))
    SqED = sumSqBEx + sumSqAEx - 2*vecProd
    SqED[SqED<0]=0.0
    ED = np.sqrt(SqED)
    return ED
def GeneSpecificity(x):
    Entro = [];x.sort()
    for i in range(x.shape[0]):
        n = math.ceil(max(x[i])-min(x[i]))+1;
        count = [0]*n;
        k = 0;cut = min(x[i]) + k + 1
        for j in range(x.shape[1]):
            if x[i][j] >= cut:
                k += 1; cut = min(x[i]) + k + 1
                count[k] += 1
            else:
                count[k] += 1
        while 0 in count:
            count.remove(0)
        entro = 0;s = sum(count);
        if s != x.shape[1]:
            print(i,s)
        count = np.array(count)/s
        entro = -np.sum(count * np.log2(count))
        Entro.append(entro)
    eps = 0.0000001;x += eps
    x = x/np.array([list(np.sum(x,1))]).T
    FS = []
    for i in range(x.shape[0]):
        FS.append(get_gs(x[i])*Entro[i])
    return np.array(FS)
    