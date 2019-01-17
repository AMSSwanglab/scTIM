# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:26:20 2018

@author: Zhanying Feng
"""

import numpy as np
import random
eps = 0.00000000001

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
def TM(x,z,y):
    xs = (x * (np.array([y]).T)).T
    Q = 1 / (1 + np.asarray(EuclideanDistances(xs,xs))**2)
    dom = np.sum(Q) - np.sum(np.diag(Q))
    Q = Q / dom;
    Q[range(Q.shape[0]),range(Q.shape[0])] = 1
    return np.sum(-z*np.log2(Q))/x.shape[1]

def C(x,y,d,z,alpha,beta,gamma,w):
    n = np.sum(w)
    C1 = alpha * np.dot(np.dot(np.array([w]),y),np.array([w]).T)[0][0]/(n*(n-1))
    C2 = beta * TM(x,d,w)
    C3 = gamma * np.dot(np.array([w]),z)[0]/n
    return (C1 - C2 + C3 ) * 500000

def ExtractGene(x,F,D,R,alpha,beta,gamma):
    F = F/(max(F)-min(F))
    R = R/(np.max(R - np.diag([100]*R.shape[0]))-np.min(R + np.diag([100]*R.shape[0])))
    w = []
    for i in range(x.shape[0]):
        w.append(random.sample([0,1],1)[0])

    ite = 100  # iteration number

    tc = C(x,R,D,F,alpha,beta,gamma,w);
    t = 100; tmin = 1/100000000; delta = 0.98
    choiceset = list(range(0, x.shape[0]))
    while t > tmin:
        for i in range(ite):
            # Get a new solution from neighbor
            change = random.sample(choiceset, 1)[0]
            w1 = w.copy(); w1[change] = 1 - w[change]

            # Decide if current solution should be replaced

            nc = C(x,R,D,F,alpha,beta,gamma,w1)
            if nc > tc:
                tc = nc; w = w1
            else:
                if np.random.rand() < np.exp(-(np.abs(nc-tc))/t):
                    tc = nc; w = w1
        t = t * delta
    return np.array(w)