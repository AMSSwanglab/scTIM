#!/usr/bin/env python 3.6.5
# encoding: utf-8
"""
TIM.py

Created by Zhanying Feng on 2018-07-24.
Copyright (c) 2018 Zhanying Feng All rights reserved.
"""

import sys, getopt
from multiprocessing import Process, Queue
import numpy as np
import random, math, datetime


def Hbeta(D = np.array([]), beta = 1/10000000):
    """Compute the perplexity and the P-row for a specific value of the precision of a Gaussian distribution."""

    # Compute P-row and corresponding perplexity
    D = D.astype('float64');P = np.exp(-D.copy() * beta);
    sumP = sum(P);
    if sumP==0:
        print("!!")
    H = np.log(sumP) + beta * np.sum(D * P) / sumP;
    P = P / sumP;
    return H, P;

def x2p(var, X = np.array([]), tol = 1e-5, perplexity = 30.0):
    """Performs a binary search to get P-values in such a way that each conditional Gaussian has the same perplexity."""

    # Initialize some variables
    # Computing pairwise distances
    (n, d) = X.shape;X = X.astype('float64');
    sum_X = np.sum(np.square(X), 1);
    D = np.add(np.add(-2 * np.dot(X, X.T), sum_X).T, sum_X);
    
    P = np.zeros((n, n));
    beta = np.ones((n, 1))/1000000000;
    logU = np.log(perplexity);

    # Loop over all datapoints
    for i in range(n):

        # Compute the Gaussian kernel and entropy for the current precision
        betamin = -np.inf;
        betamax =  np.inf;
        Di = D[i, np.concatenate((np.r_[0:i], np.r_[i+1:n]))];
        (H, thisP) = Hbeta(Di, beta[i]);

        # Evaluate whether the perplexity is within tolerance
        Hdiff = H - logU;
        tries = 0;
        while np.abs(Hdiff) > tol and tries < 1000:

            # If not, increase or decrease precision
            if Hdiff > 0:
                betamin = beta[i].copy();
                if betamax == np.inf or betamax == -np.inf:
                    beta[i] = beta[i] * 2;
                else:
                    beta[i] = (beta[i] + betamax) / 2;
            else:
                betamax = beta[i].copy();
                if betamin == np.inf or betamin == -np.inf:
                    beta[i] = beta[i] / 2;
                else:
                    beta[i] = (beta[i] + betamin) / 2;

            # Recompute the values
            (H, thisP) = Hbeta(Di, beta[i]);
            Hdiff = H - logU;
            tries = tries + 1;

        # Set the final row of P
        P[i, np.concatenate((np.r_[0:i], np.r_[i+1:n]))] = thisP;

    # Return final P-matrix
    for j in range(X.shape[0]):
        var[j] = 1/beta[j]
    return P;

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
def FS(x):
    g = open('WarningMesg','a')
    g.write(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + ": Feature Specificity Start!\n")
    g.close()
    Entro = [];x = x.astype('float64');x.sort()
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
    g = open('WarningMesg','a')
    g.write(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + ": Feature Specificity Done!\n")
    g.close()
    return np.array(FS)
def Dis(x):
    g = open('WarningMesg','a')
    g.write(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + ": Sample DisMatrix Start!\n")
    g.close()
    x = x.astype('float64');xt = x.T;eps = 0.0000001;xt += eps
    var = [0]*xt.shape[0]
    x2p(var,xt)
    pb = np.zeros((xt.shape[0],xt.shape[0]))
    for k in range(xt.shape[0]):
        num = 0;
        for m in range(xt.shape[0]):
            if m != k:
                num += np.exp(-np.square(np.linalg.norm(xt[m]-xt[k]))/(2*var[k][0]))
        for l in range(xt.shape[0]):
            if l != k :
                pb[l][k] = np.exp(-np.square(np.linalg.norm(xt[l]-xt[k]))/(2*var[k][0]))/num
    pa = np.zeros((xt.shape[0],xt.shape[0]))
    for i in range(xt.shape[0]):
        for j in range(xt.shape[0]):
            if i != j :
                pa[i][j] = (pb[i][j] + pb[j][i])/2
    g = open('WarningMesg','a')
    g.write(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + ": Sample DisMatrix Done!\n")
    g.close()
    return pa
def Red(x):
    g = open('WarningMesg','a')
    g.write(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + ": Feature Redundancy Start!\n")
    g.close()
    x = x.astype('float64');xt = x.T + eps
    xt = xt/np.sum(xt,0)
    Redu = []
    for i in range(xt.shape[1]):
        xtc = np.array([list((xt.T)[i])]).T;
        c = xt + xtc
        JSD_i = 0.5 * np.sum(xtc * np.log(2 * xtc / c) + xt * np.log( 2 * xt / c),0)
        Redu.append(JSD_i)
    g = open('WarningMesg','a')
    g.write(datetime.datetime.now().strftime('%Y-%m-%d-%H-%M-%S') + ": Feature Redundancy Done!\n")
    g.close()
    return np.array(Redu)

def TM(x,z,y):
    x = x.astype('float64');xs = (x * (np.array([y]).T)).T
    Q = 1 / (1 + np.asarray(EuclideanDistances(xs,xs))**2)
    dom = np.sum(Q) - np.sum(np.diag(Q))
    Q = Q / dom;
    Q[range(Q.shape[0]),range(Q.shape[0])] = 1
    return np.sum(-z*np.log2(Q))/x.shape[1]

def C(x,y,d,z,alpha,beta,gamma,w):
    x = x.astype('float64');n = np.sum(w)
    C1 = alpha * np.dot(np.dot(np.array([w]),y),np.array([w]).T)[0][0]/(n*(n-1))
    C2 = beta * TM(x,d,w)
    C3 = gamma * np.dot(np.array([w]),z)[0]/n
    return (C1 - C2 + C3 ) * 500000

def SimAnealing(x,F,D,R,alpha,beta,gamma,o):
    x = x.astype('float64');F = F.astype('float64');D = D.astype('float64');R = R.astype('float64');
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
    o.put(np.array(w))

def TIM(inputFile,outputFile,norm,parameter=[0.1,0.4,0.5]):
    alpha = float(parameter[0]);beta = float(parameter[1]);gamma = float(parameter[2])
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
    data = np.array(data);data = data.astype('float64')
    if norm == 'y':
        data = np.log2(np.array(data)+1)
  
   
    fs = FS(data);sample_dis = Dis(data);red = Red(data);
    # selecting landmark gene -- 10 processes
    Q = Queue()
    SimInput = (data,fs,sample_dis,red,alpha,beta,gamma,Q)
    P1 = Process(target = SimAnealing, args = SimInput)
    P2 = Process(target = SimAnealing, args = SimInput)
    P3 = Process(target = SimAnealing, args = SimInput)
    P4 = Process(target = SimAnealing, args = SimInput)
    P5 = Process(target = SimAnealing, args = SimInput)
    P6 = Process(target = SimAnealing, args = SimInput)
    P7 = Process(target = SimAnealing, args = SimInput)
    P8 = Process(target = SimAnealing, args = SimInput)
    P9 = Process(target = SimAnealing, args = SimInput)
    P10 = Process(target = SimAnealing, args = SimInput)
    P1.start();P2.start();P3.start();P4.start();P5.start();
    P6.start();P7.start();P8.start();P9.start();P10.start();
    P1.join();P2.join();P3.join();P4.join();P5.join();
    P6.join();P7.join();P8.join();P9.join();P10.join();
    Results = []
    for i in range(10):
        Results.append(Q.get())
    output = np.sum(Results,0)

    g = open(outputFile,'w')
    for i in range(output.shape[0]):
        if output[i] == 2:#10:
            g.write(GeneSet[i]+'\n')
    g.close()
def main(argv=None):
    if argv is None:
        argv = sys.argv
    parameter = ""
    inputFile = ""
    outputFile = ""
    norm = ""
    try:
        opts, args = getopt.getopt(argv[1:], "hp:i:o:a:", ["help","parameter", "input", "output","normalization"])
    except (getopt.GetoptError):
        print('USAGE: TIM -p <alpha,beta,gamma> -i <input> -o <output> -a <normalization>')

    # option processing
        
    for option, value in opts:
        if option in ("-h", "--help"):
            print('USAGE: TIM -p <alpha,beta,gamma> -i <input> -o <output> -a <normalization>')
            sys.exit()
        if option in ("-p", "--parameter"):
            parameter = value.split(",")
        if option in ("-i", "--input"):
            inputFile = value
        if option in ("-o", "--output"):
            outputFile = value
        if option in ("-a", "--normalization"):
            norm = value
    try:
        f = open(inputFile, 'r')
    except (IOError):
        f = open(inputFile, 'w');f.close()

    # input parameters, data and output file.
    TIM(inputFile,outputFile,norm,parameter)

    # run model
    
if __name__ == "__main__":
    eps = 0.00000000001;global fs;global sample_dis;global red; global data
    sys.exit(main())
