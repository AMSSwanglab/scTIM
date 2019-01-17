# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:25:13 2018

@author: Zhanying Feng
"""
import numpy as np
import sys
eps = 0.00000000001

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
def CellRedMatrix(x):
    xt = x.T;xt += eps
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
    return pa