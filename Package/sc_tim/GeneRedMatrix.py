# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 19:24:05 2018

@author: Zhanying Feng
"""
import numpy as np
import multiprocessing

eps = 0.00000000001

def Partition(n,x,m):
	s = 1; e = int(n/m)
	Par = []
	for i in range(m-1):
		Par.append((s+i*e-1,e+i*e-1,x))
	Par.append((s+(m-1)*e-1,n-1,x))
	return Par

def ProcessPartition(arg):
	begin, end, whole_x = arg
	Redu = []
	for i in range(begin,end+1):
		x_i = whole_x[:,i].reshape(whole_x.shape[0],1)
		c = whole_x + x_i
		JSD_i = 0.5 * np.sum(x_i * np.log(2 * x_i / c) + whole_x * np.log( 2 * whole_x / c),0)
		Redu.append(JSD_i)
	return Redu

def GeneRedMatrix(x,m=10):
	xt = x.T + eps
	xt = xt/np.sum(xt,0)
	Par = Partition(xt.shape[1],xt,m)
	pool = multiprocessing.Pool()
	out=pool.map(ProcessPartition,Par)
	pool.close();pool.join()
	R = out[0]
	for i in range(len(out)-1):
		R = np.vstack((R,out[i+1]))
	return R   
