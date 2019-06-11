# -*- coding: utf-8 -*-
"""
Created on Mon Oct 22 20:12:03 2018

@author: Zhanying Feng
"""
import numpy as np
import scTIM

if __name__ == "__main__":
    file_name = '/home/fengzhanying/Project2/software/scTIM/Package/data.txt'
    alpha = 0.1;beta = 0.4;gamma = 0.5;

    data,gene = scTIM.PreProcess(file_name,'y')

    p = scTIM.CellRedMatrix(data)

    fs = scTIM.GeneSpecificity(data)

    red = scTIM.GeneRedMatrix(data)
    np.savetxt('R.txt',red)

    w = []
    for i in range(10):
        w.append(scTIM.ExtractGene(data,p,red,alpha,beta,gamma))
    w = np.array(w)

    s = (np.sum(w,0) == 10).astype('int')
    g = open('LandmarkGene.txt','w')
    for i in range(len(gene)):
        if s[i] == 10:
            g.write(gene[i] + '\n')
    g.close()
