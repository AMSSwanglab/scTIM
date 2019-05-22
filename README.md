# scTIM
## Introduction
A convenient tool for marker detection based on single cell RNA-seq data. <br>

## Package usage example:
First open a python IDE:
```bash
cd Package
python
```
Then run the following python script
```python
>>> import numpy as np
>>> import scTIM

>>> file_name = 'scTIM/Package/data.txt' ### Defining file name
>>> alpha = 0.1;beta = 0.4;gamma = 0.5;                                       ### Setting Parameters
>>> data,gene = scTIM.PreProcess(file_name,'y')                               ### Preprocessing data
>>> p = scTIM.CellRedMatrix(data)                                             ### Computing cell-cell distance matrix
>>> fs = scTIM.GeneSpecificity(data)                                          ### Computing gene specificity
>>> red = scTIM.GeneRedMatrix(data)                                           ### Computing gene-gene redundancy matrix
>>> w = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)                        ### Identifying markers by simulating annealing
>>> marker = [gene[i] for i in range(data.shape[0]) if w[i] == 1]             ### Output the marker set
```
For more robust solution, we repeat the simulating annealing for 10 times and use the inersection of 10 outcomes as final result and these 10 repeats can be conducted by parallel computing:
```python
>>> w1 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w2 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w3 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w4 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w5 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w6 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w7 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w8 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w9 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma)
>>> w10 = scTIM.ExtractGene(data,p,red,alpha,beta,gamma) 
>>> w = (np.sum([w1,w2,w3,w4,w5,w6,w7,w8,w9,w10],0)==10)                       ### Intersection
>>> marker = [gene[i] for i in range(data.shape[0]) if w[i] == 1]              ### Output the marker set
```
## Requirements:
Python environment: python 3 <br>
numpy <br>
Memory >= 3.0 Gb
