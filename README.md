# scTIM
## Introduction
A convenient tool for marker detection based on single cell RNA-seq data. <br>

## Installation

```bash
pip install sc_tim
```

## Usage example:

Then run the following python script
```python
import numpy as np
import sc_tim

if __name__ == "__main__":                                                  ### This command is necessary for Windows System
  file_name = 'scTIM-master/Package/data.txt'                               ### Defining file name
  alpha = 0.1;beta = 0.4;gamma = 0.5;                                       ### Setting Parameters
  data,gene = sc_tim.PreProcess(file_name,'y')                              ### Preprocessing data
  p = sc_tim.CellRedMatrix(data)                                            ### Computing cell-cell distance matrix
  fs = sc_tim.GeneSpecificity(data)                                         ### Computing gene specificity
  red = sc_tim.GeneRedMatrix(data)                                          ### Computing gene-gene redundancy matrix
  w = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)                    ### Identifying markers by simulating annealing
  marker = [gene[i] for i in range(data.shape[0]) if w[i] == 1]             ### Output the marker set
```
For more robust solution, we repeat the simulating annealing for 10 times and use the inersection of 10 outcomes as final result and these 10 repeats can be conducted by parallel computing:
```python
w1 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w2 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w3 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w4 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w5 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w6 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w7 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w8 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w9 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma)
w10 = sc_tim.ExtractGene(data,fs,p,red,alpha,beta,gamma) 
w = (np.sum([w1,w2,w3,w4,w5,w6,w7,w8,w9,w10],0)==10)                       ### Intersection
marker = [gene[i] for i in range(data.shape[0]) if w[i] == 1]              ### Output the marker set
```

## Requirements
Operating system: Linux (strongly recommended but not necessary) <br>
Python environment: python 3 <br>
Python package: numpy <br>
Memory: >= 3.0 Gb

## Citation
If you use scTIM or scTIM associated concepts, please cite

[Zhanying, et al. scTIM: seeking cell-type-indicative marker from single cell RNA-seq data by consensus optimization. Bioinformatics, 2020.] (https://academic.oup.com/bioinformatics/article-abstract/36/8/2474/5679774?redirectedFrom=fulltext)
