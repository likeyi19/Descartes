# Descartes

#### [DEtection of Spatial Chromatin Accessibility patteRns with inTEr-cellular correlationS]

abstract

![image]()

## Installation  

### Environment setup

1. We recommend you to build a python virtual environment with [Anaconda](https://docs.anaconda.com/free/anaconda/install/linux/).  If Anaconda (or miniconda) is already installed with Python3, skip to 2.

2. Create and activate a new virtual environment:

```
$ conda create -n descartes python=3.8
$ conda activate descartes
```

### Package installation

Python packages required by Cofea are listed below:

```
1. Python 3.8.18
2. Packages for Descartes and tutorial
  anndata==0.8.0
  episcanpy==0.3.2
  h5py==3.7.0
  hdf5storage==0.1.18
  jupyter-contrib-core==0.4.0
  loess==2.1.2
  louvain==0.7.1
  matplotlib==3.5.2
  matplotlib-inline==0.1.3
  memory-profiler==0.61.0
  numba==0.55.2
  numpy==1.22.4
  pandas==1.4.3
  patsy==0.5.2
  progress==1.6
  scanpy==1.9.1
  scikit-learn==1.1.1
  scipy==1.8.1
  seaborn==0.11.2
  statsmodels==0.13.2
  umap-learn==0.5.3
```

Install the package and other requirements:

```  
Package installation:
$ git clone https://github.com/likeyi19/Descartes   
$ cd Descartes   
$ pip install -r requirements.txt
```

## Tutorial

### Demo

We provide a [quick-start notebook](https://github.com/likeyi19/Descartes/blob/main/code/demo.ipynb) which describes the fundamentals in detail and reproduces the results of Cofea.

### Descartes

Sixteen parameters are necessary, including the path of dataset, the save path for results, the chosen number of peaks, the random seed, the TF-IDF computation method, the number of principal components (PC), the quantity of K means, the similarity calculation method, the iteration count, the spatial neighborhood selection approach, the number of neighbors, the peak filtering method, the quantity of peak filtering, the distance calculation method, and the data synthesis ratio.

For exsample:
```
$ cd code/
$ python descartes.py -fp ../data/scanpy.h5ad -sp ../result -n 10000 -sb 1 -pc 10 -k 20 -iter 4 -nb 5 -r 0.4
$ cd ..
```

Or you can get help in this way:
```  
$ python code/descartes.py -h
usage: descartes.py [-h] 
                [-fp file_path] 
                [-sp save_path] 
                [-n PC] 
                [-c CORR] 
                [-n SELECT_NUMBER] 
                [-s SEED_BASE]

optional arguments:
  -h, --help           show this help message and exit
  -l, --load_path      str, default=None
                       storage path of the h5ad file, which contains the peak-by-cell data to be processed.
  -t, --TFIDF          str, default='tfidf2', options=['tfidf1', 'tfidf2', 'tfidf3']
                       TF-IDF implementation, 'tfidf1' represents the original version of TF-IDF transformation, 'tfidf2' indicates the TF-IDF transformation used by Signac, and 'tfidf3' represents the TF-IDF transformation used by scOpen.
  -p, --PC             int, default=100
                       Dimension of cell-wise PCA.
  -c, --corr           str, default='PCC', options=['PCC', 'CSC', 'SPCC']
                       Correlation coefficient calculation method. PCC represents Pearson correlation coefficient, CSC indicates Cosine correlation coefficient and SPCC represents Spearman correlation coefficient.
  -n, --select_number  int, default=20000
                       Number of selected features.
  -s, --seed_base      int, default=2
                       Random seed.
```  

### Analysis

From the perspective of revealing biological insights, we evaluated Cofea via several numeric experiments including cell type-specific peaks annotation and candidate enhancers identification.

![image]()

We provide a [notebook]() to analyze the sample dataset.

## Contact 
If you have any questions, you can contact me from the email: <lky23@mails.tsinghua.edu.cn>