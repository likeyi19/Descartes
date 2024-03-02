# Descartes

We present Descartes, a graph-based model, for DEtection of Spatial Chromatin Accessibility patteRns with inTEr-cellular correlationS. Leveraging the graph of inter-cellular correlations, Descartes adeptly evaluates and identifies SV peaks by analyzing the self-correlations of peaks within the graph. To navigate the inherent challenge of highly dispersed accessibility patterns in spATAC-seq data, Descartes incorporate chromatin accessibility information with spatial locations during graph construction, and iteratively updates the graph to capture the intricate relationships between neighboring cells.

<p align="center">
  <img src="https://github.com/likeyi19/Descartes/blob/main/inst/model.png" width="600" height="450" alt="image">
</p>


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
  anndata >= 0.9.2
  matplotlib >= 3.7.4
  numpy >= 1.22.4
  pandas >= 1.4.3
  scanpy == 1.9.6
  scikit-learn >= 1.3.0
  scipy >= 1.8.0
  seaborn >= 0.12.2
```

Install the package and other requirements:

```  
Package installation:
$ git clone https://github.com/likeyi19/Descartes   
$ cd Descartes   
$ pip install -r requirements.txt
```

Install descartes:

```
$ pip install decare
```

### Sample dataset download

We provide a slice of mouse brain as a sample dataset, which can be downloaded at: [https://cloud.tsinghua.edu.cn/d/71c9840593464e9fa122/](https://cloud.tsinghua.edu.cn/d/71c9840593464e9fa122/)

## Tutorial

We provide four Jupyter notebooks to demonstrate the functions of our methods, including [SV peaks selection](https://github.com/likeyi19/Descartes/blob/main/code/SV%20peaks%20selection.ipynb), [Data imputation](https://github.com/likeyi19/Descartes/blob/main/code/Data%20imputation.ipynb), [Peak module identification](https://github.com/likeyi19/Descartes/blob/main/code/Peak%20module%20identification.ipynb), and [Gene-peak interaction detection](https://github.com/likeyi19/Descartes/blob/main/code/Gene-peak%20interaction%20detection.ipynb).

### Descartes

Sixteen parameters are necessary, including the path of dataset, the save path for results, the chosen number of peaks, the random seed, the TF-IDF computation method, the number of principal components (PC), the quantity of K means, the similarity calculation method, the iteration count, the spatial neighborhood selection approach, the number of neighbors, the spatial strategy for score calculation, the peak filtering method, the quantity of peak filtering, the distance calculation method, and the data synthesis ratio.

For exsample:
```
$ cd code/
$ python descartes.py -fp ../data/scanpy.h5ad -sp ../result -n 10000 -sb 1 -pc 10 -k 20 -iter 4 -nb 5 -r 0.4
$ cd ..
```

Or you can get help in this way:
```  
$ python code/descartes.py -h
usage: descartes.py [-h] [-fp FILE_PATH] [-sp SAVE_PATH] [-n NUM_SELECT_PEAK]
                    [-sb SEED_BASE] [-tf TF_IDF] [-pc PC_NUMBER] [-k K_NUMBER]
                    [-s SIMILARITY] [-iter ITER_TIME] [-spm SP_METHOD]
                    [-nb NEIGHBOR] [-spd SP_DIST] [-ps PRE_SELECT]
                    [-pn PEAKS_NUM] [-d DISTANCE] [-r RATIO]

optional arguments:
  -h, --help            show this help message and exit
  -fp FILE_PATH, --file_path FILE_PATH
                        The path of dataset
  -sp SAVE_PATH, --save_path SAVE_PATH
                        The save path for results
  -n NUM_SELECT_PEAK, --num_select_peak NUM_SELECT_PEAK
                        The chosen number of peaks, defaults to 10000
  -sb SEED_BASE, --seed_base SEED_BASE
                        The random seed
  -tf TF_IDF, --TF_IDF TF_IDF
                        The TF-IDF computation method
  -pc PC_NUMBER, --pc_number PC_NUMBER
                        The number of principal components
  -k K_NUMBER, --k_number K_NUMBER
                        The quantity of K means
  -s SIMILARITY, --similarity SIMILARITY
                        The similarity calculation method
  -iter ITER_TIME, --iter_time ITER_TIME
                        The iteration count, defaults to 4
  -spm SP_METHOD, --sp_method SP_METHOD
                        The spatial neighborhood selection approach
  -nb NEIGHBOR, --neighbor NEIGHBOR
                        The number of neighbors
  -spd SP_DIST, --sp_dist SP_DIST
                        The spatial strategy for score calculation
  -ps PRE_SELECT, --pre_select PRE_SELECT
                        Peak filtering method
  -pn PEAKS_NUM, --peaks_num PEAKS_NUM
                        The quantity of peak filtering
  -d DISTANCE, --distance DISTANCE
                        The distance calculation method
  -r RATIO, --ratio RATIO
                        Data synthesis ratio
```  

## Contact 
If you have any questions, you can contact me from the email: <lky23@mails.tsinghua.edu.cn>
