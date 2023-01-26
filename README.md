# diffDomain 

## A short description

diffDomain is a new computational method for identifying reorganized TADs using chromatin contact maps from two biological conditions.   
  
## A long description diffDomain

The workflow of diffDomain is illustrated down below.  

The goal is to test if a TAD identified in one biological condition has structural changes in another biological condition.  

The core of diffDomain is formulating the problem as a hypothesis testing problem where the null hypothesis is that the TAD doesn't undergo significant structural reorganization at later condition.
The input are Hi-C contact matrices of the TAD region in the two biological conditions (*A*).
The Hi-C contact matrices are  log-transformed to adjust for the exponential decay of Hi-C contacts between chromosome bins with increased distances.  

Their entry-wise difference is calculated (*B*).  

The difference matrix *D* is normalized by iteratively standardizing its *k*-off diagonal parts, *-N+2 <= k <= N-2*, adjusting absolute differences in contact frequencies due to different sequencing depths in the two biological conditions (*C*).  

Note that, standardization is TAD-specific. Each TAD has its own parameters that are only estimated from its contact matrices in a pair of biological conditions.  

Intuitively, if a TAD is not significantly reorganized, normalized *D* would resemble a random matrix with white noise entries, enabling us to borrow theoretical results in random matrix theory.
Indeed, normalized *D* is a generalized Wigner matrix (D), a well studied high-dimensional random matrices.  

Its largest singular value is proved to be fluctuating around 2 under the null hypothesis.
Armed with the fact, diffDomain reformulates the reorganized TAD identification problem into a hypothesis testing problem:  
1. H0: the largest singular value equals to 2;  
2. H1: the largest singular value is greater than  2.  

For a user given set of TADs, *P* values are adjusted for multiple comparisons using *BH* method as default.  
Once we identify the subset of reorganized TADs, we classify them into six subtypes to aid biological analysis and interpretations (F).  
A few examples of reorganized TADs identified by diffDomain in two datasets are shown in (G).  


![workflow](/figures/workflow.jpg)

## Installation instructions

diffDomain is tested on MacOS & Linux (Centos).   

## Dependences

diffDomain-py2 is dependent on 
- Python 2.7
- hic-straw==0.0.6 

diffDomain-py3 is dependent on
- Python 3
- hic-straw==1.3.1

and
- TracyWidom 
- pandas 
- numpy 
- docopt 
- tqdm
- matplotlib 
- statsmodels
- h5py 
- seaborn


## Installation

First of all, we recommend you to have a packages manager, such as [conda](https://docs.conda.io/en/latest/miniconda.html), and create a new independent environment for diffDomain.

Note: If you encounter error when installing hicstraw that diffDomain relies on. You can use [conda](https://docs.conda.io/en/latest/miniconda.html) to install it easily:

```
conda install -c bioconda hic-straw
```

Then you can choose one of following methods. 

### to install python3 version from pypi

```
pip install diffDomain-py3
```

### to clone this repository


```
git clone https://github.com/Tian-Dechao/diffDomain.git
cd diffDomain
pip install requirements.txt

```  

## Documentation
Please see the [wiki](https://github.com/Tian-Dechao/diffDomain/wiki/0.Usage) for extensive documentation and example tutorials.

# Contact information

More information please contact Dunming Hua at huadm@mail2.sysu.edu.cn, Ming Gu at guming5@mail2.sysu.edu.cn or Dechao Tian at tiandch@mail.sysu.edu.cn.
