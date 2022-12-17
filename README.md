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
- matplotlib 
- statsmodels
- h5py 
- seaborn

Note: You should have installed libgcc.

## Installation

You can choose one of following methods. 

### to install python3 version from pypi

```
pip install diffDomain-py3
```

### to clone this repository

First of all, you should have a packages manager, such as [conda](https://docs.conda.io/en/latest/miniconda.html). Then you should create a new independent environment for diffDomain, and download diffDomain source package by running following command in a terminal:

```
git clone https://github.com/Tian-Dechao/diffDomain.git
cd diffDomain
pip install requirements.txt

```  

## Questions
If you encounter the following question, please don't  be too worried.
- **AttributeError: 'function' object has no attribute 'straw'** :
You can open the \_\_init\_\_.py of straw ( its pathway will be reported in the error, for example "/home/gum/.conda/envs/diffDomain/lib/python2.7/site-packages/straw/__ init_.py" ) and then deleted the sentence “straw = straw_module.straw”

  
# Usage

## Main method
**Usage:**  
    python diffdomains.py dvsd one \<chr> \<start> \<end> \<hic0> \<hic1> [options]  
    python diffdomains.py dvsd multiple \<hic0> \<hic1> \<tadlist_of_condition1.bed> [options]  
    python diffdomains.py visualization \<chr> \<start> \<end> \<hic0> \<hic1> [options]  
    python diffdomains.py adjustment \<method> \<input> \<output> [options] 

**Options:**  
    --reso resolution for hicfile  [default: 100000]   
    --min_nbin effective number of bin  [default: 10]    
    --f parameter for filtering the null values of the matrix[0~1)  [default: 0.5] For example, when setting ‘--f 0.6’, in the contact matrix of a TAD, if the number of the columns, whose proportions of missing values is higher than 40%, is smaller than min_nbin, DiffDomain will skip comparing this TAD anymore and set its result (statistics, the 5th column ; P value, the 6th column) as NAN.   
    --ofile filepath for output file  [default: stdout]  
    --oprefix prefix for output files  
    --oprefixFig prefix for output figures  
    --sep deliminator for hicfile  [default: \t]  
    --hicnorm hic matrix normalization method  [default: KR]  
    --chrn chromosome number  [default: ALL]   
    --ncore number of parallel process  [default: 10]  
    --filter As long as the pvalue of TADs is less than 0.05 after adjustment if argument is true  [default: false]  
  
**Note:**
  1. for most of the bulk Hi-C data, such as hic data in Adiden [Reference], results is not sensitive to the exact value of --f.  
  2. For single-cell Hi-C data, recommend users try multiple values of --f and choose one with acceptable number of TADs compared. Due to high sparisity in single-cell Hi-C data and variation in imputation methods (such as scHiCluster, Higashi, scVI-3D), we did not set a default value of --f.  

## Classification
**Usage:**  
    python classification.py -d \<result_of_diffdomains.py_multiple> -t \ <tadlist_of_condition2> [options]

**options**    
    --limit length of bases, within which the boundaries will be judged as common boundaries [default: 30000]  
    --out the filename of output [default: name_of_-d_types.txt] .  
    --kpercent the common boundareis are within max(l*bin,k% TAD's length) [default: 10] . 
    --remote the limitation of the biggeset region [default: 1000000]  
    --s1 int, to skip the first s1 rows in -d [default: 0]  
    --s2 int, to skip the first s2 rows in -t [default: 0]   
    --sep1 the separater of -d [default: \t]    
    --sep2 the separater of -t [default: \t]   
    
Note:   
You can set the --limit to adjust the 'common boundary'.  
As said in paper,we use '3bin' as the filter of common boundaies.  
That means if we use the 10kb resolution, we will set --limit as 30000, and if 25kb resolution, --limit will be 75000.


  

# Contact information

More information please contact Dunming Hua at huadm@mail2.sysu.edu.cn, Ming Gu at guming5@mail2.sysu.edu.cn or Dechao Tian at tiandch@mail.sysu.edu.cn.
