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


diffDomain is dependent on 
- Python 2.7
- hic-straw 0.0.6
- TracyWidom 0.3.0
- pandas 0.18.1
- numpy 1.15.0
- docopt 0.6.2
- matplotlib 1.5.1
- statsmodels 0.6.1
- h5py 2.9.0
- seaborn 0.9.0

## Installation

You can choose one of following three methods. 

### Method 1

First of all, you should have a packages manager, such as [conda](https://docs.conda.io/en/latest/miniconda.html). Then you should create a new independent environment for diffDomain, and download diffDomain source package by running following command in a terminal:

```
conda create -n diffDomain # the name of the new environment
conda activate diffDomain # to activate the new environment
conda install python=2.7 
# in this step, pip 20.1.1 should be installed automatically

git clone https://github.com/Tian-Dechao/diffDomain.git
cd diffDomain
pip install requirements.txt

```  

### Method 2
```
conda create -n diffDomain # the name of the new environment
conda activate diffDomain # to activate the new environment
conda install python=2.7 
# in this step, pip 20.1.1 should be installed automatically
pip install diffDomain
# in this steps, requirements will be installed automatically
```
### Method 3
You can use diffDomain in a Docker's container by pulling diffDomain's image based on centos7: 
```
docker pull guming5/centos7:diffDomain
```
## Questions
If you don't use the 3rd installation method, you will encounter the following question when you import diffdomain. Please don't  be too worried.
- **AttributeError: 'function' object has no attribute 'straw'** :
You can open the \_\_init\_\_.py of straw ( its pathway will be reported in the error, for example "/home/gum/.conda/envs/diffDomain/lib/python2.7/site-packages/straw/__ init_.py" ) and then deleted the sentence “straw = straw_module.straw”

## Get started with example usage
  
Note:  
If you download diffDomain by github, the flowing "diffdomain/……" pathway means "src/……"   
  
We downloaded data [GEO:GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) from Rao et al [(2014)](https://www.sciencedirect.com/science/article/pii/S0092867414014974) for standalone example usage of diffDomain.
Example data saved in `<data/>`:
1. GM12878 TADs. 
2. GM12878 combined Hi-C data on Chr1 that is extracted by [Juicebox](https://github.com/aidenlab/Juicebox) with resolution at 10 kb and normalization method at KR. The produced Hi-C data is 3-column: column 1 and column 2 are chromosomal bins, column 3 is KR normalized contact frequencies between the two bins.
3. K562 combined Hi-C data on Chr1. Settings are the same as GM12878.

### Testing if one TAD is reorganized

In this example, we tested the GM12878 TAD that is reorganized in K562 (Chr1:163500000-165000000, [Ref](http://dx.doi.org/10.1016/j.molcel.2017.07.022)). 
Data are saved in `<data/single-TAD/>`.

Running the command 

- Usage: scriptname dvsd one \<chr> \<start> \<end> \<hic0> \<hic1> [options]

```
python diffdomain/diffdomains.py dvsd one 1 163500000 165000000 data/single-TAD/GM12878_chr1_163500000_165000000_res_10k.txt data/single-TAD/K562_chr1_163500000_165000000_res_10k.txt --reso 10000 --ofile res/chr1_163500000_165000000.txt
```
In Python:

```
from diffdomain import pydiff   
# from diffdomain not diffDomain or DiffDomain
"""
 If you want to use diffDomain in Python,
 please use APIs in diffdomain.pydiff
"""
result = pydiff.diffdomain_one(chr = 1,start = 163500000,
    end = 165000000,
    fhic0 = 'data/single-TAD/GM12878_chr1_163500000_165000000_res_10k.txt', 
    fhic1 = 'data/single-TAD/K562_chr1_163500000_165000000_res_10k.txt',
    reso = 10000
    )
```
In Python, pydiff.diffdomain_one will return the dataframe of result.

diffDomain also provide visualization function to visualize Hi-C matrices side-by-side.

- Usage: scriptname visualization \<chr> \<start> \<end> \<hic0> \<hic1> [options]

Figure are saved in `<res/images/>`.

```
python diffdomain/diffdomains.py visualization 1 163500000 165000000 data/single-TAD/GM12878_chr1_163500000_165000000_res_10k.txt data/single-TAD/K562_chr1_163500000_165000000_res_10k.txt --reso 10000 --ofile res/images/side_by_side
```
In Python:
```
pydiff.visualization(chr=1,start=163500000,end=165000000,
    fhic0 = 'data/single-TAD/GM12878_chr1_163500000_165000000_res_10k.txt',
    fhic1 = 'data/single-TAD/K562_chr1_163500000_165000000_res_10k.txt',
    reso = 10000
    ofile = 'res/images/side_by_side'
    )
```

Note: in this example, there is no need to do multiple comparison adjustment. 
Multiple comparisons adjustment by *BH* will be demonstrated in the next example. 

### Identifying the reorganized TADs on a 50 Mb region (Chr1:1-50,000,000)

In this example, multiple comparison adjustment is requried to adjust the *P*-values.
chr1_50M_domainlist are saved in `<data/TADs_chr1/>`.

- Usage: scriptname dvsd multiple \<hic0> \<hic1> \<bed> [options]

```
python diffdomain/diffdomains.py dvsd multiple data/TADs_chr1/chr1_50M_GM12878.h5 data/TADs_chr1/chr1_50M_K562.h5 data/TADs_chr1/GM12878_chr1_50M_domainlist.txt --reso 10000 --ofile res/temp/GM12878_vs_K562_chr1_50M_temp.txt
```
In Python:

```
result_mul = pydiff.diffdomain_multiple(fhic0='data/TADs_chr1/chr1_50M_GM12878.h5',
    fhic1 = 'data/TADs_chr1/chr1_50M_K562.h5',
    fbed = 'data/TADs_chr1/GM12878_chr1_50M_domainlist.txt',
    reso = 10000
    )
```
The function pydiff.diffdomain_multiple will return the dataframe of result_mul.

- Adjusting multiple comparisons by *BH* method (default, Optional parameters: *fdr_by*, *bonferroni*, *holm*, *hommel* etc.) and Filtering out reorganized TADs with *BH < 0.05*
- Usage: scriptname adjustment \<method> \<input> \<output> 

```
python diffdomain/diffdomains.py adjustment fdr_bh res/temp/GM12878_vs_K562_chr1_50M_temp.txt res/GM12878_vs_K562_chr1_50M_adjusted_filter.tsv --filter true
```

For interactive integrative analysis, we recommend using the [Nucleome Browser](http://www.nucleome.org/).
Example visualization outputs are shown below. 

![reorganized TADs on chr1](/figures/TADs_chr1.png)

### Identifying GM12878 TADs that are reorganized in  K562, using all TADs.

Data is using Amazon.

- Identify TADs in multiple chromosomes simultaneously. 

```
python diffdomain/diffdomains.py dvsd multiple https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic data/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt --ofile res/temp/temp.txt
```
In Python:
```
result_mul = pydiff.diffdomain_multiple(fhic0='https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic',
    fhic1 = 'https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic',
    fbed = 'data/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt'
    )
```

- MultiComparison adjustment.

```
python diffdomain/diffdomains.py adjustment fdr_bh res/temp/GM12878_vs_K562_chr1_50M_temp.txt res/adjusted_TADs2.txt 
```
In Python:
```
result_adj = pydiff.adjustment(inputdf = result_mul,method='fdr_bh')
```
The function of pydiff.adjustment will return the dataframe of result_adj (adjusted result_mul by BH).

- optional parameter **[--filter]**, Filtering out reorganized TADs with *BH < 0.05*.

```
python diffdomain/diffdomains.py adjustment fdr_bh res/temp/GM12878_vs_K562_chr1_50M_temp.txt res/reorganized_TADs_GM12878_K562.tsv --filter true
```
The final output is saved to `<res/reorganized_TADs_GM12878_K562.tsv>`.

In Python:
```
result_adj_filter = pydiff.adjustment(inputdf = result_mul,method='fdr_bh',
    Filter=True
    )
```
In Python, the result will be returned as a dataframe.

- Classification of TADs

In this step, you will need the tadlist of the second hic file.

Running the command:

```
python diffdomain/classificattion.py -d adjusted_TADs2.txt -t GSE63525_K562_Arrowhead_domainlist.txt 
```
In Python:

```
"""
 You should ensure that you have not 
  taken a row of TAD as the header of the dataframe.
"""
tadlist = pd.read_table('data/GSE63525_K562_Arrowhead_domainlist.txt',header=None)

types = pydiff.classification(result_adj,tadlist)

```

# Summary
## Command line

Usage:
    scriptname dvsd one <chr> <start> <end> <hic0> <hic1> [options]
    scriptname dvsd multiple <hic0> <hic1> <bed> [options]
    scriptname visualization <chr> <start> <end> <hic0> <hic1> [options]
    scriptname adjustment <method> <input> <output> [options]

Options:
    --ofile filepath for output file  [default: stdout]
    --oprefix prefix for output files
    --oprefixFig prefix for output figures
    --sep deliminator for hicfile  [default: \t]
    --hicnorm hic matrix normalization method  [default: KR]
    --chrn chromosome number  [default: ALL]
    --reso resolution for hicfile  [default: 100000]
    --ncore number of parallel process  [default: 10]
    --min_nbin effective number of bin  [default: 10]
    --f parameters for filtering the null values of the matrix[0~1)  [default: 0.5]
    --filter As long as the pvalue of TADs is less than 0.05 after adjustment if argument is true  [default: false]
  
## Python
  
- diffdomain_one(<chr>,<start>,<end>,<hic0>,<hic1>,[options])
    
    --chr chromosome number
    --start the start position of the domain tested by this func
    --end the end position of the domain
    --fhic0 the path of the first hic file
    --fhic1 the path of the second hic file
    
    Opts:
    --reso resolution for hicfile  [default: 100000]
    --hicnorm hic matrix normalization method [default: 'KR']
    --min_nbin effective number of bin  [default: 10]
    --f parameters for filtering the null values of the matrix[0~1)  [default: 0.5]
 
- diffdomain_multiple(<fhic0>,<fhic1>,<fbed>,[options])

    --fhic0 the filepath of the first hic file
    --fhic1 the filepath of the second hic file
    --fbed the filepath of TADs' list that you want to test,it usually is the tadlist of hic0(the first hic file)

    Options:
    --sep deliminator for hicfile  [default: '\t']
    --hicnorm hic matrix normalization method  [default: 'KR']
    --chrn chromosome number  [default: 'ALL']
    --reso resolution for hicfile  [default: 100000]
    --ncore number of parallel process  [default: 10]
    --min_nbin effective number of bin  [default: 10]
    --f parameters for filtering the null values of the matrix[0~1)  [default: 0.5]

- adjustment(<inputdf>,<Filter>,[options])
    
    --inputdf the result of diffdomain_multiple (pd.DataFrame)
    
    Options:
    --alpha the threshold of adjusted pvalue [default: 0.05]
    --Filter As long as the pvalue of TADs is less than alpha after adjustment if argument is true [True/False, default: False]
    --method adjustment method you want to use [default: 'fdr_bh']
    --sort wheter to sort the result [default: False]
  
- visualizing(<chr>,<start>,<end>,<fhic0>,<fhic1>,[options])
    
    --chr chromosome number
    --start the start position of the domain visulized by the func
    --end the end position of the domain visulized by the func
    --fhic0 the filepath of the first hic file
    --fhic1 the filepath of the second hic file
    --ofile filepath for output file
   
    Options:
    --reso resolution for hicfile  [default: 100000]
    --hicnorm hic matrix normalization method [default: 'KR']
  
- classification(<result_adj_df>,<tadlist_df>,[options])
    --result_adj_df the dateframe of adjusted outcome of diffdomain_multiple
    --tadlist_df the tadlist(dataframe) of the second hic file
    
    Options:
    --alpha the threshold of adjusted p-value
    --limit length of bases, within which the boundaries will be judged as common boundaries [default: 40000]
    --kpercent the common boundareis are within max(l*bin,k% TAD's length) [default: 10]
    --remote the limitation of the biggeset region
    --ofile the filepath of output file
  
# Contact information

More information please contact Dunming Hua at huadm@mail2.sysu.edu.cn, Ming Gu at guming5@mail2.sysu.edu.cn or Dechao Tian at tiandch@mail.sysu.edu.cn.
