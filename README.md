# Work in progress
The source code for diffDoamin will be uploaded by Nov 12, 2021.

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

# Installation instructions
diffDomain is tested on MacOS & Linux (Centos). 
## Dependences
diffDomain is dependent on xx. 

## Installation
Download diffDomain source package by running following command in a terminal:

`<git clone https://github.com/Tian-Dechao/diffDomain.git>`  

# Get started with example usage
We downloaded data [GEO:GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) from Rao et al [(2014)](https://www.sciencedirect.com/science/article/pii/S0092867414014974) for standalone example usage of diffDomain.
Example data saved in `<data/>`:
1. GM12878 TADs. 
2. GM12878 combined Hi-C data on Chr1 that is extracted by [Juicerbox](https://github.com/aidenlab/Juicebox) with resolution at 10 kb and normalization method at KR. The produced Hi-C data is 3-column: column 1 and column 2 are chromosomal bins, column 3 is KR normalized contact frequencies between the two bins.
3. K562 combined Hi-C data on Chr1. Settings are the same as GM12878.

## Testing if one TAD is reorganized
In this example, we tested the TAD (ChrX:xx-xx, use one GM12878 TAD that is clearly reorganized in K562). 
Data are saved in `<data/single-TAD/>`
Running the command 

`<python3 xxx>`

diffDomain also provide visualization function to visualize Hi-C matrices side-by-side.

`<python3 visualization step>`

Note: in this example, there is no need to do multiple comparison adjustment. 
Multiple comparisons adjustment by *BH* will be demonstrated in the next example. 

### Identifying the reorganized TADs on a 50 Mb region (Chr1:1-50,000,000)
In this example, multiple comparison adjustment is requried to adjust the *P*-values.
Data are saved in `<data/TADs_chr1/>`

`<python3 compare xx xx xx>`

Adjusting multiple comparisons by *BH* method

`<python3 adjustment BH >`

Visualization of GM12878 TADs that are reorganized in K562

`<python3 xx xxx >`

For interactive integrative analysis, we recommend using the [Nucleome Browser](http://www.nucleome.org/).
Example visualization outputs are shown below. 

![reorganized TADs on chr1](/figures/TADs_chr1.jpg)

### Identifying GM12878 TADs that are reorganized in  K562, using all TADs.
Data is using Amazon.

Identifying reorganized TADs on each chromosome.


Combining results into one and adjusting for multiple comparisons 

`<python3 xxx>`

Filtering out reorganized TADs with *BH < 0.05*

`<python 3 xxx>`

The output is saved to `<res/reorganized_TADs_GM12878_K562.tsv>`

# Contact information
More information please contact Dunming Hua at huadm@mail2.sysu.edu.cn or Dechao Tian at tiandch@mail.sysu.edu.cn.
