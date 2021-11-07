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
## Testing if one TAD is reorganized
In this example, there is no need to do multiple comparison adjustment.

### Identifying the reorganized TADs from a user given list of TADs
In this example, multiple comparison adjustment is requried to adjust the *P*-values.


# Contact information
More information please contact Dechao Tian at tiandch@mail.sysu.edu.cn
