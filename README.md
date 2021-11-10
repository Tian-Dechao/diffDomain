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
diffDomain is dependent on 
- Python 2.7
- hic-straw 0.0.6
- TracyWidom 0.3.0
- pandas 0.18.1
- numpy 1.15.0
- docopt 0.6.2
- matplotlib 1.5.1
- statsmodels 0.6.1

## Installation
Download diffDomain source package by running following command in a terminal:

```
git clone https://github.com/Tian-Dechao/diffDomain.git
```  

# Get started with example usage
We downloaded data [GEO:GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) from Rao et al [(2014)](https://www.sciencedirect.com/science/article/pii/S0092867414014974) for standalone example usage of diffDomain.
Example data saved in `<data/>`:
1. GM12878 TADs. 
2. GM12878 combined Hi-C data on Chr1 that is extracted by [Juicerbox](https://github.com/aidenlab/Juicebox) with resolution at 10 kb and normalization method at KR. The produced Hi-C data is 3-column: column 1 and column 2 are chromosomal bins, column 3 is KR normalized contact frequencies between the two bins.
3. K562 combined Hi-C data on Chr1. Settings are the same as GM12878.

## Testing if one TAD is reorganized
In this example, we tested the GM12878 TAD that is reorganized in K562 (Chr1:163500000-165000000, [Ref](http://dx.doi.org/10.1016/j.molcel.2017.07.022)). 
Data are saved in `<data/single-TAD/>`.

Running the command 

- Usage: scriptname dvsd one \<chr> \<start> \<end> \<hic0> \<hic1> [options]

```
python src/diffdomains.py dvsd one 1 163500000 165000000 data/single-TAD/GM12878_chr1_163500000_165000000_res_10k.txt data/single-TAD/K562_chr1_163500000_165000000_res_10k.txt --reso 10000 --ofile res/chr1_163500000_165000000.txt
```

diffDomain also provide visualization function to visualize Hi-C matrices side-by-side.

- Usage: scriptname visualization \<chr> \<start> \<end> \<hic0> \<hic1> [options]

Figure are saved in `<res/images/>`.

```
python src/diffdomains.py visualization 1 163500000 165000000 data/single-TAD/GM12878_chr1_163500000_165000000_res_10k.txt data/single-TAD/K562_chr1_163500000_165000000_res_10k.txt --reso 10000 --ofile res/images/side_by_side
```
Note: in this example, there is no need to do multiple comparison adjustment. 
Multiple comparisons adjustment by *BH* will be demonstrated in the next example. 

### Identifying the reorganized TADs in chr1
In this example, multiple comparison adjustment is requried to adjust the *P*-values.
chr1_domainlist are saved in `<data/TADs_chr1/>`, [GM12878.hic](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FGM12878%5Finsitu%5Fprimary%2Breplicate%5Fcombined%2Ehic) and [K562.hic](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525%5FK562%5Fcombined%2Ehic) comes from these two links, because the file is too large, you may need to download it to your hard drive in advance.

- Usage: scriptname dvsd multiple \<hic0> \<hic1> \<bed> [options]

```
python src/diffdomains.py dvsd multiple GM12878.hic K562.hic data/TADs_chr1/GM12878_chr1_domainlist.txt --ofile res/temp/GSE63525_GM12878_vs_K562_reso_10k_chr1.txt
```

- Adjusting multiple comparisons by *BH* method (default, Optional parameters: *fdr_by*, *bonferroni*, *holm*, *hommel* etc.)
- Usage: scriptname adjustment \<method> \<input> \<output> 

```
python src/diffdomains.py adjustment fdr_bh res/temp/GSE63525_GM12878_vs_K562_reso_10k_chr1.txt res/adjusted_chr1.tsv --filter true
```

For interactive integrative analysis, we recommend using the [Nucleome Browser](http://www.nucleome.org/).
Example visualization outputs are shown below. 

![reorganized TADs on chr1](/figures/TADs_chr1.jpg)

### Identifying GM12878 TADs that are reorganized in  K562, using all TADs.
Data is using Amazon.

- Identify TADs in multiple chromosomes simultaneously. 

```
python src/diffdomains.py dvsd multiple https://hicfiles.s3.amazonaws.com/hiseq/gm12878/in-situ/combined.hic https://hicfiles.s3.amazonaws.com/hiseq/k562/in-situ/combined.hic data/GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt --ofile res/temp/temp.txt
```

- MultiComparison adjustment.

```
python src/diffdomains.py adjustment fdr_bh res/temp/temp.txt res/adjusted_TADs2.txt 
```

- optional parameter **[--filter]**, Filtering out reorganized TADs with *BH < 0.05*.

```
python src/diffdomains.py adjustment fdr_bh res/temp/temp.txt res/reorganized_TADs_GM12878_K562.tsv --filter true
```

The final output is saved to `<res/reorganized_TADs_GM12878_K562.tsv>`.

# Contact information
More information please contact Dunming Hua at huadm@mail2.sysu.edu.cn or Dechao Tian at tiandch@mail.sysu.edu.cn.
