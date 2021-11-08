"""
Test if one domain is significantly rewried in another biological condition
Usage:
    scriptname dvsd one <chr> <start> <end> <hic0> <hic1> [options]
    scriptname dvsd multiple <hic0> <hic1> <bed> [options]
    scriptname visualization <chr> <start> <end> <hic0> <hic1> [options]
    scriptname adjustment <input> <output> [options]
    scriptname adjustment <input> <output> [options]

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
    --adj MultiComparison method  [default: bh]
    --alpha  [default: 0.05]


"""

import matplotlib
matplotlib.use('Agg')
import os
import sys
import pandas as pd
import numpy as np
from docopt import docopt
from  multiprocessing import Pool
from comparedomains import comp2domins_by_twtest, loadtads, visualization

opts = __doc__
scriptname = os.path.basename(__file__)
print(scriptname)
opts = opts.replace("scriptname", scriptname)
opts = docopt(doc=opts)
# save opts to output for every output
opts_df = pd.DataFrame.from_dict(opts, orient='index')
opts_df.index = '#' + opts_df.index

if(opts['dvsd']):
    if opts['one']:
        result = comp2domins_by_twtest(chrn=opts['<chr>'], start=int(opts['<start>']), end=int(opts['<end>']), reso=int(opts['--reso']), hicnorm=opts['--hicnorm'], fhic0=opts['<hic0>'], fhic1=opts['<hic1>'], min_nbin=int(opts['--min_nbin']), f=opts['--f'])
        result = pd.DataFrame(result)
        print(result)
        opts_df.to_csv(opts['--ofile'], sep='\t', index=True, header=False)
        result.to_csv(opts['--ofile'], sep='\t', index=False, header=False, mode='a')

    if opts['multiple']:
        reso = int(opts['--reso'])
        tadb = loadtads(opts['<bed>'], sep=opts['--sep'], chrnum=opts['--chrn'], min_nbin=int(opts['--min_nbin']), reso=int(opts['--reso']))
        # Parallel implementation
        def comp2domins_by_twtest_parallel(i):
            tmp_res = comp2domins_by_twtest(chrn=tadb.iloc[i, 0], start=tadb.iloc[i, 1], end=tadb.iloc[i, 2], reso=int(opts['--reso']), hicnorm=opts['--hicnorm'], fhic0=opts['<hic0>'], fhic1=opts['<hic1>'],min_nbin=int(opts['--min_nbin']),f=opts['--f'])
            return tmp_res

        P = Pool(int(opts['--ncore']))
        result = P.map(comp2domins_by_twtest_parallel, range(tadb.shape[0]))
        P.terminate()
        result = pd.DataFrame(result)

        # first the options; then the result
        opts_df.to_csv(opts['--ofile'], sep='\t', index=True, header=False)
        result.to_csv(opts['--ofile'], sep='\t', index=False, header=False, mode='a')

if(opts['visualization']):
    def save_fig(fig_id, tight_layout=True, fig_extension="pdf", resolution=600):
        path = os.path.join("res/images", fig_id + "." + fig_extension)
        print("Saving figure", fig_id)
        if tight_layout:
            plt.tight_layout()
            plt.savefig(path, format=fig_extension, dpi=resolution)
    def plot(mat,outputfile):
        fig,ax = plt.subplots(figsize=(6,6))
        sns.heatmap(data=mat,vmax=int(round(np.nanmax(mat)))/5,cmap=cdict,cbar=False,mask=mat<0.1,square=True)
        plt.xticks([])
        plt.yticks([])
        ax.tick_params(bottom=False, top=False, left=False,right=False)
        save_fig(outputfile)

    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    cdict=LinearSegmentedColormap.from_list('mycmap',['w','#FF0000'])
    result = visualization(chrn=opts['<chr>'], start=int(opts['<start>']), end=int(opts['<end>']), reso=int(opts['--reso']), hicnorm=opts['--hicnorm'], fhic0=opts['<hic0>'], fhic1=opts['<hic1>'])
    #the lower triangle is the first input hic0, the upper triangle is the second input hic1
    matrix1 = result[0]; matrix2 = result[1]
    matrix1[np.isnan(matrix1)] = 1; matrix2[np.isnan(matrix2)] = 1
    tril_mat = (np.tril(matrix1,-1)); triu_mat = (np.triu(matrix2,1))
    mat = tril_mat + triu_mat
    plot(mat,opts['--ofile'])


