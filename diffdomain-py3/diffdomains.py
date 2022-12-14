"""
Test if one domain is significantly rewried in another biological condition
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


"""

import matplotlib
matplotlib.use('Agg')
import os
import sys
import pandas as pd
import numpy as np
from docopt import docopt
from  multiprocessing import Pool
from utils import comp2domins_by_twtest, loadtads, visualization

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
        result = pd.DataFrame(result).T
        print(result)
        opts_df.to_csv(opts['--ofile'], sep='\t', index=True, header=False)
        result.to_csv(opts['--ofile'], sep='\t', index=False, header=False, mode='a')

    if opts['multiple']:
        reso = int(opts['--reso'])
        tadb = loadtads(opts['<bed>'], sep=opts['--sep'], chrnum=opts['--chrn'], min_nbin=int(opts['--min_nbin']), reso=int(opts['--reso']))
        # Parallel implementation
        def comp2domins_by_twtest_parallel(i):
            tmp_res = comp2domins_by_twtest(chrn=tadb.iloc[i, 0], start=tadb.iloc[i, 1], 
            end=tadb.iloc[i, 2], reso=int(opts['--reso']), hicnorm=opts['--hicnorm'],
             fhic0=opts['<hic0>'], fhic1=opts['<hic1>'],min_nbin=int(opts['--min_nbin']),f=opts['--f'])
            return tmp_res
        if __name__=='__main__':
            try:
                # if the hicstraw can not find the chromosomes
                # it will directly end the script
                # without any exception !!
                comp2domins_by_twtest_parallel(0)
            finally:
                P = Pool(int(opts['--ncore']))
                note = []
                for i in range(tadb.shape[0]):
                    note.append(P.apply_async(comp2domins_by_twtest_parallel, [i]))
                P.close()
                P.join()
                result = []
                for i in note:
                    result.append(i.get())
                result = pd.DataFrame(result)

                # first the options; then the result
                opts_df.to_csv(opts['--ofile'], sep='\t', index=True, header=False)
                result.to_csv(opts['--ofile'], sep='\t', index=False, header=False, mode='a')

elif(opts['visualization']):
    #
    def save_fig(fig_id, tight_layout=True, fig_extension="pdf", resolution=600):
        path = os.path.join(fig_id + "." + fig_extension)
        print("Saving figure", fig_id)
        if tight_layout:
            plt.tight_layout()
            plt.savefig(path, format=fig_extension, dpi=resolution)

    def plot(mat,outputfile):
        fig,ax = plt.subplots(figsize=(6,6))
        sns.heatmap(data=mat,
                vmax=int(round(np.nanmax(mat)))/5,
                cmap=cdict,cbar=False,
                mask=mat<0.1,square=True)
        plt.xticks([])
        plt.yticks([])
        ax.tick_params(bottom=False, top=False, left=False,right=False)
        save_fig(outputfile)

    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import LinearSegmentedColormap
    cdict=LinearSegmentedColormap.from_list('mycmap',['w','#FF0000'])
    result = visualization(chrn=opts['<chr>'], start=int(opts['<start>']),
                           end=int(opts['<end>']), reso=int(opts['--reso']),
                           hicnorm=opts['--hicnorm'], fhic0=opts['<hic0>'],
                           fhic1=opts['<hic1>'])
    #the lower triangle is the first input hic0, the upper triangle is the second input hic1
    matrix1 = result[0] 
    matrix2 = result[1]
    matrix1[np.isnan(matrix1)] = 1
    matrix2[np.isnan(matrix2)] = 1
    tril_mat = (np.tril(matrix1,-1))
    triu_mat = (np.triu(matrix2,1))
    mat = tril_mat + triu_mat
    plot(mat,opts['--ofile'])


elif(opts['adjustment']):
    import statsmodels.stats.multitest as smm
    def read_compared(inputfile, method='fdr_bh', alpha=0.05, skiprows=25):
        df = pd.read_csv(inputfile, header=None, sep='\t', skiprows=skiprows, error_bad_lines=False)
        df.columns=['chr','start','end','region','stat','pvalue','bins']
        pvalue = df['pvalue']
        pvalue = df['pvalue'].fillna(pvalue.median())
        rej, pval_corr = smm.multipletests(pvalue, is_sorted=False, alpha=alpha, method=method)[:2]
        df['adj_pvalue'] = pval_corr
        return df


    res = read_compared(inputfile=opts['<input>'], method=opts['<method>'])
    if opts['--filter'] == "true":
        #save result
        res = res[res['adj_pvalue']<=0.05]
        print(res)
        res.to_csv(opts['<output>'], index=False, sep='\t')

    elif opts['--filter'] == "false":
        res.to_csv(opts['<output>'], index=False, sep='\t')
        print(res)

    else:
        print("Sorry, this is an invalid parameter ")





