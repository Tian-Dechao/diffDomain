"""
Test if one domain is significantly rewried in another biological condition
Usage: 
    scriptname dvsd one <chr> <start> <end> <hic0> <hic1> [options]
    scriptname dvsd multiple <hic0> <hic1> <bed> [options]

Options:
    --ofile filepath for output file  [default: stdout]
    --oprefix prefix for output files 
    --oprefixFig prefix for output figures
    --sep deliminator for hicfile  [default: \t]
    --hicnorm hic matrix normalization method  [default: KR]
    --chrn chromosome number  [default: ALL]
    --reso resolution for hicfile  [default: 100000]
    --ncore number of parallel process   [default: 10]

"""
import os
import sys
import pandas as pd
import numpy as np
from docopt import docopt
from  multiprocessing import Pool
#from matrix import comp2domins_by_twtest, loadtads
from domains import comp2domins_by_twtest, loadtads

opts = __doc__
scriptname = os.path.basename(__file__)
opts = opts.replace("scriptname", scriptname)
opts = docopt(doc=opts)
# save opts to output for every output
opts_df = pd.DataFrame.from_dict(opts, orient='index')
opts_df.index = '#' + opts_df.index

if(opts['dvsd']):
    if opts['one']:
        result= comp2domins_by_twtest(chrn=opts['<chr>'], start=int(opts['<start>']), end=int(opts['<end>']), reso=int(opts['--reso']), hicnorm=opts['--hicnorm'], fhic0=opts['<hic0>'], fhic1=opts['<hic1>'])
        print result


    if opts['multiple']:
        reso = int(opts['--reso'])
        tadb = loadtads(opts['<bed>'], sep=opts['--sep'], chrnum=opts['--chrn'], min_nbin=10, reso=reso)
        # Parallel implementation
        def comp2domins_by_twtest_parallel(i):
            tmp_res = comp2domins_by_twtest(chrn=tadb.iloc[i, 0], start=tadb.iloc[i, 1], end=tadb.iloc[i, 2], reso=reso, hicnorm=opts['--hicnorm'], fhic0=opts['<hic0>'], fhic1=opts['<hic1>'])
            return tmp_res

        P = Pool(int(opts['--ncore']))
        result = P.map(comp2domins_by_twtest_parallel, range(tadb.shape[0]))
        P.terminate()
        result = pd.DataFrame(result)

        # first the options; then the result
        opts_df.to_csv(opts['--ofile'], sep='\t', index=True, header=False)
        result.to_csv(opts['--ofile'], sep='\t', index=False, header=False, mode='a')

